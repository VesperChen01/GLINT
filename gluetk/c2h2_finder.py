# -*- coding: utf-8 -*-
"""
C2H2 Zinc Finger 发现模块
=========================
两层筛选策略：
1. 序列初筛（高召回）：HMM 域定位 + 弱规则打标签
2. 结构二筛（高精度）：全局 fold check + 局部 turn 对齐

设计原则：
- 第一层只负责"把对齐候选位置找全 + 给候选排个队"，不负责"判死刑"
- 弱规则用于优先级排序/分桶，不做 hard filter
"""
from __future__ import print_function
import os
import re
import csv
import tempfile
from collections import defaultdict
from typing import Optional, List, Tuple, Dict, Any, NamedTuple
from dataclasses import dataclass, field

# ============================================================================
# 数据结构
# ============================================================================

@dataclass
class C2H2Domain:
    """单个 C2H2 锌指域"""
    protein_id: str
    chain: str
    domain_start: int
    domain_end: int
    sequence: str
    finger_index: int = 0
    
    # HMM 相关
    hmm_score: float = 0.0
    hmm_evalue: float = 1.0
    detection_method: str = "regex"  # "hmm" or "regex"
    
    # 弱规则标签（用于优先级排序，不做 hard filter）
    has_turn_gly: bool = False          # turn 区域是否有关键 Gly
    num_gly_in_turn: int = 0            # turn 窗口内 G 数量
    spacer_cc_length: int = 0           # C..C 间隔长度
    spacer_hh_length: int = 0           # H..H 间隔长度
    spacer_valid: bool = True           # 间隔是否典型
    is_low_complexity: bool = False     # polyQ 等低复杂度
    
    # 结构二筛结果
    global_rmsd: Optional[float] = None
    global_tm_score: Optional[float] = None
    turn_rmsd: Optional[float] = None
    turn_tm_score: Optional[float] = None
    plddt_turn: Optional[float] = None  # AF2/ESM 结构质量
    
    # 综合评分
    priority_score: float = 0.0
    status: str = "candidate"  # candidate / pass / fail
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "Chain": self.chain,
            "Sequence": self.sequence,
            "Domain_Start": self.domain_start,
            "Domain_End": self.domain_end,
            "Finger_Index": self.finger_index,
            "HMM_Score": f"{self.hmm_score:.2f}" if self.hmm_score else "N/A",
            "Detection_Method": self.detection_method,
            "Turn_RMSD": f"{self.turn_rmsd:.2f}" if self.turn_rmsd is not None else "N/A",
            "Has_Turn_Gly": "Yes" if self.has_turn_gly else "No",
            "Spacer_Valid": "Yes" if self.spacer_valid else "No",
            "LowComplexity": "Yes" if self.is_low_complexity else "No",
            "Priority_Score": f"{self.priority_score:.2f}",
            "Status": self.status,
        }


# ============================================================================
# 内置模板定义
# ============================================================================

# C2H2 锌指典型模板（用于结构对齐）
BUILTIN_C2H2_TEMPLATES = {
    "Egr1 (ZF1, 1ZAA)": {
        "pdb": "1ZAA",
        "selection_fmt": "{obj} and chain A and resi 338-359 and name CA",
        "turn_resi": "346-352",  # β-turn 区域
        "description": "Egr1/Zif268 第一个锌指，经典 C2H2 结构",
    },
    "ZNF10 (ZF1, 2YT7)": {
        "pdb": "2YT7",
        "selection_fmt": "{obj} and chain A and resi 1-23 and name CA",
        "turn_resi": "8-14",
        "description": "ZNF10 第一个锌指",
    },
}

# 理想化 C2H2 β-turn 模板坐标（8×Cα）
# 基于 Egr1 ZF1 的 β-turn 区域
def _ideal_c2h2_turn_template():
    """返回理想化的 C2H2 β-turn 区域 Cα 坐标（7 个残基）"""
    # 典型 C2H2 turn 区域：C-X(2-4)-C 后的 ~7 残基形成 β-hairpin turn
    # 坐标基于 1ZAA Egr1 ZF1 归一化
    return [
        (0.0, 0.0, 0.0),     # turn 起始
        (3.8, 0.5, 0.2),
        (6.5, 2.8, 0.1),     # turn 顶点
        (5.2, 5.5, -0.3),
        (2.0, 6.2, 0.0),
        (-0.5, 4.0, 0.2),
        (-1.0, 1.0, 0.0),    # turn 结束
    ]


# ============================================================================
# 氨基酸工具
# ============================================================================

AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'MSE': 'M', 'SEC': 'U', 'PYL': 'O',
}

def _aa_3to1(resn: str) -> str:
    return AA_3TO1.get(resn.strip().upper(), 'X')

def _parse_resi(resi_str: str) -> Tuple[int, str]:
    """解析残基编号（支持插入码如 100A）"""
    if resi_str is None:
        return (0, "")
    s = str(resi_str).strip()
    m = re.match(r'(-?\d+)\s*([A-Za-z]?)', s)
    if m:
        return (int(m.group(1)), m.group(2) or "")
    try:
        return (int(s), "")
    except Exception:
        return (0, s or "")


# ============================================================================
# 序列初筛层
# ============================================================================

# C2H2 锌指的典型 pattern
# 标准格式：C-X(2-4)-C-X(3)-φ-X(5)-L-X(2)-H-X(3-5)-H
# 其中 φ = F/Y, L = L/F/Y（疏水残基）
C2H2_REGEX_PATTERN = re.compile(
    r'C'                    # 第一个 Cys
    r'[A-Z]{2,4}'          # 2-4 个任意氨基酸
    r'C'                    # 第二个 Cys
    r'[A-Z]{2,5}'          # 2-5 个氨基酸（包含 turn 区域）
    r'[FYWLIV]'            # 疏水残基（通常在 α-helix 起始）
    r'[A-Z]{3,8}'          # 3-8 个氨基酸
    r'[FYWLIV]?'           # 可选疏水残基
    r'[A-Z]{0,3}'          # 0-3 个氨基酸
    r'H'                    # 第一个 His
    r'[A-Z]{3,5}'          # 3-5 个氨基酸
    r'H'                    # 第二个 His
)

# 更宽松的 pattern（提高召回）
C2H2_REGEX_LOOSE = re.compile(
    r'C'                    # 第一个 Cys
    r'[A-Z]{1,6}'          # 1-6 个任意氨基酸
    r'C'                    # 第二个 Cys
    r'[A-Z]{8,20}'         # turn + α-helix 区域
    r'H'                    # 第一个 His
    r'[A-Z]{2,6}'          # 2-6 个氨基酸
    r'H'                    # 第二个 His
)


def _detect_low_complexity(seq: str, window: int = 10, threshold: float = 0.5) -> bool:
    """检测低复杂度区域（如 polyQ, polyN, polyS）"""
    if len(seq) < window:
        return False
    for i in range(len(seq) - window + 1):
        w = seq[i:i + window]
        # 单氨基酸重复
        for aa in "QNSEAGP":
            if w.count(aa) >= int(window * threshold):
                return True
        # 双氨基酸重复（如 GS, AG）
        for di in ["GS", "SG", "AG", "GA", "PS", "SP"]:
            count = sum(1 for j in range(0, len(w) - 1, 2) if w[j:j+2] == di)
            if count >= int(window * threshold / 2):
                return True
    return False


def _analyze_c2h2_soft_rules(seq: str, cc_start: int, cc_end: int, hh_start: int, hh_end: int) -> Dict[str, Any]:
    """
    弱规则分析（用于优先级排序，不做 hard filter）
    
    返回特征字典，每个特征都是"便宜"的计算
    """
    result = {
        "has_turn_gly": False,
        "num_gly_in_turn": 0,
        "spacer_cc_length": cc_end - cc_start - 1,
        "spacer_hh_length": hh_end - hh_start - 1,
        "spacer_valid": True,
        "is_low_complexity": False,
    }
    
    # 1. C..C 间隔检查（典型 2-4）
    cc_len = result["spacer_cc_length"]
    if cc_len < 2 or cc_len > 4:
        result["spacer_valid"] = False
    
    # 2. H..H 间隔检查（典型 3-5）
    hh_len = result["spacer_hh_length"]
    if hh_len < 3 or hh_len > 5:
        result["spacer_valid"] = False
    
    # 3. Turn 区域 Gly 检查（C2 和 H1 之间）
    # Turn 区域大约在第二个 C 后 2-6 位
    turn_region = seq[cc_end + 2:cc_end + 7] if cc_end + 7 <= len(seq) else seq[cc_end + 2:]
    result["num_gly_in_turn"] = turn_region.count('G')
    result["has_turn_gly"] = result["num_gly_in_turn"] >= 1
    
    # 4. 低复杂度检查
    result["is_low_complexity"] = _detect_low_complexity(seq)
    
    return result


def _sequence_search_regex(sequence: str, chain: str = "", protein_id: str = "") -> List[C2H2Domain]:
    """
    使用 regex 进行序列搜索（回退方案）
    
    返回所有匹配的 C2H2 域候选
    """
    domains = []
    finger_idx = 0
    
    # 先用宽松 pattern 找所有可能
    for m in C2H2_REGEX_LOOSE.finditer(sequence):
        start, end = m.start(), m.end()
        domain_seq = sequence[start:end]
        
        # 定位 C..C 和 H..H
        c_positions = [i for i, c in enumerate(domain_seq) if c == 'C']
        h_positions = [i for i, c in enumerate(domain_seq) if c == 'H']
        
        if len(c_positions) < 2 or len(h_positions) < 2:
            continue
        
        cc_start, cc_end = c_positions[0], c_positions[1]
        hh_start, hh_end = h_positions[-2], h_positions[-1]
        
        # 弱规则分析
        soft_rules = _analyze_c2h2_soft_rules(domain_seq, cc_start, cc_end, hh_start, hh_end)
        
        finger_idx += 1
        domain = C2H2Domain(
            protein_id=protein_id or "unknown",
            chain=chain,
            domain_start=start + 1,  # 1-based
            domain_end=end,
            sequence=domain_seq,
            finger_index=finger_idx,
            detection_method="regex",
            **soft_rules,
        )
        domains.append(domain)
    
    return domains


def _sequence_search_hmm(sequence: str, chain: str = "", protein_id: str = "",
                          hmm_profile_path: Optional[str] = None) -> List[C2H2Domain]:
    """
    使用 HMM（pyhmmer）进行序列搜索（首选方案）
    
    需要 pyhmmer 和 Pfam zf-C2H2.hmm
    """
    try:
        import pyhmmer
        from pyhmmer import easel, plan7
    except ImportError:
        print("[C2H2] pyhmmer 未安装，回退到 regex 模式")
        return _sequence_search_regex(sequence, chain, protein_id)
    
    # 查找 HMM profile
    if hmm_profile_path is None:
        # 尝试在插件目录查找
        here = os.path.dirname(os.path.abspath(__file__))
        possible_paths = [
            os.path.join(here, "data", "zf-C2H2.hmm"),
            os.path.join(here, "hmm", "zf-C2H2.hmm"),
            os.path.join(here, "PF00096.hmm"),  # Pfam ID
        ]
        for p in possible_paths:
            if os.path.exists(p):
                hmm_profile_path = p
                break
    
    if hmm_profile_path is None or not os.path.exists(hmm_profile_path):
        print("[C2H2] HMM profile 未找到，回退到 regex 模式")
        print("       建议下载 Pfam zf-C2H2 (PF00096) HMM 并放置于 gluetk/data/ 目录")
        return _sequence_search_regex(sequence, chain, protein_id)
    
    domains = []
    
    try:
        # 读取 HMM
        with pyhmmer.plan7.HMMFile(hmm_profile_path) as hmm_file:
            hmm = hmm_file.read()
        
        # 创建序列对象
        alphabet = easel.Alphabet.amino()
        seq = easel.TextSequence(name=b"query", sequence=sequence)
        digital_seq = seq.digitize(alphabet)
        
        # 运行 hmmsearch
        pipeline = plan7.Pipeline(alphabet)
        hits = pipeline.search_hmm(hmm, [digital_seq])
        
        finger_idx = 0
        for hit in hits:
            for domain in hit.domains:
                finger_idx += 1
                start = domain.env_from  # envelope coordinates
                end = domain.env_to
                domain_seq = sequence[start - 1:end]
                
                # 定位 C..C 和 H..H
                c_positions = [i for i, c in enumerate(domain_seq) if c == 'C']
                h_positions = [i for i, c in enumerate(domain_seq) if c == 'H']
                
                if len(c_positions) >= 2 and len(h_positions) >= 2:
                    cc_start, cc_end = c_positions[0], c_positions[1]
                    hh_start, hh_end = h_positions[-2], h_positions[-1]
                    soft_rules = _analyze_c2h2_soft_rules(domain_seq, cc_start, cc_end, hh_start, hh_end)
                else:
                    soft_rules = {
                        "has_turn_gly": False,
                        "num_gly_in_turn": 0,
                        "spacer_cc_length": 0,
                        "spacer_hh_length": 0,
                        "spacer_valid": False,
                        "is_low_complexity": _detect_low_complexity(domain_seq),
                    }
                
                c2h2 = C2H2Domain(
                    protein_id=protein_id or "unknown",
                    chain=chain,
                    domain_start=start,
                    domain_end=end,
                    sequence=domain_seq,
                    finger_index=finger_idx,
                    hmm_score=domain.score,
                    hmm_evalue=domain.c_evalue,
                    detection_method="hmm",
                    **soft_rules,
                )
                domains.append(c2h2)
        
        print(f"[C2H2] HMM 检测到 {len(domains)} 个锌指域")
        
    except Exception as e:
        print(f"[C2H2] HMM 搜索失败: {e}，回退到 regex")
        return _sequence_search_regex(sequence, chain, protein_id)
    
    # 如果 HMM 没找到任何结果，也用 regex 兜底
    if not domains:
        print("[C2H2] HMM 无结果，使用 regex 兜底")
        return _sequence_search_regex(sequence, chain, protein_id)
    
    return domains


# ============================================================================
# 结构二筛层
# ============================================================================

def _kabsch_rmsd(P, Q):
    """Kabsch 算法计算 RMSD"""
    import numpy as np
    P = np.asarray(P, float)
    Q = np.asarray(Q, float)
    Pc = P.mean(0)
    Qc = Q.mean(0)
    P0 = P - Pc
    Q0 = Q - Qc
    C = P0.T @ Q0
    V, S, Wt = np.linalg.svd(C)
    if np.linalg.det(V @ Wt) < 0.0:
        V[:, -1] *= -1.0
    R = V @ Wt
    P_aln = (R @ P0.T).T + Qc
    diff = P_aln - Q
    return float((diff ** 2).sum() / len(P)) ** 0.5


def _get_ca_coords_for_region(obj_name: str, chain: str, start: int, end: int) -> List[Tuple[float, float, float]]:
    """从 PyMOL 对象获取指定区域的 Cα 坐标"""
    try:
        from pymol import cmd
    except ImportError:
        return []
    
    sel = f"{obj_name} and chain {chain} and resi {start}-{end} and name CA"
    model = cmd.get_model(sel)
    
    coords = []
    for a in model.atom:
        alt = (a.alt or '').strip().upper()
        if alt not in ("", "A"):
            continue
        coords.append((a.coord[0], a.coord[1], a.coord[2]))
    
    return coords


def _get_plddt_for_region(obj_name: str, chain: str, start: int, end: int) -> Optional[float]:
    """获取指定区域的平均 pLDDT（存储在 B-factor 中）"""
    try:
        from pymol import cmd
    except ImportError:
        return None
    
    sel = f"{obj_name} and chain {chain} and resi {start}-{end} and name CA"
    model = cmd.get_model(sel)
    
    b_factors = []
    for a in model.atom:
        b = getattr(a, 'b', None)
        if b is not None:
            b_factors.append(b)
    
    if not b_factors:
        return None
    
    avg_b = sum(b_factors) / len(b_factors)
    
    # pLDDT 通常在 0-100 范围；如果 B-factor 在这个范围，可能就是 pLDDT
    if 0 <= avg_b <= 100:
        return avg_b
    
    return None


def _structural_filter(domain: C2H2Domain, obj_name: str,
                       template_coords: Optional[List[Tuple[float, float, float]]] = None,
                       global_rmsd_cutoff: float = 3.0,
                       turn_rmsd_cutoff: float = 2.0,
                       plddt_threshold: float = 70.0) -> C2H2Domain:
    """
    结构二筛：全局 fold check + 局部 turn 对齐
    
    返回更新后的 domain（填充结构相关字段）
    """
    if template_coords is None:
        template_coords = _ideal_c2h2_turn_template()
    
    # 1. 获取域的 Cα 坐标
    domain_coords = _get_ca_coords_for_region(
        obj_name, domain.chain, domain.domain_start, domain.domain_end
    )
    
    if len(domain_coords) < 7:
        domain.status = "fail"
        return domain
    
    # 2. 全局 fold check（取中间部分，避免端部不稳定）
    try:
        # 取域的前 7 个 Cα 做粗略比对
        global_coords = domain_coords[:min(7, len(domain_coords))]
        if len(global_coords) >= len(template_coords):
            global_coords = global_coords[:len(template_coords)]
            domain.global_rmsd = _kabsch_rmsd(global_coords, template_coords)
        
        # 如果全局 RMSD 太大，直接标记为 fail
        if domain.global_rmsd is not None and domain.global_rmsd > global_rmsd_cutoff * 2:
            domain.status = "fail"
            return domain
        
    except Exception as e:
        print(f"[C2H2] 全局 RMSD 计算失败: {e}")
    
    # 3. 局部 turn 对齐（更敏感）
    # Turn 区域：C2 后 2-8 位（约 7 个残基）
    try:
        # 估算 turn 起始位置（第二个 C 后约 2 位）
        seq = domain.sequence
        c_positions = [i for i, c in enumerate(seq) if c == 'C']
        if len(c_positions) >= 2:
            turn_start_offset = c_positions[1] + 2
            turn_end_offset = turn_start_offset + 7
            
            # 转换为绝对残基编号
            turn_start = domain.domain_start + turn_start_offset
            turn_end = min(domain.domain_start + turn_end_offset, domain.domain_end)
            
            turn_coords = _get_ca_coords_for_region(obj_name, domain.chain, turn_start, turn_end)
            
            if len(turn_coords) >= 5:
                # 截取与模板相同长度
                n = min(len(turn_coords), len(template_coords))
                domain.turn_rmsd = _kabsch_rmsd(turn_coords[:n], template_coords[:n])
    except Exception as e:
        print(f"[C2H2] Turn RMSD 计算失败: {e}")
    
    # 4. pLDDT 检查（若有 AF2/ESM 结构）
    plddt = _get_plddt_for_region(obj_name, domain.chain, domain.domain_start, domain.domain_end)
    if plddt is not None:
        domain.plddt_turn = plddt
    
    # 5. 综合状态判定
    if domain.turn_rmsd is not None and domain.turn_rmsd <= turn_rmsd_cutoff:
        domain.status = "pass"
    elif domain.global_rmsd is not None and domain.global_rmsd <= global_rmsd_cutoff:
        domain.status = "pass"
    else:
        domain.status = "candidate"  # 保持候选，不直接 fail
    
    return domain


# ============================================================================
# 优先级评分
# ============================================================================

def _calculate_priority_score(domain: C2H2Domain) -> float:
    """
    计算综合优先级评分
    
    评分原则：
    - HMM 检测 > regex 检测
    - 有 turn Gly > 无
    - 间隔典型 > 间隔异常
    - 低复杂度 降权
    - 结构 RMSD 好 > 差
    - pLDDT 高 > 低
    """
    score = 50.0  # 基础分
    
    # 检测方法
    if domain.detection_method == "hmm":
        score += 20.0
        # HMM 分数加成（通常 10-30）
        if domain.hmm_score > 0:
            score += min(domain.hmm_score, 30)
    
    # 弱规则加分
    if domain.has_turn_gly:
        score += 10.0
    score += domain.num_gly_in_turn * 3.0
    
    if domain.spacer_valid:
        score += 10.0
    else:
        score -= 10.0
    
    # 低复杂度惩罚
    if domain.is_low_complexity:
        score -= 20.0
    
    # 结构评分
    if domain.turn_rmsd is not None:
        if domain.turn_rmsd < 1.5:
            score += 30.0
        elif domain.turn_rmsd < 2.5:
            score += 15.0
        elif domain.turn_rmsd > 4.0:
            score -= 15.0
    
    if domain.plddt_turn is not None:
        if domain.plddt_turn >= 90:
            score += 15.0
        elif domain.plddt_turn >= 70:
            score += 5.0
        elif domain.plddt_turn < 50:
            score -= 15.0
    
    return max(0.0, score)


# ============================================================================
# 主函数
# ============================================================================

def find_c2h2_domains(obj_name: Optional[str] = None,
                       pdb_file: Optional[str] = None,
                       hmm_profile: Optional[str] = None,
                       turn_rmsd_cutoff: float = 2.0,
                       global_rmsd_cutoff: float = 3.5,
                       require_turn_gly: bool = False,
                       skip_low_complexity: bool = False,
                       out_csv: Optional[str] = None,
                       auto_highlight: bool = True,
                       topk_debug: int = 20) -> List[C2H2Domain]:
    """
    C2H2 锌指蛋白识别主函数
    
    两层筛选策略：
    1. 序列初筛（高召回）：HMM + 弱规则打标签
    2. 结构二筛（高精度）：全局 fold + 局部 turn 对齐
    
    参数:
        obj_name: PyMOL 对象名
        pdb_file: PDB 文件路径
        hmm_profile: HMM profile 路径（默认查找 Pfam zf-C2H2）
        turn_rmsd_cutoff: 局部 turn RMSD 阈值（默认 2.0 Å）
        global_rmsd_cutoff: 全局 RMSD 阈值（默认 3.5 Å）
        require_turn_gly: 是否要求 turn 区域有 Gly（False = 不做 hard filter）
        skip_low_complexity: 是否跳过低复杂度区域（False = 保留但降权）
        out_csv: 输出 CSV 路径
        auto_highlight: 是否自动高亮
        topk_debug: 显示前 N 个结果（调试）
    
    返回:
        List[C2H2Domain]: 所有候选域（按优先级排序）
    """
    try:
        from pymol import cmd
    except ImportError:
        print("[C2H2] PyMOL 未加载，无法进行结构分析")
        return []
    
    # 1. 载入对象
    tmp_obj = None
    if pdb_file:
        tmp_obj = "_c2h2_tmp_obj"
        cmd.load(pdb_file, tmp_obj, quiet=1)
        obj = tmp_obj
    else:
        try:
            objs = cmd.get_names("objects")
        except AttributeError:
            objs = cmd.get_object_list() if hasattr(cmd, "get_object_list") else []
        obj = obj_name or (objs[0] if objs else None)
    
    if not obj:
        print("[C2H2] 无可用对象；请加载结构或提供 pdb_file")
        return []
    
    print(f"[C2H2] 分析对象: {obj}")
    
    # 2. 提取序列（按链）
    all_domains = []
    model = cmd.get_model(f"{obj} and polymer.protein and name CA")
    
    # 按链组织残基
    chains = defaultdict(list)
    for a in model.atom:
        chain = (a.chain or '').strip() or "_"
        resn = (a.resn or '').strip().upper()
        resi = (a.resi or '').strip()
        chains[chain].append((resi, resn))
    
    # 3. 对每条链进行序列搜索
    for chain, residues in chains.items():
        # 去重并排序
        seen = set()
        unique_residues = []
        for resi, resn in residues:
            if resi not in seen:
                seen.add(resi)
                unique_residues.append((resi, resn))
        
        # 按残基编号排序
        unique_residues.sort(key=lambda x: _parse_resi(x[0]))
        
        # 构建序列
        sequence = ''.join(_aa_3to1(resn) for _, resn in unique_residues)
        
        if len(sequence) < 20:  # C2H2 至少需要 ~23 aa
            continue
        
        print(f"[C2H2] Chain {chain}: {len(sequence)} aa")
        
        # 序列搜索（优先 HMM，回退 regex）
        domains = _sequence_search_hmm(sequence, chain, obj, hmm_profile)
        
        # 调整域的绝对位置（相对于整条链）
        resi_map = {i: resi for i, (resi, _) in enumerate(unique_residues)}
        for d in domains:
            # 将序列索引转换为残基编号
            if d.domain_start - 1 in resi_map:
                d.domain_start = int(resi_map[d.domain_start - 1])
            if d.domain_end - 1 in resi_map:
                d.domain_end = int(resi_map[d.domain_end - 1])
        
        all_domains.extend(domains)
    
    print(f"[C2H2] 序列初筛: {len(all_domains)} 个候选域")
    
    # 4. 结构二筛
    template_coords = _ideal_c2h2_turn_template()
    
    for i, domain in enumerate(all_domains):
        all_domains[i] = _structural_filter(
            domain, obj, template_coords,
            global_rmsd_cutoff=global_rmsd_cutoff,
            turn_rmsd_cutoff=turn_rmsd_cutoff,
        )
    
    # 5. 计算优先级评分
    for domain in all_domains:
        domain.priority_score = _calculate_priority_score(domain)
    
    # 6. 应用筛选（可选）
    filtered_domains = all_domains
    if require_turn_gly:
        filtered_domains = [d for d in filtered_domains if d.has_turn_gly]
    if skip_low_complexity:
        filtered_domains = [d for d in filtered_domains if not d.is_low_complexity]
    
    # 7. 按优先级排序
    filtered_domains.sort(key=lambda x: x.priority_score, reverse=True)
    
    # 8. 调试输出
    print(f"\n[C2H2] 结果汇总: {len(filtered_domains)} 个域")
    print("=" * 80)
    n_show = min(topk_debug, len(filtered_domains))
    for i, d in enumerate(filtered_domains[:n_show]):
        status_icon = "✅" if d.status == "pass" else ("⚠️" if d.status == "candidate" else "❌")
        turn_rmsd_str = f"{d.turn_rmsd:.2f}" if d.turn_rmsd else "N/A"
        print(f"  #{i+1:02d} {status_icon} Chain {d.chain} {d.domain_start}-{d.domain_end} "
              f"[{d.detection_method}] seq={d.sequence[:15]}... "
              f"turn_RMSD={turn_rmsd_str} score={d.priority_score:.1f}")
    print("=" * 80)
    
    # 9. 写 CSV
    if out_csv is None:
        fd, out_csv = tempfile.mkstemp(suffix="_c2h2.csv")
        os.close(fd)
    
    with open(out_csv, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "Chain", "Sequence", "Domain_Start", "Domain_End", "Finger_Index",
            "HMM_Score", "Detection_Method", "Turn_RMSD", "Has_Turn_Gly",
            "Spacer_Valid", "LowComplexity", "Priority_Score", "Status"
        ])
        writer.writeheader()
        for d in filtered_domains:
            writer.writerow(d.to_dict())
    
    print(f"[C2H2] 结果已保存: {out_csv}")
    
    # 10. 自动高亮（可选）
    if auto_highlight and filtered_domains:
        try:
            _highlight_c2h2_domains(obj, filtered_domains)
        except Exception as e:
            print(f"[C2H2] 高亮失败: {e}")
    
    # 清理临时对象
    if tmp_obj:
        try:
            cmd.delete(tmp_obj)
        except:
            pass
    
    return filtered_domains


def _highlight_c2h2_domains(obj_name: str, domains: List[C2H2Domain], color: str = "magenta"):
    """高亮 C2H2 域"""
    try:
        from pymol import cmd
    except ImportError:
        return
    
    for i, d in enumerate(domains[:20]):  # 最多高亮 20 个
        sel_name = f"c2h2_zf_{i+1}"
        sel_expr = f"{obj_name} and chain {d.chain} and resi {d.domain_start}-{d.domain_end}"
        
        cmd.select(sel_name, sel_expr)
        cmd.show("cartoon", sel_name)
        
        # 根据状态选择颜色
        if d.status == "pass":
            cmd.color("tv_green", sel_name)
        elif d.status == "fail":
            cmd.color("tv_red", sel_name)
        else:
            cmd.color(color, sel_name)
        
        # 显示 Cys 和 His 为 sticks
        cmd.show("sticks", f"{sel_name} and resn CYS+HIS")
        cmd.color("yellow", f"{sel_name} and resn CYS and name SG")
        cmd.color("blue", f"{sel_name} and resn HIS and name ND1+NE2")
    
    print(f"[C2H2] 已高亮 {min(len(domains), 20)} 个锌指域")


# ============================================================================
# PyMOL 命令注册
# ============================================================================

try:
    from pymol import cmd
    cmd.extend("find_c2h2_domains", find_c2h2_domains)
except ImportError:
    pass
