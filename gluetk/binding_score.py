# -*- coding: utf-8 -*-
"""
binding_score.py
分子胶/PROTAC三元复合物结合能评分模块

基于经验评分函数,快速估算结合亲和力
- 支持单蛋白-配体二元复合物
- 支持分子胶/PROTAC三元复合物(蛋白-配体-蛋白)
- 零外部依赖,跨平台,适合初步筛选

理论基础:
- 氢键贡献: -0.5 to -3.0 kcal/mol
- 盐桥: -1.0 to -5.0 kcal/mol  
- 疏水接触: -0.15 kcal/mol per contact
- 去溶剂化惩罚: +2-5 kcal/mol
- 协同效应: ΔΔG_cooperativity

参考文献:
- Böhm, H.-J. (1994). J. Comput.-Aided Mol. Des., 8, 243-256.
- Tripos SYBYL scoring functions
"""

from __future__ import print_function
import math
from collections import defaultdict

# ========== 评分参数 ==========
SCORE_PARAMS = {
    # 氢键能量 (kcal/mol) - 按距离分段
    "hbond": {
        "excellent": (-2.5, 2.5),   # ≤2.5Å: -2.5 kcal/mol
        "good": (-2.0, 2.8),        # 2.5-2.8Å: -2.0
        "moderate": (-1.2, 3.2),    # 2.8-3.2Å: -1.2
        "weak": (-0.5, 3.5)         # 3.2-3.5Å: -0.5
    },
    
    # 盐桥能量 (kcal/mol)
    "ionic": {
        "strong": (-4.0, 3.5),      # ≤3.5Å: -4.0 kcal/mol
        "moderate": (-2.5, 4.0),    # 3.5-4.0Å: -2.5
        "weak": (-1.0, 4.5)         # 4.0-4.5Å: -1.0
    },
    
    # 疏水接触
    "hydrophobic": -0.15,  # kcal/mol per contact
    
    # π-π堆积
    "pi_pi": -1.2,  # kcal/mol
    
    # π-阳离子
    "pi_cation": -1.5,  # kcal/mol
    
    # 金属配位
    "metal": -3.5,  # kcal/mol
    
    # 卤素键
    "halogen": -0.8,  # kcal/mol
    
    # 去溶剂化惩罚(基于接触数估算)
    "desolvation_penalty_per_contact": 0.05,  # kcal/mol
    "desolvation_base": 2.0,  # kcal/mol baseline
    
    # 三元复合物协同系数
    "cooperativity": {
        "balanced_threshold": 0.3,  # 相互作用平衡度阈值
        "bonus_max": -3.0,  # 最大协同奖励 (kcal/mol)
        "penalty_imbalance": 1.5  # 不平衡惩罚系数
    }
}


def _score_distance_based(dist, tiers):
    """
    根据距离分段计算能量贡献
    
    参数:
        dist: 距离(Å)
        tiers: 分段字典 {label: (energy, max_dist), ...}
    
    返回:
        energy: kcal/mol
    """
    for label in ["excellent", "good", "moderate", "weak"]:
        if label in tiers:
            energy, max_dist = tiers[label]
            if dist <= max_dist:
                return energy
    return 0.0


def calculate_binary_score(interactions_result):
    """
    计算蛋白-配体二元复合物结合能
    
    参数:
        interactions_result: analyze_protein_ligand_interactions() 返回的字典
    
    返回:
        dict: {
            'total': float,  # 总估算结合能 (kcal/mol)
            'components': {
                'hbond': float,
                'ionic': float,
                'hydrophobic': float,
                'pi_pi': float,
                'pi_cation': float,
                'metal': float,
                'halogen': float,
                'desolvation': float  # 正值(惩罚项)
            },
            'interaction_counts': {...}
        }
    """
    if not interactions_result or 'interactions' not in interactions_result:
        return None
    
    components = {
        'hbond': 0.0,
        'ionic': 0.0,
        'hydrophobic': 0.0,
        'pi_pi': 0.0,
        'pi_cation': 0.0,
        'metal': 0.0,
        'halogen': 0.0,
        'desolvation': 0.0
    }
    
    counts = defaultdict(int)
    
    # 解析相互作用列表
    for inter in interactions_result['interactions']:
        inter_type = inter.get('Interaction', '').lower()
        dist = float(inter.get('Distance', 999))
        
        # 氢键
        if '氢键' in inter_type or 'hydrogen' in inter_type or 'hbond' in inter_type:
            components['hbond'] += _score_distance_based(dist, SCORE_PARAMS['hbond'])
            counts['hbond'] += 1
        
        # 盐桥
        elif '盐桥' in inter_type or 'salt' in inter_type or 'ionic' in inter_type:
            components['ionic'] += _score_distance_based(dist, SCORE_PARAMS['ionic'])
            counts['ionic'] += 1
        
        # 疏水
        elif '疏水' in inter_type or 'hydrophobic' in inter_type:
            components['hydrophobic'] += SCORE_PARAMS['hydrophobic']
            counts['hydrophobic'] += 1
        
        # π-π
        elif 'π-π' in inter_type or 'pi-pi' in inter_type or 'pi_pi' in inter_type:
            components['pi_pi'] += SCORE_PARAMS['pi_pi']
            counts['pi_pi'] += 1
        
        # π-阳离子
        elif 'π-阳离子' in inter_type or 'pi-cation' in inter_type or 'pi_cation' in inter_type:
            components['pi_cation'] += SCORE_PARAMS['pi_cation']
            counts['pi_cation'] += 1
        
        # 金属配位
        elif '金属' in inter_type or 'metal' in inter_type:
            components['metal'] += SCORE_PARAMS['metal']
            counts['metal'] += 1
        
        # 卤素键
        elif '卤素' in inter_type or 'halogen' in inter_type:
            components['halogen'] += SCORE_PARAMS['halogen']
            counts['halogen'] += 1
    
    # 去溶剂化惩罚(基于总接触数)
    total_contacts = sum(counts.values())
    components['desolvation'] = (SCORE_PARAMS['desolvation_base'] + 
                                  total_contacts * SCORE_PARAMS['desolvation_penalty_per_contact'])
    
    # 总分
    total_score = sum(components.values())
    
    return {
        'total': round(total_score, 2),
        'components': {k: round(v, 2) for k, v in components.items()},
        'interaction_counts': dict(counts),
        'contact_count': total_contacts
    }


def calculate_ternary_score(ternary_result, ppi_result=None, neo_epitope_result=None):
    """
    计算分子胶/PROTAC三元复合物结合能(含协同效应)
    
    ✨ 新增分子胶特异分析:
    - ppi_result: 蛋白-蛋白界面分析结果 (来自 ppi_analyzer.analyze_protein_protein_interface)
    - neo_epitope_result: Neo-表位识别结果 (来自 ppi_analyzer.identify_neo_epitope)
    
    参数:
        ternary_result: analyze_ternary_complex() 返回的字典,包含:
            - protein1_result: 蛋白1-配体相互作用
            - protein2_result: 蛋白2-配体相互作用
            - bridging_analysis: 桥接分析
        ppi_result: (可选) PPI分析结果,用于分子胶判定
        neo_epitope_result: (可选) Neo-表位分析结果
    
    返回:
        dict: {
            'total': float,  # 总结合能 (kcal/mol)
            'protein1_score': {...},  # 蛋白1评分
            'protein2_score': {...},  # 蛋白2评分
            'cooperativity': float,  # 协同效应 (kcal/mol)
            'cooperativity_breakdown': {...},  # 协同效应分解
            'balance_factor': float,  # 相互作用平衡度 [0-1]
            'mechanism': str,  # "Molecular Glue" 或 "PROTAC" 或 "Unknown"
            'recommendation': str  # 设计建议
        }
    """
    # 获取两个蛋白的相互作用结果
    p1_result = ternary_result.get('protein1_result')
    p2_result = ternary_result.get('protein2_result')
    
    if not p1_result or not p2_result:
        print("[calculate_ternary_score] ⚠️ Missing protein interaction data")
        return None
    
    # 分别计算两个蛋白的评分
    p1_score = calculate_binary_score(p1_result)
    p2_score = calculate_binary_score(p2_result)
    
    if not p1_score or not p2_score:
        return None
    
    # 计算相互作用平衡度
    # balance_factor: 0 = 完全不平衡, 1 = 完全平衡
    p1_contacts = p1_score['contact_count']
    p2_contacts = p2_score['contact_count']
    total_contacts = p1_contacts + p2_contacts
    
    if total_contacts == 0:
        balance_factor = 0.0
    else:
        # 使用Shannon entropy衡量平衡度
        p1_frac = p1_contacts / total_contacts
        p2_frac = p2_contacts / total_contacts
        
        if p1_frac == 0 or p2_frac == 0:
            balance_factor = 0.0
        else:
            entropy = -(p1_frac * math.log2(p1_frac) + p2_frac * math.log2(p2_frac))
            balance_factor = entropy / 1.0  # normalize to [0, 1]
    
    # ========== 基础协同效应 (PROTAC模式) ==========
    cooperativity_base = 0.0
    coop_params = SCORE_PARAMS['cooperativity']
    
    if balance_factor >= coop_params['balanced_threshold']:
        # 平衡的相互作用 → 协同奖励
        cooperativity_base = coop_params['bonus_max'] * balance_factor * min(1.0, total_contacts / 10.0)
    else:
        # 不平衡 → 惩罚
        cooperativity_base = coop_params['penalty_imbalance'] * (1.0 - balance_factor)
    
    # ========== 分子胶特异协同因子 ==========
    cooperativity_factors = {
        'base': cooperativity_base,
        'ppi_bonus': 0.0,           # 蛋白-蛋白接触奖励
        'neo_epitope_bonus': 0.0,   # Neo-表位奖励
        'interface_strength': 0.0,  # 界面强度奖励
        'bsa_bonus': 0.0            # 埋藏表面积奖励
    }
    
    mechanism = "Unknown"  # "Molecular Glue" / "PROTAC" / "Unknown"
    
    # 如果提供了PPI分析结果
    if ppi_result:
        interface_contacts = ppi_result.get('interface_contacts', 0)
        is_strong_interface = ppi_result.get('is_strong_interface', False)
        interface_strength = ppi_result.get('interface_strength', 0.0)
        bsa = ppi_result.get('bsa', 0.0)
        
        # PPI接触奖励 (强界面显著奖励)
        if interface_contacts >= 10:
            cooperativity_factors['ppi_bonus'] = -5.0  # 强PPI → 大奖励
            mechanism = "Molecular Glue"
        elif interface_contacts >= 5:
            cooperativity_factors['ppi_bonus'] = -2.5  # 中等PPI
        
        # 界面强度奖励
        if interface_strength > 6.0:
            cooperativity_factors['interface_strength'] = -2.0
        
        # BSA奖励 (BSA > 800 Ų 认为是强界面)
        if bsa and bsa > 800.0:
            cooperativity_factors['bsa_bonus'] = -3.0
            mechanism = "Molecular Glue"  # 高BSA强烈提示分子胶
        elif bsa and bsa > 400.0:
            cooperativity_factors['bsa_bonus'] = -1.5
    
    # 如果提供了Neo-表位分析结果
    if neo_epitope_result:
        neo_epitope_count = neo_epitope_result.get('neo_epitope_count', 0)
        is_molecular_glue = neo_epitope_result.get('is_molecular_glue', False)
        confidence = neo_epitope_result.get('confidence', 0.0)
        
        # Neo-表位奖励 (分子胶的核心特征)
        if is_molecular_glue:
            cooperativity_factors['neo_epitope_bonus'] = -4.0 * confidence  # 最高-4.0 kcal/mol
            mechanism = "Molecular Glue"  # Neo-表位是分子胶的决定性证据
        elif neo_epitope_count >= 2:
            cooperativity_factors['neo_epitope_bonus'] = -2.0 * confidence
    
    # 如果没有明显PPI/Neo-表位特征，判定为PROTAC
    if mechanism == "Unknown":
        if ppi_result and ppi_result.get('interface_contacts', 0) < 3:
            mechanism = "PROTAC"  # 弱/无PPI通常是PROTAC
    
    # 总协同效应
    cooperativity_total = sum(cooperativity_factors.values())
    
    # 总结合能 = 蛋白1 + 蛋白2 + 协同效应
    total_score = p1_score['total'] + p2_score['total'] + cooperativity_total
    
    # 生成设计建议
    recommendation = _generate_recommendation_v2(
        p1_score, p2_score, balance_factor, cooperativity_total,
        mechanism, ppi_result, neo_epitope_result
    )
    
    return {
        'total': round(total_score, 2),
        'protein1_score': p1_score,
        'protein2_score': p2_score,
        'cooperativity': round(cooperativity_total, 2),
        'cooperativity_breakdown': {k: round(v, 2) for k, v in cooperativity_factors.items()},
        'balance_factor': round(balance_factor, 3),
        'mechanism': mechanism,
        'recommendation': recommendation,
        'breakdown': {
            'protein1_total': p1_score['total'],
            'protein2_total': p2_score['total'],
            'cooperativity_bonus': round(cooperativity_total, 2)
        }
    }


def _generate_recommendation(p1_score, p2_score, balance_factor, cooperativity):
    """
    基于评分生成分子设计建议（旧版，保留向后兼容）
    """
    recommendations = []
    
    # 平衡度建议
    if balance_factor < 0.3:
        weak_side = "Protein 1" if p1_score['contact_count'] < p2_score['contact_count'] else "Protein 2"
        recommendations.append(f"⚠️ 不平衡设计: {weak_side} 接触太少,考虑增强该侧相互作用")
    elif balance_factor > 0.7:
        recommendations.append("✅ 良好的平衡设计")
    
    # 氢键建议
    p1_hb = p1_score['interaction_counts'].get('hbond', 0)
    p2_hb = p2_score['interaction_counts'].get('hbond', 0)
    
    if p1_hb < 2:
        recommendations.append("💡 Protein 1: 氢键数量偏少(<2),考虑引入氢键受/供体")
    if p2_hb < 2:
        recommendations.append("💡 Protein 2: 氢键数量偏少(<2),考虑引入氢键受/供体")
    
    # 协同效应
    if cooperativity < -1.0:
        recommendations.append(f"✨ 检测到显著协同效应({cooperativity:.1f} kcal/mol)")
    elif cooperativity > 1.0:
        recommendations.append(f"⚠️ 不平衡惩罚较大(+{cooperativity:.1f} kcal/mol)")
    
    # 总分评估
    total = p1_score['total'] + p2_score['total'] + cooperativity
    if total < -10.0:
        recommendations.append("🎯 预测结合能优秀,建议实验验证")
    elif total > -5.0:
        recommendations.append("⚠️ 预测结合能偏弱,建议优化设计")
    
    return " | ".join(recommendations) if recommendations else "无特殊建议"


def _generate_recommendation_v2(p1_score, p2_score, balance_factor, cooperativity,
                                mechanism, ppi_result, neo_epitope_result):
    """
    基于评分生成分子设计建议（分子胶特异版本）
    
    参数:
        mechanism: "Molecular Glue" / "PROTAC" / "Unknown"
        ppi_result: PPI分析结果
        neo_epitope_result: Neo-表位分析结果
    """
    recommendations = []
    
    # 机制判定
    if mechanism == "Molecular Glue":
        recommendations.append("✨ 检测到分子胶机制 (Molecular Glue)")
    elif mechanism == "PROTAC":
        recommendations.append("🔗 检测到PROTAC机制 (Linker-based)")
    else:
        recommendations.append("❓ 未能确定机制类型")
    
    # PPI相关建议
    if ppi_result:
        interface_contacts = ppi_result.get('interface_contacts', 0)
        bsa = ppi_result.get('bsa')
        
        if interface_contacts >= 10:
            recommendations.append(f"💪 强蛋白-蛋白界面 ({interface_contacts}个接触)")
        elif interface_contacts >= 5:
            recommendations.append(f"⚡ 中等蛋白-蛋白界面 ({interface_contacts}个接触)")
        elif interface_contacts < 3 and mechanism == "Unknown":
            recommendations.append("⚠️ 蛋白-蛋白接触弱,可能是PROTAC或需优化界面")
        
        if bsa:
            if bsa > 800.0:
                recommendations.append(f"🌟 高埋藏表面积 (BSA={bsa:.1f}Ų),界面稳定")
            elif bsa < 400.0:
                recommendations.append(f"⚠️ 低埋藏表面积 (BSA={bsa:.1f}Ų),考虑增强界面")
    
    # Neo-表位相关建议
    if neo_epitope_result:
        neo_count = neo_epitope_result.get('neo_epitope_count', 0)
        confidence = neo_epitope_result.get('confidence', 0.0)
        is_glue = neo_epitope_result.get('is_molecular_glue', False)
        
        if is_glue:
            recommendations.append(f"🎯 识别到{neo_count}个Neo-表位残基 (置信度={confidence:.2f})")
        elif neo_count >= 2:
            recommendations.append(f"💡 检测到{neo_count}个候选Neo-表位,可能具有分子胶特征")
        elif neo_count == 0:
            recommendations.append("⚠️ 未检测到Neo-表位,不符合典型分子胶特征")
    
    # 平衡度建议
    if balance_factor < 0.3:
        weak_side = "E3 Ligase" if p1_score['contact_count'] < p2_score['contact_count'] else "Substrate"
        recommendations.append(f"⚠️ {weak_side}侧接触偏少,建议增强")
    elif balance_factor > 0.7:
        recommendations.append("✅ 相互作用平衡良好")
    
    # 氢键建议
    p1_hb = p1_score['interaction_counts'].get('hbond', 0)
    p2_hb = p2_score['interaction_counts'].get('hbond', 0)
    
    if p1_hb < 2:
        recommendations.append("💡 E3侧氢键偏少(<2),考虑引入极性基团")
    if p2_hb < 2:
        recommendations.append("💡 Substrate侧氢键偏少(<2),考虑引入极性基团")
    
    # 协同效应
    if cooperativity < -3.0:
        recommendations.append(f"🚀 强协同效应 ({cooperativity:.1f} kcal/mol),设计优秀")
    elif cooperativity < -1.0:
        recommendations.append(f"✨ 显著协同效应 ({cooperativity:.1f} kcal/mol)")
    elif cooperativity > 2.0:
        recommendations.append(f"❌ 大幅不平衡惩罚 (+{cooperativity:.1f} kcal/mol),需优化")
    
    # 总分评估
    total = p1_score['total'] + p2_score['total'] + cooperativity
    if total < -15.0:
        recommendations.append("🏆 预测结合能极优 (<-15 kcal/mol),强烈建议实验验证")
    elif total < -10.0:
        recommendations.append("🎯 预测结合能优秀,建议实验验证")
    elif total > -5.0:
        recommendations.append("⚠️ 预测结合能偏弱,建议优化设计")
    
    # 针对分子胶的特殊建议
    if mechanism == "Molecular Glue":
        if ppi_result and ppi_result.get('interface_contacts', 0) < 8:
            recommendations.append("💡 分子胶提示: 考虑进一步增强蛋白-蛋白界面")
        if neo_epitope_result and neo_epitope_result.get('neo_epitope_count', 0) < 3:
            recommendations.append("💡 分子胶提示: Neo-表位较少,考虑优化底物结合")
    
    return "\n  ".join(recommendations) if recommendations else "无特殊建议"


def format_score_report(score_dict, mode="binary"):
    """
    格式化评分报告(用于打印或GUI显示)
    
    参数:
        score_dict: calculate_binary_score() 或 calculate_ternary_score() 返回值
        mode: "binary" 或 "ternary"
    
    返回:
        str: 格式化的报告文本
    """
    if not score_dict:
        return "评分计算失败"
    
    lines = []
    lines.append("=" * 60)
    lines.append("结合能评分报告 (Binding Energy Score)")
    lines.append("=" * 60)
    
    if mode == "binary":
        lines.append(f"\n总估算结合能: {score_dict['total']:.2f} kcal/mol")
        lines.append("\n能量组成 (Components):")
        lines.append("-" * 40)
        
        comp = score_dict['components']
        counts = score_dict['interaction_counts']
        
        lines.append(f"  氢键 (H-bond):        {comp['hbond']:>7.2f} kcal/mol  (n={counts.get('hbond', 0)})")
        lines.append(f"  盐桥 (Ionic):         {comp['ionic']:>7.2f} kcal/mol  (n={counts.get('ionic', 0)})")
        lines.append(f"  疏水 (Hydrophobic):   {comp['hydrophobic']:>7.2f} kcal/mol  (n={counts.get('hydrophobic', 0)})")
        lines.append(f"  π-π 堆积:             {comp['pi_pi']:>7.2f} kcal/mol  (n={counts.get('pi_pi', 0)})")
        lines.append(f"  π-阳离子:             {comp['pi_cation']:>7.2f} kcal/mol  (n={counts.get('pi_cation', 0)})")
        lines.append(f"  金属配位 (Metal):     {comp['metal']:>7.2f} kcal/mol  (n={counts.get('metal', 0)})")
        lines.append(f"  去溶剂化惩罚:         +{comp['desolvation']:>6.2f} kcal/mol")
        
    elif mode == "ternary":
        lines.append(f"\n三元复合物总结合能: {score_dict['total']:.2f} kcal/mol")
        
        # 机制判定
        mechanism = score_dict.get('mechanism', 'Unknown')
        if mechanism == "Molecular Glue":
            lines.append("机制类型: ✨ 分子胶 (Molecular Glue)")
        elif mechanism == "PROTAC":
            lines.append("机制类型: 🔗 PROTAC (Linker-based)")
        else:
            lines.append("机制类型: ❓ 未确定")
        
        lines.append(f"相互作用平衡度: {score_dict['balance_factor']:.3f}")
        lines.append(f"协同效应: {score_dict['cooperativity']:+.2f} kcal/mol")
        
        # 协同效应分解（如果有）
        if 'cooperativity_breakdown' in score_dict:
            breakdown = score_dict['cooperativity_breakdown']
            lines.append("\n协同效应分解:")
            lines.append(f"  基础: {breakdown.get('base', 0):+.2f} kcal/mol")
            if breakdown.get('ppi_bonus', 0) != 0:
                lines.append(f"  PPI奖励: {breakdown['ppi_bonus']:+.2f} kcal/mol")
            if breakdown.get('neo_epitope_bonus', 0) != 0:
                lines.append(f"  Neo-表位奖励: {breakdown['neo_epitope_bonus']:+.2f} kcal/mol")
            if breakdown.get('interface_strength', 0) != 0:
                lines.append(f"  界面强度奖励: {breakdown['interface_strength']:+.2f} kcal/mol")
            if breakdown.get('bsa_bonus', 0) != 0:
                lines.append(f"  BSA奖励: {breakdown['bsa_bonus']:+.2f} kcal/mol")
        
        lines.append("\n" + "-" * 60)
        lines.append("Protein 1 (E3 Ligase) - Ligand:")
        p1 = score_dict['protein1_score']
        lines.append(f"  总分: {p1['total']:.2f} kcal/mol  (接触数: {p1['contact_count']})")
        lines.append(f"  氢键: {p1['components']['hbond']:.2f}  盐桥: {p1['components']['ionic']:.2f}  疏水: {p1['components']['hydrophobic']:.2f}")
        
        lines.append("\n" + "-" * 60)
        lines.append("Protein 2 (Substrate) - Ligand:")
        p2 = score_dict['protein2_score']
        lines.append(f"  总分: {p2['total']:.2f} kcal/mol  (接触数: {p2['contact_count']})")
        lines.append(f"  氢键: {p2['components']['hbond']:.2f}  盐桥: {p2['components']['ionic']:.2f}  疏水: {p2['components']['hydrophobic']:.2f}")
        
        lines.append("\n" + "=" * 60)
        lines.append("设计建议:")
        lines.append("  " + score_dict['recommendation'])
    
    lines.append("\n" + "=" * 60)
    lines.append("⚠️ 注意: 此评分为快速估算,适合初步排序,不能替代精确计算")
    lines.append("=" * 60)
    
    return "\n".join(lines)


# ========== PyMOL命令封装 ==========
def score_protein_ligand(obj_name=None, ligand_resname=None, show_report=True):
    """
    PyMOL命令: 计算蛋白-配体结合能评分
    
    用法:
        score_protein_ligand protein_obj, LIG
    """
    from .interaction_analyzer import analyze_protein_ligand_interactions
    
    result = analyze_protein_ligand_interactions(obj_name, ligand_resname)
    if not result:
        print("[score_protein_ligand] 分析失败")
        return None
    
    score = calculate_binary_score(result)
    
    if show_report:
        print(format_score_report(score, mode="binary"))
    
    return score


def score_ternary_complex(obj_name=None, ligand_resname=None,
                         protein1_chains=None, protein2_chains=None,
                         show_report=True):
    """
    PyMOL命令: 计算三元复合物结合能评分
    
    用法:
        score_ternary_complex complex_obj, PROTAC, protein1_chains=[A], protein2_chains=[B]
    """
    from .interaction_analyzer import analyze_ternary_complex
    
    result = analyze_ternary_complex(
        obj_name=obj_name,
        ligand_resname=ligand_resname,
        protein1_chains=protein1_chains,
        protein2_chains=protein2_chains
    )
    
    if not result:
        print("[score_ternary_complex] 分析失败")
        return None
    
    # 重新组织结果格式以匹配评分函数
    ternary_input = {
        'protein1_result': {
            'interactions': result.get('protein1_interactions', [])
        },
        'protein2_result': {
            'interactions': result.get('protein2_interactions', [])
        }
    }
    
    score = calculate_ternary_score(ternary_input)
    
    if show_report:
        print(format_score_report(score, mode="ternary"))
    
    return score


if __name__ == "__main__":
    print("[binding_score] 这是一个PyMOL插件模块,请在PyMOL中加载使用")
    print("[binding_score] 命令: score_protein_ligand, score_ternary_complex")
