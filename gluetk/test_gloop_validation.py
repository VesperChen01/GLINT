# -*- coding: utf-8 -*-
"""
测试 G-loop CRBN 氢键验证功能
Test G-loop CRBN H-bond Validation Functions

基于已知文献结构验证功能准确性
"""

from pymol import cmd

# 测试案例（来自 Annual Review 2023）
TEST_CASES = [
    {
        'name': 'GSPT1',
        'pdb': '6H0G',
        'substrate_chain': 'A',
        'substrate_range': '60-67',
        'crbn_chain': 'E',
        'glue_resname': 'CC9',  # Lenalidomide derivative
        'expected_hbonds': 3,
        'expected_seq': 'IKIIQGAK'
    },
    {
        'name': 'CK1α',
        'pdb': '5FQD',
        'substrate_chain': 'A',
        'substrate_range': '36-43',
        'crbn_chain': 'B',
        'glue_resname': 'LEN',  # Lenalidomide
        'expected_hbonds': 3,
        'expected_seq': 'INTGSEQR'
    },
    {
        'name': 'IKZF3 (Aiolos)',
        'pdb': '6H0F',
        'substrate_chain': 'A',
        'substrate_range': '143-150',
        'crbn_chain': 'E',
        'glue_resname': 'CC9',
        'expected_hbonds': 3,
        'expected_seq': 'QAVQGDGP'
    }
]


def test_single_case(case):
    """测试单个案例"""
    print("\n" + "="*80)
    print(f"测试案例: {case['name']} (PDB: {case['pdb']})")
    print("="*80)
    
    # 1. Fetch PDB
    print(f"\n1️⃣ 下载 PDB: {case['pdb']}")
    try:
        cmd.fetch(case['pdb'], async_=0)
        print(f"   ✅ {case['pdb']} 已加载")
    except Exception as e:
        print(f"   ❌ 下载失败: {e}")
        return False
    
    # 2. 检测 G-motif
    print(f"\n2️⃣ 检测 G-motif")
    try:
        from g_motif_analyzer import find_crbn_g_motif
        hits = find_crbn_g_motif(
            obj_name=case['pdb'],
            rmsd_cutoff=3.5,
            require_gly_pos6=True,
            topk_debug=3
        )
        
        # 查找匹配的 hit
        target_hit = None
        for chain, start, end, seq, rmsd in hits:
            if chain == case['substrate_chain']:
                if seq == case['expected_seq']:
                    target_hit = (chain, start, end, seq, rmsd)
                    print(f"   ✅ 找到目标 G-loop: {chain}:{start}-{end} {seq} (RMSD={rmsd:.2f}Å)")
                    break
        
        if not target_hit:
            print(f"   ⚠️ 未找到预期的 G-loop {case['substrate_range']}")
            print(f"      预期序列: {case['expected_seq']}")
    except Exception as e:
        print(f"   ❌ G-motif 检测失败: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # 3. 验证 CRBN 氢键
    print(f"\n3️⃣ 验证 CRBN 氢键")
    try:
        from g_motif_analyzer import validate_crbn_hbonds
        hb_result = validate_crbn_hbonds(
            obj_name=case['pdb'],
            g_motif_chain=case['substrate_chain'],
            g_motif_resi_range=case['substrate_range'],
            crbn_chain=case['crbn_chain'],
            max_hbond_dist=3.5
        )
        
        if hb_result:
            passed = hb_result['total_hbonds'] >= case['expected_hbonds'] - 1  # 允许±1误差
            status = "✅" if passed else "⚠️"
            print(f"   {status} 氢键数: {hb_result['total_hbonds']}/{case['expected_hbonds']} (预期)")
            print(f"   {'✅' if hb_result['is_canonical_gloop'] else '❌'} 标准 G-loop: {hb_result['is_canonical_gloop']}")
            
            return passed
        else:
            print("   ❌ 氢键验证返回 None")
            return False
            
    except Exception as e:
        print(f"   ❌ 氢键验证失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_with_glue(case):
    """测试包含 Glue 的完整分析"""
    print("\n" + "="*80)
    print(f"完整分析: {case['name']} with Glue ({case['glue_resname']})")
    print("="*80)
    
    try:
        from g_motif_analyzer import analyze_g_motif_glue_binding
        result = analyze_g_motif_glue_binding(
            obj_name=case['pdb'],
            g_motif_chain=case['substrate_chain'],
            g_motif_resi_range=case['substrate_range'],
            glue_resname=case['glue_resname'],
            crbn_chain=case['crbn_chain'],
            validate_hbonds=True
        )
        
        if result:
            print(f"\n📊 分析结果:")
            print(f"   结合模式: {result['binding_mode']}")
            print(f"   分子胶底物: {result['is_glue_substrate']}")
            print(f"   置信度: {result['confidence']}")
            
            if result.get('crbn_hbond_validation'):
                hb = result['crbn_hbond_validation']
                print(f"   CRBN 氢键: {hb['total_hbonds']}/3")
            
            if result['key_residues_engaged']:
                print(f"   关键残基参与: {len(result['key_residues_engaged'])} 个")
            
            return result['is_glue_substrate']
        else:
            print("   ❌ 分析返回 None")
            return False
            
    except Exception as e:
        print(f"   ❌ 完整分析失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def run_all_tests():
    """运行所有测试"""
    print("\n" + "🧪"*40)
    print("G-loop CRBN 氢键验证测试套件")
    print("G-loop CRBN H-bond Validation Test Suite")
    print("🧪"*40)
    
    results = []
    
    for i, case in enumerate(TEST_CASES):
        print(f"\n\n{'='*80}")
        print(f"测试 {i+1}/{len(TEST_CASES)}: {case['name']}")
        print('='*80)
        
        # 测试氢键验证
        hbond_pass = test_single_case(case)
        
        # 测试完整分析（如果 Glue 存在）
        glue_pass = False
        try:
            # 检查 Glue 是否存在
            glue_sel = f"{case['pdb']} and resn {case['glue_resname']}"
            model = cmd.get_model(glue_sel)
            if len(model.atom) > 0:
                glue_pass = test_with_glue(case)
            else:
                print(f"\n⚠️ Glue ({case['glue_resname']}) 未在结构中找到，跳过完整分析")
        except Exception as e:
            print(f"\n⚠️ 完整分析跳过: {e}")
        
        results.append({
            'case': case['name'],
            'hbond': hbond_pass,
            'glue': glue_pass if glue_pass else 'N/A'
        })
    
    # 输出总结
    print("\n\n" + "="*80)
    print("测试总结 (Test Summary)")
    print("="*80)
    print(f"{'案例':<20} {'氢键验证':<15} {'Glue分析':<15}")
    print("-"*80)
    for r in results:
        hb_status = "✅ PASS" if r['hbond'] else "❌ FAIL"
        glue_status = "✅ PASS" if r['glue'] is True else ("⚠️ N/A" if r['glue'] == 'N/A' else "❌ FAIL")
        print(f"{r['case']:<20} {hb_status:<15} {glue_status:<15}")
    
    passed = sum(1 for r in results if r['hbond'])
    total = len(results)
    print("="*80)
    print(f"总计: {passed}/{total} 通过")
    print("="*80 + "\n")
    
    return results


def quick_test_6h0g():
    """快速测试 GSPT1 (6H0G) - 最经典的案例"""
    print("\n🚀 快速测试: GSPT1 (6H0G) - 分子胶经典案例")
    
    # Fetch
    cmd.fetch('6H0G', async_=0)
    
    # 验证氢键
    from g_motif_analyzer import validate_crbn_hbonds
    result = validate_crbn_hbonds(
        obj_name='6H0G',
        g_motif_chain='A',
        g_motif_resi_range='60-67',
        crbn_chain='E'
    )
    
    if result:
        print(f"\n✅ 测试完成!")
        print(f"   标准 G-loop: {result['is_canonical_gloop']}")
        print(f"   CRBN 氢键: {result['total_hbonds']}/3")
        return result
    else:
        print("\n❌ 测试失败")
        return None


# 主入口
if __name__ == "__main__":
    # 可以单独运行快速测试或完整测试
    import sys
    
    if len(sys.argv) > 1 and sys.argv[1] == 'quick':
        quick_test_6h0g()
    else:
        run_all_tests()
