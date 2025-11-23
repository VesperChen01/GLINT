# -*- coding: utf-8 -*-
"""
突变分析 GUI 组件 - 待集成到 unified_gui.py
"""

# 将以下代码添加到 unified_gui.py 的 GlueTKDialog 类中

def _create_mutation_analysis_card(self) -> QWidget:
    """创建突变分析卡片"""
    card = QGroupBox("Protein Mutation & ΔΔG Analysis")
    layout = QVBoxLayout(card)
    layout.setSpacing(8)
    
    # 输入区域
    input_grid = QGridLayout()
    input_grid.setSpacing(6)
    
    # 对象选择
    input_grid.addWidget(QLabel("Object:"), 0, 0)
    self.mut_obj_combo = QComboBox()
    self.mut_obj_combo.setFixedHeight(28)
    self.mut_refresh_btn = QPushButton("Refresh")
    self.mut_refresh_btn.setObjectName("refresh_btn")
    self.mut_refresh_btn.setFixedHeight(28)
    self.mut_refresh_btn.clicked.connect(self.refresh_objects)
    obj_row = QHBoxLayout()
    obj_row.addWidget(self.mut_obj_combo, 1)
    obj_row.addWidget(self.mut_refresh_btn)
    input_grid.addLayout(obj_row, 0, 1)
    
    # 突变输入
    input_grid.addWidget(QLabel("Mutations:"), 1, 0)
    self.mut_input = QLineEdit()
    self.mut_input.setPlaceholderText("e.g., A:23:ALA, B:45:GLY")
    self.mut_input.setFixedHeight(28)
    input_grid.addWidget(self.mut_input, 1, 1)
    
    # 方法选择
    input_grid.addWidget(QLabel("Method:"), 2, 0)
    method_row = QHBoxLayout()
    self.mut_method_combo = QComboBox()
    self.mut_method_combo.addItems(["Auto", "FoldX", "PyRosetta"])
    self.mut_method_combo.setFixedHeight(28)
    self.mut_method_combo.setFixedWidth(120)
    method_row.addWidget(self.mut_method_combo)
    method_row.addStretch()
    input_grid.addLayout(method_row, 2, 1)
    
    layout.addLayout(input_grid)
    
    # 按钮行
    btn_row = QHBoxLayout()
    btn_row.setSpacing(6)
    
    self.mut_perform_btn = QPushButton("Perform Mutation")
    self.mut_perform_btn.setObjectName("primary_btn")
    self.mut_perform_btn.setFixedHeight(32)
    self.mut_perform_btn.clicked.connect(self.run_mutation)
    
    self.mut_minimize_btn = QPushButton("Minimize Energy")
    self.mut_minimize_btn.setObjectName("highlight_btn")
    self.mut_minimize_btn.setFixedHeight(32)
    self.mut_minimize_btn.clicked.connect(self.run_minimize)
    
    self.mut_analyze_btn = QPushButton("Full Analysis")
    self.mut_analyze_btn.setObjectName("highlight_btn")
    self.mut_analyze_btn.setFixedHeight(32)
    self.mut_analyze_btn.clicked.connect(self.run_mutation_analysis)
    
    btn_row.addWidget(self.mut_perform_btn)
    btn_row.addWidget(self.mut_minimize_btn)
    btn_row.addWidget(self.mut_analyze_btn)
    btn_row.addStretch()
    
    layout.addLayout(btn_row)
    
    return card


# 事件处理函数

def run_mutation(self):
    """执行突变"""
    obj_name = self.mut_obj_combo.currentText()
    mutations_str = self.mut_input.text().strip()
    
    if not obj_name or not mutations_str:
        self.log("❌ 请选择对象并输入突变")
        return
    
    try:
        # 解析突变字符串
        mutations = []
        for mut in mutations_str.split(','):
            mut = mut.strip()
            if ':' in mut:
                parts = mut.split(':')
                if len(parts) == 3:
                    mutations.append((parts[0], parts[1], parts[2]))
        
        if not mutations:
            self.log("❌ 突变格式错误，示例: A:23:ALA, B:45:GLY")
            return
        
        self.log(f"🧬 执行突变: {obj_name}")
        
        # 调用突变分析模块
        from pymol import cmd
        cmd.perform_mutation(obj_name, mutations, method='pymol')
        
        self.log("✅ 突变完成")
        
    except Exception as e:
        self.log(f"❌ 突变失败: {e}")
        import traceback
        traceback.print_exc()


def run_minimize(self):
    """能量最小化"""
    obj_name = self.mut_obj_combo.currentText()
    
    if not obj_name:
        self.log("❌ 请选择对象")
        return
    
    try:
        self.log(f"⚡ 能量最小化: {obj_name}")
        
        from pymol import cmd
        cmd.minimize_energy(obj_name, cycles=100)
        
        self.log("✅ 最小化完成")
        
    except Exception as e:
        self.log(f"❌ 最小化失败: {e}")
        import traceback
        traceback.print_exc()


def run_mutation_analysis(self):
    """完整突变分析"""
    obj_name = self.mut_obj_combo.currentText()
    mutations_str = self.mut_input.text().strip()
    method = self.mut_method_combo.currentText().lower()
    
    if not obj_name or not mutations_str:
        self.log("❌ 请选择对象并输入突变")
        return
    
    try:
        # 解析突变
        mutations = []
        for mut in mutations_str.split(','):
            mut = mut.strip()
            if ':' in mut:
                parts = mut.split(':')
                if len(parts) == 3:
                    mutations.append((parts[0], parts[1], parts[2]))
        
        if not mutations:
            self.log("❌ 突变格式错误")
            return
        
        self.log(f"🧬 开始完整分析: {obj_name}")
        self.log(f"突变数量: {len(mutations)}")
        self.log(f"方法: {method}")
        
        # 调用完整分析
        from pymol import cmd
        cmd.analyze_mutation_effects(obj_name, mutations, method=method)
        
        self.log("✅ 分析完成")
        
    except Exception as e:
        self.log(f"❌ 分析失败: {e}")
        import traceback
        traceback.print_exc()
