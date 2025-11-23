# -*- coding: utf-8 -*-
"""
Disease Analysis GUI Component
疾病分析 GUI 组件

为 GlueTK 提供疾病-靶点分析的图形化界面

Author: Vesper
"""

from __future__ import annotations
import os
from typing import List, Dict, Optional

# Qt 兼容导入
try:
    from PyQt5.QtCore import Qt, QThread, pyqtSignal
    from PyQt5.QtWidgets import (
        QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
        QComboBox, QTableWidget, QTableWidgetItem, QGroupBox, QProgressBar,
        QMessageBox, QFileDialog, QCheckBox, QSpinBox, QTextEdit, QTabWidget,
        QHeaderView, QAbstractItemView
    )
    QT_LIB = "PyQt5"
except ImportError:
    from PyQt6.QtCore import Qt, QThread, pyqtSignal
    from PyQt6.QtWidgets import (
        QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
        QComboBox, QTableWidget, QTableWidgetItem, QGroupBox, QProgressBar,
        QMessageBox, QFileDialog, QCheckBox, QSpinBox, QTextEdit, QTabWidget,
        QHeaderView, QAbstractItemView
    )
    QT_LIB = "PyQt6"

try:
    import pandas as pd
except ImportError:
    pd = None

# 导入 Open Targets API 模块
try:
    from .open_targets_api import (
        search_disease,
        get_disease_targets,
        enrich_targets_with_e3_scores
    )
    from .disease_config import DEFAULT_OUTPUT_DIR, E3_LIGASES
except ImportError:
    # 如果相对导入失败，尝试绝对导入
    try:
        from open_targets_api import (
            search_disease,
            get_disease_targets,
            enrich_targets_with_e3_scores
        )
        from disease_config import DEFAULT_OUTPUT_DIR, E3_LIGASES
    except ImportError:
        print("⚠️ Warning: Open Targets API modules not found")
        search_disease = None
        get_disease_targets = None
        enrich_targets_with_e3_scores = None
        DEFAULT_OUTPUT_DIR = "./disease_data"
        E3_LIGASES = ["CRBN", "VHL", "MDM2", "XIAP"]


# ============================================================================
# Worker 线程：异步疾病查询
# ============================================================================

class DiseaseQueryWorker(QThread):
    """
    疾病查询工作线程（避免 GUI 卡顿）
    """
    progress = pyqtSignal(str)
    finished = pyqtSignal(object)  # DataFrame
    error = pyqtSignal(str)
    
    def __init__(self, disease_name: str, top_n: int, include_e3: bool, e3_symbol: str):
        super().__init__()
        self.disease_name = disease_name
        self.top_n = top_n
        self.include_e3 = include_e3
        self.e3_symbol = e3_symbol
    
    def run(self):
        try:
            # 步骤 1: 搜索疾病
            self.progress.emit(f"Searching for disease: {self.disease_name}...")
            candidates = search_disease(self.disease_name, max_results=5)
            
            if not candidates:
                self.error.emit(f"No disease found for '{self.disease_name}'")
                return
            
            # 使用第一个候选
            disease_id = candidates[0]['id']
            disease_name = candidates[0]['name']
            
            # 步骤 2: 获取靶点
            self.progress.emit(f"Fetching targets for {disease_name}...")
            df = get_disease_targets(disease_id, top_n=self.top_n)
            
            if df is None or df.empty:
                self.error.emit(f"No targets found for {disease_name}")
                return
            
            # 步骤 3: 添加 E3 评分
            if self.include_e3:
                self.progress.emit(f"Calculating E3 compatibility scores ({self.e3_symbol})...")
                df = enrich_targets_with_e3_scores(df, e3_symbol=self.e3_symbol)
            
            # 添加疾病信息
            df['disease_name'] = disease_name
            df['disease_id'] = disease_id
            
            self.progress.emit(f"Query complete: {len(df)} targets found")
            self.finished.emit(df)
            
        except Exception as e:
            self.error.emit(f"Query failed: {str(e)}")


# ============================================================================
# 主标签页：疾病分析
# ============================================================================

class DiseaseAnalysisTab(QWidget):
    """
    疾病分析主标签页
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        self.current_results = None  # 当前查询结果 DataFrame
        self.query_worker = None
        
        self.init_ui()
    
    def init_ui(self):
        """初始化界面"""
        main_layout = QVBoxLayout(self)
        main_layout.setSpacing(12)
        main_layout.setContentsMargins(12, 12, 12, 12)
        
        # ========== 疾病搜索区域 ==========
        search_group = QGroupBox("Disease Search")
        search_layout = QVBoxLayout(search_group)
        search_layout.setSpacing(10)
        
        # 疾病名称输入
        disease_row = QHBoxLayout()
        disease_row.addWidget(QLabel("Disease Name:"))
        self.disease_input = QLineEdit()
        self.disease_input.setPlaceholderText("e.g., multiple myeloma, breast cancer")
        self.disease_input.setMinimumHeight(32)
        disease_row.addWidget(self.disease_input, 1)
        search_layout.addLayout(disease_row)
        
        # 参数设置
        params_layout = QHBoxLayout()
        
        # Top N
        params_layout.addWidget(QLabel("Top N Targets:"))
        self.top_n_spin = QSpinBox()
        self.top_n_spin.setRange(5, 100)
        self.top_n_spin.setValue(30)
        self.top_n_spin.setMinimumHeight(32)
        params_layout.addWidget(self.top_n_spin)
        
        params_layout.addSpacing(20)
        
        # E3 连接酶
        params_layout.addWidget(QLabel("E3 Ligase:"))
        self.e3_combo = QComboBox()
        self.e3_combo.addItems(E3_LIGASES)
        self.e3_combo.setMinimumHeight(32)
        params_layout.addWidget(self.e3_combo)
        
        params_layout.addSpacing(20)
        
        # E3 评分选项
        self.e3_checkbox = QCheckBox("Include E3 Score")
        self.e3_checkbox.setChecked(True)
        params_layout.addWidget(self.e3_checkbox)
        
        params_layout.addStretch()
        search_layout.addLayout(params_layout)
        
        # 查询按钮
        btn_row = QHBoxLayout()
        self.query_btn = QPushButton("🔍 Query Targets")
        self.query_btn.setMinimumHeight(36)
        self.query_btn.setObjectName("primary_btn")
        self.query_btn.clicked.connect(self.on_query_clicked)
        btn_row.addWidget(self.query_btn)
        btn_row.addStretch()
        search_layout.addLayout(btn_row)
        
        # 进度条
        self.progress_bar = QProgressBar()
        self.progress_bar.setTextVisible(True)
        self.progress_bar.setVisible(False)
        search_layout.addWidget(self.progress_bar)
        
        main_layout.addWidget(search_group)
        
        # ========== 结果展示区域 ==========
        results_group = QGroupBox("Results")
        results_layout = QVBoxLayout(results_group)
        results_layout.setSpacing(10)
        
        # 结果表格
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(6)
        self.results_table.setHorizontalHeaderLabels([
            "Symbol", "Name", "Disease Score", "E3 Score", "Composite Score", "Actions"
        ])
        
        # 设置表格属性
        self.results_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.results_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.results_table.setSortingEnabled(True)
        
        # 设置列宽
        header = self.results_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)  # Symbol
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)  # Name
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)  # Disease Score
        header.setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)  # E3 Score
        header.setSectionResizeMode(4, QHeaderView.ResizeMode.ResizeToContents)  # Composite Score
        header.setSectionResizeMode(5, QHeaderView.ResizeMode.Fixed)  # Actions
        self.results_table.setColumnWidth(5, 120)
        
        results_layout.addWidget(self.results_table)
        
        # 导出按钮
        export_row = QHBoxLayout()
        self.export_csv_btn = QPushButton("📄 Export CSV")
        self.export_csv_btn.setMinimumHeight(32)
        self.export_csv_btn.clicked.connect(self.on_export_csv)
        self.export_csv_btn.setEnabled(False)
        export_row.addWidget(self.export_csv_btn)
        export_row.addStretch()
        results_layout.addLayout(export_row)
        
        main_layout.addWidget(results_group, 1)
        
        # ========== 状态栏 ==========
        self.status_label = QLabel("Ready")
        self.status_label.setStyleSheet("color: #666; font-size: 11px;")
        main_layout.addWidget(self.status_label)
    
    def on_query_clicked(self):
        """查询按钮点击事件"""
        disease_name = self.disease_input.text().strip()
        
        if not disease_name:
            QMessageBox.warning(self, "Input Required", "Please enter a disease name")
            return
        
        # 检查模块是否可用
        if search_disease is None or get_disease_targets is None:
            QMessageBox.critical(
                self,
                "Module Not Available",
                "Open Targets API modules are not available.\n"
                "Please ensure the required modules are installed."
            )
            return
        
        # 禁用查询按钮
        self.query_btn.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # 不确定进度
        
        # 创建并启动 Worker 线程
        top_n = self.top_n_spin.value()
        include_e3 = self.e3_checkbox.isChecked()
        e3_symbol = self.e3_combo.currentText()
        
        self.query_worker = DiseaseQueryWorker(disease_name, top_n, include_e3, e3_symbol)
        self.query_worker.progress.connect(self.on_query_progress)
        self.query_worker.finished.connect(self.on_query_finished)
        self.query_worker.error.connect(self.on_query_error)
        self.query_worker.start()
    
    def on_query_progress(self, message: str):
        """查询进度更新"""
        self.status_label.setText(message)
    
    def on_query_finished(self, df):
        """查询完成"""
        self.current_results = df
        self.populate_results_table(df)
        
        # 恢复 UI
        self.query_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        self.export_csv_btn.setEnabled(True)
        
        disease_name = df['disease_name'].iloc[0] if not df.empty else "Unknown"
        self.status_label.setText(f"Query complete: {len(df)} targets found for {disease_name}")
    
    def on_query_error(self, error_msg: str):
        """查询错误"""
        self.query_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        
        QMessageBox.warning(
            self,
            "Query Error",
            f"Failed to query disease targets:\n\n{error_msg}\n\n"
            "Please check:\n"
            "• Network connection\n"
            "• Disease name spelling\n"
            "• API availability"
        )
        
        self.status_label.setText(f"Error: {error_msg}")
    
    def populate_results_table(self, df):
        """填充结果表格"""
        if df is None or df.empty:
            self.results_table.setRowCount(0)
            return
        
        self.results_table.setRowCount(len(df))
        self.results_table.setSortingEnabled(False)  # 填充时禁用排序
        
        for i, row in df.iterrows():
            # Symbol
            self.results_table.setItem(i, 0, QTableWidgetItem(row['symbol']))
            
            # Name
            name = row.get('name', '')
            if len(name) > 50:
                name = name[:47] + "..."
            self.results_table.setItem(i, 1, QTableWidgetItem(name))
            
            # Disease Score
            disease_score = f"{row['score']:.4f}"
            self.results_table.setItem(i, 2, QTableWidgetItem(disease_score))
            
            # E3 Score
            if 'e3_score' in row and pd.notna(row['e3_score']):
                e3_score = f"{row['e3_score']:.4f}"
                self.results_table.setItem(i, 3, QTableWidgetItem(e3_score))
            else:
                self.results_table.setItem(i, 3, QTableWidgetItem("N/A"))
            
            # Composite Score
            if 'composite_score' in row and pd.notna(row['composite_score']):
                composite_score = f"{row['composite_score']:.4f}"
                self.results_table.setItem(i, 4, QTableWidgetItem(composite_score))
            else:
                self.results_table.setItem(i, 4, QTableWidgetItem("N/A"))
            
            # Actions - Load Structure 按钮
            load_btn = QPushButton("Load")
            load_btn.setMinimumHeight(28)
            load_btn.clicked.connect(lambda checked, symbol=row['symbol']: self.on_load_structure(symbol))
            self.results_table.setCellWidget(i, 5, load_btn)
        
        self.results_table.setSortingEnabled(True)  # 重新启用排序
    
    def on_load_structure(self, symbol: str):
        """加载蛋白结构到 PyMOL"""
        try:
            from pymol import cmd
            
            # 尝试从 PDB 加载
            # 这里可以扩展为从 AlphaFold 或其他数据库加载
            reply = QMessageBox.question(
                self,
                "Load Structure",
                f"Load structure for {symbol}?\n\n"
                "This will attempt to fetch the structure from PDB.\n"
                "You can also manually load from AlphaFold or other sources.",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
            )
            
            if reply == QMessageBox.StandardButton.Yes:
                # 简单示例：尝试 fetch
                # 实际应用中可以查询 UniProt ID 并从 AlphaFold 下载
                QMessageBox.information(
                    self,
                    "Load Structure",
                    f"To load structure for {symbol}:\n\n"
                    f"1. In PyMOL console, run:\n"
                    f"   fetch <PDB_ID>, {symbol}\n\n"
                    f"2. Or load AlphaFold structure:\n"
                    f"   load AF_<UniProt_ID>.pdb, {symbol}\n\n"
                    f"3. Then visualize hotspots:\n"
                    f"   ot_glue_insight protein_obj=\"{symbol}\", gene_symbol=\"{symbol}\""
                )
        except ImportError:
            QMessageBox.warning(
                self,
                "PyMOL Not Available",
                "PyMOL is not available in this environment"
            )
    
    def on_export_csv(self):
        """导出 CSV"""
        if self.current_results is None or self.current_results.empty:
            QMessageBox.warning(self, "No Data", "No results to export")
            return
        
        # 文件对话框
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export CSV",
            os.path.join(DEFAULT_OUTPUT_DIR, "disease_targets.csv"),
            "CSV Files (*.csv)"
        )
        
        if file_path:
            try:
                # 确保目录存在
                os.makedirs(os.path.dirname(file_path), exist_ok=True)
                
                # 导出
                self.current_results.to_csv(file_path, index=False)
                
                QMessageBox.information(
                    self,
                    "Export Successful",
                    f"Results exported to:\n{file_path}"
                )
                
                self.status_label.setText(f"Exported to: {file_path}")
                
            except Exception as e:
                QMessageBox.critical(
                    self,
                    "Export Failed",
                    f"Failed to export CSV:\n{str(e)}"
                )


# ============================================================================
# 测试函数
# ============================================================================

def test_gui():
    """测试 GUI（独立运行）"""
    import sys
    from PyQt5.QtWidgets import QApplication
    
    app = QApplication(sys.argv)
    
    window = DiseaseAnalysisTab()
    window.setWindowTitle("Disease Analysis - Test")
    window.resize(1000, 700)
    window.show()
    
    sys.exit(app.exec_())


if __name__ == "__main__":
    test_gui()
