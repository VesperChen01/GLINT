# -*- coding: utf-8 -*-
"""
MolStruct 现代化样式表
灵感来自 Warp Plus 设计 - 时尚感与科研风格并重
"""

# 配色方案 - 蓝灰主题 + 科研配色
COLORS = {
    # 主色调
    "primary": "#3b82f6",        # 明亮蓝
    "primary_hover": "#2563eb",  # 深蓝
    "primary_pressed": "#1d4ed8",
    
    # 次要色
    "secondary": "#8b5cf6",      # 紫色
    "accent": "#06b6d4",         # 青色
    
    # 中性色 - 蓝灰系
    "bg_main": "#0f172a",        # 深蓝灰背景
    "bg_card": "#1e293b",        # 卡片背景
    "bg_hover": "#334155",       # 悬停背景
    "bg_light": "#f8fafc",       # 浅色背景
    
    # 文字
    "text_primary": "#f1f5f9",   # 主文字（深色模式）
    "text_secondary": "#94a3b8", # 次要文字
    "text_dark": "#0f172a",      # 深色文字（浅色模式）
    "text_muted": "#64748b",     # 减弱文字
    
    # 功能色
    "success": "#10b981",        # 成功/分析
    "warning": "#f59e0b",        # 警告/渲染
    "danger": "#ef4444",         # 危险/关闭
    "info": "#3b82f6",           # 信息
    
    # 边框
    "border": "#334155",
    "border_light": "#e2e8f0",
    
    # 特殊
    "glow_blue": "rgba(59, 130, 246, 0.3)",
    "glow_purple": "rgba(139, 92, 246, 0.3)",
}

def get_modern_stylesheet(dark_mode=True):
    """
    返回现代化样式表
    dark_mode: True=深色主题（默认）, False=浅色主题
    """
    
    if dark_mode:
        return get_dark_theme()
    else:
        return get_light_theme()


def get_dark_theme():
    """深色主题 - 科技感"""
    return f"""
    /* ========== 全局样式 ========== */
    QDialog {{
        background: qlineargradient(
            x1:0, y1:0, x2:0, y2:1,
            stop:0 {COLORS['bg_main']},
            stop:1 #0a0f1e
        );
        color: {COLORS['text_primary']};
        font-family: -apple-system, BlinkMacSystemFont, "SF Pro Display", 
                     "Segoe UI", "Microsoft YaHei", sans-serif;
        font-size: 13px;
    }}
    
    /* ========== 导航列表 ========== */
    QListWidget {{
        background: {COLORS['bg_card']};
        border: 1px solid {COLORS['border']};
        border-radius: 12px;
        padding: 8px;
        outline: none;
    }}
    
    QListWidget::item {{
        background: transparent;
        color: {COLORS['text_secondary']};
        padding: 12px 16px;
        margin: 4px 0;
        border-radius: 8px;
        font-weight: 500;
        font-size: 13px;
        border-left: 3px solid transparent;
    }}
    
    QListWidget::item:hover {{
        background: {COLORS['bg_hover']};
        color: {COLORS['text_primary']};
        border-left: 3px solid {COLORS['primary']};
    }}
    
    QListWidget::item:selected {{
        background: qlineargradient(
            x1:0, y1:0, x2:1, y2:0,
            stop:0 {COLORS['primary']},
            stop:1 {COLORS['secondary']}
        );
        color: white;
        font-weight: 600;
        border-left: 3px solid {COLORS['accent']};
    }}
    
    /* ========== 标题标签 ========== */
    QLabel {{
        color: {COLORS['text_primary']};
        padding: 2px 0;      /* 适中的垂直间距 */
        min-height: 20px;    /* 适中的最小高度 */
    }}
    
    QGroupBox {{
        background: {COLORS['bg_card']};
        border: 1px solid {COLORS['border']};
        border-radius: 10px;
        margin-top: 10px;
        padding: 18px 12px 10px 12px; /* 上/左右/下 */
        font-weight: 500;
        font-size: 13px;
        color: {COLORS['text_primary']};
    }}
    
    QGroupBox::title {{
        subcontrol-origin: margin;
        subcontrol-position: top left;
        padding: 0 8px;
        margin-left: 12px;
        color: {COLORS['primary']};
    }}
    
    /* ========== 输入框 ========== */
    QLineEdit {{
        background: {COLORS['bg_main']};
        border: 1px solid {COLORS['border']};
        border-radius: 6px;
        padding: 3px 8px;   /* 进一步减小垂直内边距 */
        min-height: 24px;    /* 进一步减小最小高度 */
        max-height: 24px;    /* 限制最大高度 */
        color: {COLORS['text_primary']};
        selection-background-color: {COLORS['primary']};
        font-size: 12px;     /* 减小字体 */
    }}
    
    QLineEdit:focus {{
        border: 2px solid {COLORS['primary']};
        background: #1a1f2e;
    }}
    
    QLineEdit:disabled {{
        background: #141824;
        color: {COLORS['text_muted']};
    }}
    
    /* ========== 按钮样式 ========== */
    QPushButton {{
        background: {COLORS['bg_card']};
        color: {COLORS['text_primary']};
        border: 1px solid {COLORS['border']};
        border-radius: 8px;
        padding: 4px 14px; /* 进一步减小按钮垂直内边距 */
        min-height: 24px;    /* 限制按钮高度 */
        max-height: 26px;
        font-weight: 500;
        font-size: 12px;     /* 减小字体 */
    }}
    
    QPushButton:hover {{
        background: {COLORS['bg_hover']};
        border: 1px solid {COLORS['primary']};
    }}
    
    QPushButton:pressed {{
        background: {COLORS['bg_main']};
    }}
    
    QPushButton:disabled {{
        background: #141824;
        color: {COLORS['text_muted']};
        border: 1px solid #1e293b;
    }}
    
    /* 主要按钮 - 分析 */
    QPushButton#analyze_btn,
    QPushButton#gm_btn,
    QPushButton#pl_analyze_btn,
    QPushButton#tc_analyze_btn,
    QPushButton#ap_analyze_btn {{
        background: qlineargradient(
            x1:0, y1:0, x2:1, y2:0,
            stop:0 {COLORS['success']},
            stop:1 #059669
        );
        color: white;
        border: none;
        font-weight: 600;
        padding: 8px 20px; /* 减小主要按钮垂直内边距 */
    }}
    
    QPushButton#analyze_btn:hover,
    QPushButton#gm_btn:hover,
    QPushButton#pl_analyze_btn:hover,
    QPushButton#tc_analyze_btn:hover,
    QPushButton#ap_analyze_btn:hover {{
        background: qlineargradient(
            x1:0, y1:0, x2:1, y2:0,
            stop:0 #059669,
            stop:1 #047857
        );
        box-shadow: 0 0 20px {COLORS['glow_blue']};
    }}
    
    /* 渲染按钮 - 警告色 */
    QPushButton#render_btn,
    QPushButton#gm_btn_render,
    QPushButton#pl_vis_btn,
    QPushButton#ap_vis_btn {{
        background: qlineargradient(
            x1:0, y1:0, x2:1, y2:0,
            stop:0 {COLORS['warning']},
            stop:1 #d97706
        );
        color: white;
        border: none;
        font-weight: 600;
        padding: 12px 24px;
    }}
    
    QPushButton#render_btn:hover,
    QPushButton#gm_btn_render:hover,
    QPushButton#pl_vis_btn:hover,
    QPushButton#ap_vis_btn:hover {{
        background: qlineargradient(
            x1:0, y1:0, x2:1, y2:0,
            stop:0 #d97706,
            stop:1 #b45309
        );
    }}
    
    /* 网络/图表按钮 - 紫色 */
    QPushButton#pl_net_btn,
    QPushButton#tc_net_btn {{
        background: qlineargradient(
            x1:0, y1:0, x2:1, y2:0,
            stop:0 {COLORS['secondary']},
            stop:1 #7c3aed
        );
        color: white;
        border: none;
        font-weight: 600;
        padding: 12px 24px;
    }}
    
    QPushButton#pl_net_btn:hover,
    QPushButton#tc_net_btn:hover {{
        background: qlineargradient(
            x1:0, y1:0, x2:1, y2:0,
            stop:0 #7c3aed,
            stop:1 #6d28d9
        );
    }}
    
    /* 浏览按钮 - 蓝色 */
    QPushButton[text="Browse…"],
    QPushButton[text="浏览…"] {{
        background: {COLORS['primary']};
        color: white;
        border: none;
        font-weight: 500;
    }}
    
    QPushButton[text="Browse…"]:hover,
    QPushButton[text="浏览…"]:hover {{
        background: {COLORS['primary_hover']};
    }}
    
    /* 关闭按钮 - 红色 */
    QPushButton#close_btn {{
        background: {COLORS['danger']};
        color: white;
        border: none;
        font-weight: 600;
    }}
    
    QPushButton#close_btn:hover {{
        background: #dc2626;
    }}
    
    /* 高亮/清除按钮 */
    QPushButton#highlight_btn {{
        background: {COLORS['primary']};
        color: white;
        border: none;
        font-weight: 600;
    }}
    
    QPushButton#clear_btn {{
        background: {COLORS['text_muted']};
        color: white;
        border: none;
    }}
    
    /* 主题切换按钮 */
    QPushButton#theme_toggle_btn {{
        background: qlineargradient(
            x1:0, y1:0, x2:1, y2:0,
            stop:0 {COLORS['accent']},
            stop:1 {COLORS['primary']}
        );
        color: white;
        border: none;
        font-weight: 600;
        padding: 10px 16px;
        font-size: 13px;
    }}
    
    QPushButton#theme_toggle_btn:hover {{
        background: qlineargradient(
            x1:0, y1:0, x2:1, y2:0,
            stop:0 {COLORS['primary']},
            stop:1 {COLORS['secondary']}
        );
        box-shadow: 0 0 15px {COLORS['glow_blue']};
    }}
    
    /* ========== 下拉框 ========== */
    QComboBox {{
        background: {COLORS['bg_main']};
        border: 1px solid {COLORS['border']};
        border-radius: 6px;
        padding: 3px 8px;   /* 进一步减小垂直内边距 */
        min-height: 24px;    /* 进一步减小最小高度 */
        max-height: 24px;    /* 限制最大高度 */
        color: {COLORS['text_primary']};
        font-size: 12px;     /* 减小字体 */
        min-width: 120px;
    }}
    
    QComboBox:hover {{
        border: 1px solid {COLORS['primary']};
    }}
    
    QComboBox::drop-down {{
        border: none;
        width: 30px;
    }}
    
    QComboBox::down-arrow {{
        image: none;
        border-left: 5px solid transparent;
        border-right: 5px solid transparent;
        border-top: 5px solid {COLORS['text_secondary']};
        margin-right: 8px;
    }}
    
    QComboBox QAbstractItemView {{
        background: {COLORS['bg_card']};
        border: 1px solid {COLORS['border']};
        border-radius: 8px;
        selection-background-color: {COLORS['primary']};
        selection-color: white;
        padding: 4px;
    }}
    
    /* ========== 复选框 ========== */
    QCheckBox {{
        color: {COLORS['text_primary']};
        spacing: 8px;
    }}
    
    QCheckBox::indicator {{
        width: 18px;
        height: 18px;
        border: 2px solid {COLORS['border']};
        border-radius: 4px;
        background: {COLORS['bg_main']};
    }}
    
    QCheckBox::indicator:hover {{
        border: 2px solid {COLORS['primary']};
    }}
    
    QCheckBox::indicator:checked {{
        background: {COLORS['primary']};
        border: 2px solid {COLORS['primary']};
        image: url(data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iMTIiIGhlaWdodD0iOSIgdmlld0JveD0iMCAwIDEyIDkiIGZpbGw9Im5vbmUiIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyI+PHBhdGggZD0iTTEgNEw0LjUgNy41TDExIDEiIHN0cm9rZT0id2hpdGUiIHN0cm9rZS13aWR0aD0iMiIgc3Ryb2tlLWxpbmVjYXA9InJvdW5kIiBzdHJva2UtbGluZWpvaW49InJvdW5kIi8+PC9zdmc+);
    }}
    
    /* ========== 表格 ========== */
    QTableWidget {{
        background: {COLORS['bg_card']};
        alternate-background-color: {COLORS['bg_main']};
        gridline-color: {COLORS['border']};
        border: 1px solid {COLORS['border']};
        border-radius: 8px;
        color: {COLORS['text_primary']};
        selection-background-color: {COLORS['primary']};
        selection-color: white;
    }}
    
    QTableWidget::item {{
        padding: 10px 8px;  /* 适中的垂直内边距 */
        min-height: 34px;   /* 适中的最小行高 */
        color: {COLORS['text_primary']};
    }}
    
    QTableWidget::item:selected {{
        background: {COLORS['primary']};
        color: white;
    }}
    
    QHeaderView::section {{
        background: {COLORS['bg_main']};
        color: {COLORS['text_secondary']};
        padding: 11px 10px;  /* 适中的内边距 */
        min-height: 36px;    /* 适中的表头高度 */
        border: none;
        border-bottom: 2px solid {COLORS['primary']};
        font-weight: 600;
        font-size: 12px;
        text-transform: uppercase;
        letter-spacing: 0.5px;
    }}
    
    /* ========== 文本编辑器 ========== */
    QTextEdit {{
        background: {COLORS['bg_main']};
        border: 1px solid {COLORS['border']};
        border-radius: 8px;
        color: {COLORS['text_primary']};
        padding: 12px;
        font-family: 'SF Mono', 'Consolas', 'Monaco', monospace;
        font-size: 12px;
        selection-background-color: {COLORS['primary']};
    }}
    
    /* ========== 进度条 ========== */
    QProgressBar {{
        background: {COLORS['bg_main']};
        border: 1px solid {COLORS['border']};
        border-radius: 6px;
        text-align: center;
        color: {COLORS['text_primary']};
        height: 20px;
    }}
    
    QProgressBar::chunk {{
        background: qlineargradient(
            x1:0, y1:0, x2:1, y2:0,
            stop:0 {COLORS['primary']},
            stop:1 {COLORS['secondary']}
        );
        border-radius: 6px;
    }}
    
    /* ========== 滚动条 ========== */
    QScrollBar:vertical {{
        background: {COLORS['bg_card']};
        width: 12px;
        border-radius: 6px;
        margin: 0;
    }}
    
    QScrollBar::handle:vertical {{
        background: {COLORS['bg_hover']};
        border-radius: 6px;
        min-height: 30px;
    }}
    
    QScrollBar::handle:vertical:hover {{
        background: {COLORS['primary']};
    }}
    
    QScrollBar::add-line:vertical,
    QScrollBar::sub-line:vertical {{
        height: 0px;
    }}
    
    QScrollBar:horizontal {{
        background: {COLORS['bg_card']};
        height: 12px;
        border-radius: 6px;
    }}
    
    QScrollBar::handle:horizontal {{
        background: {COLORS['bg_hover']};
        border-radius: 6px;
        min-width: 30px;
    }}
    
    QScrollBar::handle:horizontal:hover {{
        background: {COLORS['primary']};
    }}
    
    /* ========== 标签页 ========== */
    QTabWidget::pane {{
        border: 1px solid {COLORS['border']};
        border-radius: 12px;
        background: {COLORS['bg_card']};
        padding: 12px;
    }}
    
    QTabBar::tab {{
        background: {COLORS['bg_main']};
        color: {COLORS['text_secondary']};
        padding: 12px 24px;
        margin-right: 4px;
        border-top-left-radius: 8px;
        border-top-right-radius: 8px;
        font-weight: 500;
    }}
    
    QTabBar::tab:selected {{
        background: {COLORS['bg_card']};
        color: {COLORS['primary']};
        border-bottom: 2px solid {COLORS['primary']};
    }}
    
    QTabBar::tab:hover {{
        background: {COLORS['bg_hover']};
        color: {COLORS['text_primary']};
    }}
    
    /* ========== 分隔线 ========== */
    QFrame[frameShape="4"],
    QFrame[frameShape="5"] {{
        color: {COLORS['border']};
    }}
    """


def get_light_theme():
    """浅色主题 - 清新专业"""
    return f"""
    /* ========== 全局样式 ========== */
    QDialog {{
        background: {COLORS['bg_light']};
        color: {COLORS['text_dark']};
        font-family: -apple-system, BlinkMacSystemFont, "SF Pro Display", 
                     "Segoe UI", "Microsoft YaHei", sans-serif;
        font-size: 13px;
    }}
    
    /* ========== Label ========== */
    QLabel {{
        color: {COLORS['text_dark']};
    }}
    
    /* ========== 导航列表 ========== */
    QListWidget {{
        background: white;
        border: 1px solid {COLORS['border_light']};
        border-radius: 12px;
        padding: 8px;
        outline: none;
    }}
    
    QListWidget::item {{
        background: transparent;
        color: {COLORS['text_muted']};
        padding: 12px 16px;
        margin: 4px 0;
        border-radius: 8px;
        font-weight: 500;
        font-size: 13px;
        border-left: 3px solid transparent;
    }}
    
    QListWidget::item:hover {{
        background: #f1f5f9;
        color: {COLORS['text_dark']};
        border-left: 3px solid {COLORS['primary']};
    }}
    
    QListWidget::item:selected {{
        background: qlineargradient(
            x1:0, y1:0, x2:1, y2:0,
            stop:0 {COLORS['primary']},
            stop:1 {COLORS['secondary']}
        );
        color: white;
        font-weight: 600;
        border-left: 3px solid {COLORS['accent']};
    }}
    
    /* ========== GroupBox ========== */
    QGroupBox {{
        background: white;
        border: 1px solid {COLORS['border_light']};
        border-radius: 10px;
        margin-top: 10px;
        padding: 18px 12px 10px 12px; /* 上/左右/下 */
        font-weight: 500;
        font-size: 13px;
        color: {COLORS['text_dark']};
    }}
    
    QGroupBox::title {{
        subcontrol-origin: margin;
        subcontrol-position: top left;
        padding: 0 8px;
        margin-left: 12px;
        color: {COLORS['primary']};
    }}
    
    /* ========== 输入框 ========== */
    QLineEdit {{
        background: white;
        border: 1px solid {COLORS['border_light']};
        border-radius: 6px;
        padding: 7px 12px;   /* 更紧凑的垂直内边距 */
        min-height: 28px;    /* 更小的最小高度 */
        color: {COLORS['text_dark']};
        selection-background-color: {COLORS['primary']};
        font-size: 13px;
    }}
    
    QLineEdit:focus {{
        border: 2px solid {COLORS['primary']};
        background: white;
    }}
    
    QLineEdit:disabled {{
        background: #f1f5f9;
        color: {COLORS['text_muted']};
    }}
    
    /* ========== 按钮 ========== */
    QPushButton {{
        background: white;
        color: {COLORS['text_dark']};
        border: 1px solid {COLORS['border_light']};
        border-radius: 8px;
        padding: 10px 20px;
        font-weight: 500;
        font-size: 13px;
    }}
    
    QPushButton:hover {{
        background: #f1f5f9;
        border: 1px solid {COLORS['primary']};
    }}
    
    QPushButton:pressed {{
        background: #e2e8f0;
    }}
    
    QPushButton:disabled {{
        background: #f8fafc;
        color: #cbd5e1;
        border: 1px solid #e2e8f0;
    }}
    
    /* 主要按钮 - 分析 */
    QPushButton#analyze_btn,
    QPushButton#gm_btn,
    QPushButton#pl_analyze_btn,
    QPushButton#tc_analyze_btn,
    QPushButton#ap_analyze_btn {{
        background: qlineargradient(
            x1:0, y1:0, x2:1, y2:0,
            stop:0 {COLORS['success']},
            stop:1 #059669
        );
        color: white;
        border: none;
        font-weight: 600;
        padding: 12px 24px;
    }}
    
    QPushButton#analyze_btn:hover,
    QPushButton#gm_btn:hover,
    QPushButton#pl_analyze_btn:hover,
    QPushButton#tc_analyze_btn:hover,
    QPushButton#ap_analyze_btn:hover {{
        background: qlineargradient(
            x1:0, y1:0, x2:1, y2:0,
            stop:0 #059669,
            stop:1 #047857
        );
    }}
    
    /* 渲染按钮 */
    QPushButton#render_btn,
    QPushButton#gm_btn_render,
    QPushButton#pl_visualize_btn,
    QPushButton#ap_visualize_btn {{
        background: qlineargradient(
            x1:0, y1:0, x2:1, y2:0,
            stop:0 {COLORS['warning']},
            stop:1 #d97706
        );
        color: white;
        border: none;
        font-weight: 600;
        padding: 12px 24px;
    }}
    
    QPushButton#render_btn:hover,
    QPushButton#gm_btn_render:hover,
    QPushButton#pl_visualize_btn:hover,
    QPushButton#ap_visualize_btn:hover {{
        background: qlineargradient(
            x1:0, y1:0, x2:1, y2:0,
            stop:0 #d97706,
            stop:1 #b45309
        );
    }}
    
    /* 网络按钮 */
    QPushButton#pl_network_btn,
    QPushButton#tc_network_btn {{
        background: qlineargradient(
            x1:0, y1:0, x2:1, y2:0,
            stop:0 {COLORS['secondary']},
            stop:1 #7c3aed
        );
        color: white;
        border: none;
        font-weight: 600;
        padding: 12px 24px;
    }}
    
    /* 浏览按钮 */
    QPushButton#browse_btn {{
        background: {COLORS['primary']};
        color: white;
        border: none;
        font-weight: 500;
        min-width: 70px;
    }}
    
    QPushButton#browse_btn:hover {{
        background: {COLORS['primary_hover']};
    }}
    
    /* 刷新按钮 */
    QPushButton#refresh_btn {{
        background: #f5f5f5;
        color: #666;
        border: 1px solid #ddd;
        padding: 6px 10px;
        font-size: 11px;
        min-width: 60px;
    }}
    
    QPushButton#refresh_btn:hover {{
        background: #e8e8e8;
    }}
    
    /* 保存按钮 */
    QPushButton#save_btn {{
        background: {COLORS['accent']};
        color: white;
        border: none;
        font-weight: 500;
    }}
    
    /* 高亮按钮 */
    QPushButton#highlight_btn {{
        background: {COLORS['primary']};
        color: white;
        border: none;
        font-weight: 600;
    }}
    
    QPushButton#highlight_btn:hover {{
        background: {COLORS['primary_hover']};
    }}
    
    /* 小主题切换按钮 */
    QPushButton#theme_toggle_btn_small {{
        background: #f5f5f5;
        color: #666;
        border: 1px solid #ddd;
        border-radius: 4px;
        padding: 2px;
        font-size: 16px;
        font-weight: normal;
    }}
    
    QPushButton#theme_toggle_btn_small:hover {{
        background: #e8e8e8;
        border-color: #999;
    }}
    
    /* ========== 下拉框 ========== */
    QComboBox {{
        background: white;
        border: 1px solid {COLORS['border_light']};
        border-radius: 6px;
        padding: 7px 12px;   /* 更紧凑的垂直内边距 */
        min-height: 28px;    /* 更小的最小高度 */
        color: {COLORS['text_dark']};
        font-size: 13px;
        min-width: 120px;
    }}
    
    QComboBox:hover {{
        border: 1px solid {COLORS['primary']};
    }}
    
    QComboBox::drop-down {{
        border: none;
        width: 30px;
    }}
    
    QComboBox::down-arrow {{
        image: none;
        border-left: 5px solid transparent;
        border-right: 5px solid transparent;
        border-top: 5px solid {COLORS['text_muted']};
        margin-right: 8px;
    }}
    
    QComboBox QAbstractItemView {{
        background: white;
        border: 1px solid {COLORS['border_light']};
        border-radius: 8px;
        selection-background-color: {COLORS['primary']};
        selection-color: white;
        padding: 4px;
        color: {COLORS['text_dark']};
    }}
    
    /* ========== 复选框 ========== */
    QCheckBox {{
        color: {COLORS['text_dark']};
        spacing: 8px;
    }}
    
    QCheckBox::indicator {{
        width: 18px;
        height: 18px;
        border: 2px solid {COLORS['border_light']};
        border-radius: 4px;
        background: white;
    }}
    
    QCheckBox::indicator:hover {{
        border: 2px solid {COLORS['primary']};
    }}
    
    QCheckBox::indicator:checked {{
        background: {COLORS['primary']};
        border: 2px solid {COLORS['primary']};
    }}
    
    /* ========== 表格 ========== */
    QTableWidget {{
        background: white;
        alternate-background-color: #f8fafc;
        gridline-color: {COLORS['border_light']};
        border: 1px solid {COLORS['border_light']};
        border-radius: 8px;
        color: {COLORS['text_dark']};
        selection-background-color: {COLORS['primary']};
        selection-color: white;
    }}
    
    QTableWidget::item {{
        padding: 10px 8px;  /* 适中的垂直内边距 */
        min-height: 34px;   /* 适中的最小行高 */
        color: {COLORS['text_dark']};
        font-size: 13px;
    }}
    
    QTableWidget::item:selected {{
        background: {COLORS['primary']};
        color: white;
    }}
    
    QHeaderView::section {{
        background: #f1f5f9;
        color: {COLORS['text_dark']};
        padding: 11px 10px;  /* 适中的内边距 */
        min-height: 36px;    /* 适中的表头高度 */
        border: none;
        border-bottom: 2px solid {COLORS['primary']};
        font-weight: 600;
        font-size: 12px;
        text-transform: uppercase;
        letter-spacing: 0.5px;
    }}
    
    /* ========== 文本编辑器 ========== */
    QTextEdit {{
        background: white;
        border: 1px solid {COLORS['border_light']};
        border-radius: 8px;
        color: {COLORS['text_dark']};
        padding: 12px;
        font-family: 'SF Mono', 'Consolas', 'Monaco', monospace;
        font-size: 12px;
        selection-background-color: {COLORS['primary']};
    }}
    
    /* ========== 进度条 ========== */
    QProgressBar {{
        background: #e2e8f0;
        border: 1px solid {COLORS['border_light']};
        border-radius: 6px;
        text-align: center;
        color: {COLORS['text_dark']};
        height: 20px;
    }}
    
    QProgressBar::chunk {{
        background: qlineargradient(
            x1:0, y1:0, x2:1, y2:0,
            stop:0 {COLORS['primary']},
            stop:1 {COLORS['secondary']}
        );
        border-radius: 6px;
    }}
    
    /* ========== 滚动条 ========== */
    QScrollBar:vertical {{
        background: #f1f5f9;
        width: 12px;
        border-radius: 6px;
        margin: 0;
    }}
    
    QScrollBar::handle:vertical {{
        background: #cbd5e1;
        border-radius: 6px;
        min-height: 30px;
    }}
    
    QScrollBar::handle:vertical:hover {{
        background: {COLORS['primary']};
    }}
    
    QScrollBar::add-line:vertical,
    QScrollBar::sub-line:vertical {{
        height: 0px;
    }}
    
    QScrollBar:horizontal {{
        background: #f1f5f9;
        height: 12px;
        border-radius: 6px;
    }}
    
    QScrollBar::handle:horizontal {{
        background: #cbd5e1;
        border-radius: 6px;
        min-width: 30px;
    }}
    
    QScrollBar::handle:horizontal:hover {{
        background: {COLORS['primary']};
    }}
    """


# 使用示例
if __name__ == "__main__":
    print("=== 深色主题 ===")
    print(get_modern_stylesheet(dark_mode=True)[:500])
    print("\n=== 浅色主题 ===")
    print(get_modern_stylesheet(dark_mode=False)[:500])
