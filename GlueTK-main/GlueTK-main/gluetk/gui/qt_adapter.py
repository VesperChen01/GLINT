# -*- coding: utf-8 -*-
"""
GlueTK Qt Adapter
Provides a unified interface for PyQt5, PyQt6, PySide2, and PySide6.
Optimized for PyMOL environment to avoid segmentation faults.
"""

import sys
import os

# Define the order of preference
# In PyMOL 3.x, PySide6 or PyQt6 is preferred.
# In PyMOL 2.x, PyQt5 is common.
BINDING_ORDER = ['PyQt6', 'PySide6', 'PyQt5', 'PySide2']

# Globals to store modules and classes
QtCore = None
QtWidgets = None
QtGui = None
Qt = None
Signal = None
Slot = None
Property = None

# Commonly used classes will be added to this module's namespace
_COMMON_CLASSES = [
    'QDialog', 'QVBoxLayout', 'QHBoxLayout', 'QListWidget', 'QStackedWidget',
    'QWidget', 'QPushButton', 'QLabel', 'QFrame', 'QTextEdit', 'QProgressBar', 'QMessageBox',
    'QIcon', 'QPixmap', 'QColor', 'QBrush', 'QRadialGradient', 'QLinearGradient',
    'QTimer', 'QSize', 'QSettings', 'QPainter', 'QGraphicsDropShadowEffect', 'QGridLayout',
    'QScrollArea', 'QTextBrowser', 'QLineEdit', 'QComboBox', 'QSpinBox', 'QDoubleSpinBox', 'QGroupBox',
    'QTabWidget', 'QTabBar', 'QAction', 'QMainWindow', 'QFileDialog', 'QCheckBox',
    'QFormLayout', 'QAbstractItemView', 'QListWidgetItem', 'QThread'
]

def _get_loaded_binding():
    """Check if any Qt binding is already loaded in the process."""
    # Specifically check sys.modules for any of the bindings
    for binding in ['PyQt6', 'PySide6', 'PyQt5', 'PySide2']:
        if binding + '.QtWidgets' in sys.modules:
            return binding
    return None

def _try_import(binding):
    """Try to import a specific Qt binding and populate globals."""
    global QtCore, QtWidgets, QtGui, Qt, Signal, Slot, Property
    
    try:
        if binding == 'PyQt5':
            from PyQt5 import QtCore as _QtCore, QtWidgets as _QtWidgets, QtGui as _QtGui
            from PyQt5.QtCore import pyqtSignal as _Signal, pyqtSlot as _Slot, pyqtProperty as _Property
        elif binding == 'PyQt6':
            from PyQt6 import QtCore as _QtCore, QtWidgets as _QtWidgets, QtGui as _QtGui
            from PyQt6.QtCore import pyqtSignal as _Signal, pyqtSlot as _Slot, pyqtProperty as _Property
        elif binding == 'PySide2':
            from PySide2 import QtCore as _QtCore, QtWidgets as _QtWidgets, QtGui as _QtGui
            from PySide2.QtCore import Signal as _Signal, Slot as _Slot, Property as _Property
        elif binding == 'PySide6':
            from PySide6 import QtCore as _QtCore, QtWidgets as _QtWidgets, QtGui as _QtGui
            from PySide6.QtCore import Signal as _Signal, Slot as _Slot, Property as _Property
        else:
            return False
            
        QtCore = _QtCore
        QtWidgets = _QtWidgets
        QtGui = _QtGui
        Qt = _QtCore.Qt
        Signal = _Signal
        Slot = _Slot
        Property = _Property
        
        # Populate common classes into this module's namespace
        current_module = sys.modules[__name__]
        for cls_name in _COMMON_CLASSES:
            # Look in all three modules
            for mod in [QtWidgets, QtGui, QtCore]:
                if hasattr(mod, cls_name):
                    setattr(current_module, cls_name, getattr(mod, cls_name))
                    break
                    
        return True
    except ImportError:
        return False

# 1. First, check if something is already loaded
current_binding = _get_loaded_binding()

# 2. If nothing is loaded, try in order
if not current_binding:
    # Detect PyMOL version to adjust preference
    try:
        import pymol
        v = pymol.get_version()[0]
        if v.startswith('3'):
            # PyMOL 3.x prefers Qt6
            BINDING_ORDER = ['PySide6', 'PyQt6', 'PyQt5', 'PySide2']
    except:
        pass

    for binding in BINDING_ORDER:
        if _try_import(binding):
            current_binding = binding
            break
else:
    # Use the already loaded binding
    _try_import(current_binding)

if current_binding:
    print(f"[GlueTK] Using Qt binding: {current_binding}")
else:
    print("[GlueTK] Error: No Qt binding found.")

# Export
__all__ = ['QtCore', 'QtWidgets', 'QtGui', 'Qt', 'Signal', 'Slot', 'Property', 'current_binding'] + _COMMON_CLASSES
