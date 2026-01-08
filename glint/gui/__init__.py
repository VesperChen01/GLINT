# -*- coding: utf-8 -*-
"""
GLINT GUI Package

This package provides the graphical user interface for GLINT,
a molecular glue discovery and analysis toolkit.
"""

# Delayed imports to prevent segmentation faults
# These will be imported when actually needed
__all__ = [
    'GLINTDialog',
    't',
    'get_lang',
    'AnalysisWorker',
    'CommonTab',
    'GMotifWorker',
    'SurfaceAnalysisWorker',
    'SurfaceSimilarityWorker',
    'PNAnalysisWorker',
]

def __getattr__(name):
    """Lazy import to prevent crashes on module load"""
    if name == 'GLINTDialog':
        from .main_window import GLINTDialog
        return GLINTDialog
    elif name == 't':
        from .utils import t
        return t
    elif name == 'get_lang':
        from .utils import get_lang
        return get_lang
    elif name == 'CommonTab':
        from .tabs.common import CommonTab
        return CommonTab
    elif name in ['AnalysisWorker', 'PNAnalysisWorker', 'GMotifWorker',
                  'SurfaceAnalysisWorker', 'SurfaceSimilarityWorker']:
        from .workers import (
            AnalysisWorker,
            PNAnalysisWorker,
            GMotifWorker,
            SurfaceAnalysisWorker,
            SurfaceSimilarityWorker
        )
        return locals()[name]
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")