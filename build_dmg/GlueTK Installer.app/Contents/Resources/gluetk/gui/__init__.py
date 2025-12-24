# -*- coding: utf-8 -*-
"""
GlueTK GUI Package

This package provides the graphical user interface for GlueTK,
a molecular glue discovery and analysis toolkit.
"""

# Delayed imports to prevent segmentation faults
# These will be imported when actually needed
__all__ = [
    'GlueTKDialog',
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
    if name == 'GlueTKDialog':
        from .main_window import GlueTKDialog
        return GlueTKDialog
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