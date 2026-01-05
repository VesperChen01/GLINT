# -*- coding: utf-8 -*-
"""
GlueTK GUI Tabs Package

This package contains the individual tab widgets for the GlueTK GUI.
Each tab represents a different workflow stage in molecular glue discovery.
"""

from .common import CommonTab
from .target_discovery import TargetDiscoveryTab
from .hit_identification import HitIdentificationTab
from .lead_optimization import LeadOptimizationTab
from .visualization import VisualizationTab

__all__ = [
    'CommonTab',
    'TargetDiscoveryTab',
    'HitIdentificationTab',
    'LeadOptimizationTab',
    'VisualizationTab',
]