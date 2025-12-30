# -*- coding: utf-8 -*-
"""
Visualization Tab Module
"""
import os
from typing import Optional

from ..qt_adapter import (
    Qt, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
    QCheckBox, QComboBox, QGroupBox, QFormLayout, QFileDialog, QMessageBox
)

from .common import CommonTab
from ..utils import t

class VisualizationTab(CommonTab):
    def __init__(self, parent):
        super().__init__(parent)
        self._esp_maps = {}  # Store map names
        
        # Initialize UI
        self.init_ui()

    def init_ui(self):
        """Create Visualization (APBS) Tab UI"""
        layout = QVBoxLayout(self)
        layout.setSpacing(10)
        
        # APBS Group
        grp = QGroupBox(t("grp_apbs"))
        form = QFormLayout(grp)
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
        
        # Target Object
        r0 = QHBoxLayout()
        self.parent_window.obj_combo_apbs = QComboBox()
        self.parent_window.obj_combo_apbs.setMinimumHeight(32)
        self.parent_window.refresh_obj_apbs = QPushButton(t("refresh"))
        self.parent_window.refresh_obj_apbs.setObjectName("refresh_btn")
        self.parent_window.refresh_obj_apbs.setMinimumHeight(32)
        self.parent_window.refresh_obj_apbs.clicked.connect(self.refresh_objects)
        r0.addWidget(self.parent_window.obj_combo_apbs, 1)
        r0.addWidget(self.parent_window.refresh_obj_apbs)
        form.addRow(QLabel(t("apbs_target")), r0)
        
        # Parameters
        self.parent_window.apbs_grid = QLineEdit("1.0")
        self.parent_window.apbs_grid.setMinimumHeight(32)
        form.addRow(QLabel(t("apbs_grid")), self.parent_window.apbs_grid)
        
        self.parent_window.apbs_range = QLineEdit("-5,0,5")
        self.parent_window.apbs_range.setMinimumHeight(32)
        form.addRow(QLabel(t("apbs_range")), self.parent_window.apbs_range)
        
        # Run Buttons
        btn_row = QHBoxLayout()
        self.parent_window.btn_apbs_quick = QPushButton(t("btn_quick"))
        self.parent_window.btn_apbs_quick.setObjectName("highlight_btn")
        self.parent_window.btn_apbs_quick.setMinimumHeight(32)
        self.parent_window.btn_apbs_quick.clicked.connect(self.apbs_run_quick)
        
        self.parent_window.btn_apbs_true = QPushButton(t("btn_apbs"))
        self.parent_window.btn_apbs_true.setObjectName("highlight_btn")
        self.parent_window.btn_apbs_true.setMinimumHeight(32)
        self.parent_window.btn_apbs_true.clicked.connect(self.apbs_run_true)
        
        btn_row.addWidget(self.parent_window.btn_apbs_quick)
        btn_row.addWidget(self.parent_window.btn_apbs_true)
        btn_row.addStretch(1)
        
        # Export Group
        grp_exp = QGroupBox(t("grp_export"))
        exp_form = QFormLayout(grp_exp)
        
        # Image Settings
        exp_row1 = QHBoxLayout()
        self.parent_window.png_w = QLineEdit("3000")
        self.parent_window.png_w.setMinimumHeight(32)
        self.parent_window.png_h = QLineEdit("2000")
        self.parent_window.png_h.setMinimumHeight(32)
        self.parent_window.png_dpi = QLineEdit("300")
        self.parent_window.png_dpi.setMinimumHeight(32)
        
        exp_row1.addWidget(QLabel(t("img_w")))
        exp_row1.addWidget(self.parent_window.png_w)
        exp_row1.addWidget(QLabel(t("img_h")))
        exp_row1.addWidget(self.parent_window.png_h)
        exp_row1.addWidget(QLabel(t("img_dpi")))
        exp_row1.addWidget(self.parent_window.png_dpi)
        exp_form.addRow(exp_row1)
        
        # Background/Raytrace Settings
        exp_row2 = QHBoxLayout()
        self.parent_window.bg_white = QCheckBox(t("bg_white"))
        self.parent_window.bg_white.setChecked(True)
        self.parent_window.bg_trans = QCheckBox(t("bg_trans"))
        
        # Mutual exclusion logic
        self.parent_window.bg_white.toggled.connect(lambda v: (self.parent_window.bg_trans.setChecked(False) if v else None))
        self.parent_window.bg_trans.toggled.connect(lambda v: (self.parent_window.bg_white.setChecked(False) if v else None))
        
        self.parent_window.chk_ray = QCheckBox(t("raytrace"))
        self.parent_window.chk_ray.setChecked(True)
        
        self.parent_window.btn_viewport = QPushButton(t("btn_viewport"))
        self.parent_window.btn_viewport.setObjectName("refresh_btn")
        self.parent_window.btn_viewport.setMinimumHeight(32)
        self.parent_window.btn_viewport.clicked.connect(self.fill_viewport_size)
        
        exp_row2.addWidget(self.parent_window.bg_white)
        exp_row2.addWidget(self.parent_window.bg_trans)
        exp_row2.addWidget(self.parent_window.chk_ray)
        exp_row2.addStretch(1)
        exp_row2.addWidget(self.parent_window.btn_viewport)
        exp_form.addRow(exp_row2)
        
        # Export Actions
        exp_row3 = QHBoxLayout()
        self.parent_window.btn_export_png = QPushButton(t("btn_export_png"))
        self.parent_window.btn_export_png.setObjectName("save_btn")
        self.parent_window.btn_export_png.setMinimumHeight(32)
        self.parent_window.btn_export_png.clicked.connect(self.export_png)
        
        self.parent_window.btn_export_dx = QPushButton(t("btn_export_dx"))
        self.parent_window.btn_export_dx.setObjectName("save_btn")
        self.parent_window.btn_export_dx.setMinimumHeight(32)
        self.parent_window.btn_export_dx.clicked.connect(self.export_dx)
        
        exp_row3.addWidget(self.parent_window.btn_export_png)
        exp_row3.addWidget(self.parent_window.btn_export_dx)
        exp_row3.addStretch(1)
        exp_form.addRow(exp_row3)
        
        layout.addWidget(grp)
        layout.addLayout(btn_row)
        layout.addWidget(grp_exp)
        layout.addStretch(1)

    def apbs_run_quick(self):
        try:
            obj = self.parent_window.obj_combo_apbs.currentText().strip()
            if not obj or obj == t("no_object"):
                QMessageBox.warning(self, t("title"), t("no_object"))
                return
            try:
                grid = float(self.parent_window.apbs_grid.text().strip() or "1.0")
            except Exception:
                grid = 1.0
            rng_text = (self.parent_window.apbs_range.text().strip() or "-5,0,5")
            try:
                vmin, v0, vmax = [float(x) for x in rng_text.split(",")]
            except Exception:
                vmin, v0, vmax = -5.0, 0.0, 5.0

            from pymol import cmd
            map_name = f"{obj}_esp_map"
            ramp_name = f"{obj}_esp_ramp"
            cmd.map_new(map_name, "coulomb", grid, obj)
            cmd.ramp_new(ramp_name, map_name, [vmin, v0, vmax], ["blue", "white", "red"])
            cmd.show("surface", obj)
            cmd.color(ramp_name, obj)
            cmd.set("surface_quality", 1, obj)
            cmd.set("surface_color_smoothing", 1, obj)
            self._esp_maps[obj] = (map_name, ramp_name)
            self.log(f"🔷 Quick ESP: map={map_name}, ramp={ramp_name}, range=({vmin},{v0},{vmax}), grid={grid} Å")
        except Exception as e:
            self.on_error(str(e))

    def apbs_run_true(self):
        try:
            obj = self.parent_window.obj_combo_apbs.currentText().strip()
            if not obj or obj == t("no_object"):
                QMessageBox.warning(self, t("title"), t("no_object"))
                return

            apbs_tools = None
            try:
                from pymol.plugins import apbs_tools as _apbs
                apbs_tools = _apbs
            except Exception:
                try:
                    import apbs_tools as _apbs
                    apbs_tools = _apbs
                except Exception:
                    apbs_tools = None

            if apbs_tools is None or not hasattr(apbs_tools, "run_apbs"):
                self.log("APBS tools unavailable (apbs_tools.run_apbs not found). Falling back to Quick mode.")
                self.apbs_run_quick()
                return

            try:
                apbs_tools.run_apbs(selection=obj)
                self.log("APBS job submitted; if no visualization appears, check external paths in APBS Tools.")
            except Exception as ee:
                self.log(f"APBS call failed: {ee}. Falling back to Quick mode.")
                self.apbs_run_quick()
        except Exception as e:
            self.on_error(str(e))

    def fill_viewport_size(self):
        try:
            from pymol import cmd
            w, h = cmd.get_viewport()
            if w and h:
                self.parent_window.png_w.setText(str(int(w)))
                self.parent_window.png_h.setText(str(int(h)))
                self.log(f"Viewport: {w}x{h}px → filled into export settings")
            else:
                self.log("Failed to get viewport size; kept defaults")
        except Exception as e:
            self.on_error(str(e))

    def _export_png_for_object(self, obj: str):
        """Export PNG based on export panel settings."""
        from pymol import cmd
        fn, _ = QFileDialog.getSaveFileName(self, t("btn_export_png"), f"{obj}_gmotif_esp.png", "PNG (*.png)")
        if not fn: return
        if not fn.lower().endswith(".png"): fn += ".png"
        
        try:
            W = int(float(self.parent_window.png_w.text().strip()))
            H = int(float(self.parent_window.png_h.text().strip()))
        except Exception:
            W, H = 3000, 2000
        try:
            dpi = int(float(self.parent_window.png_dpi.text().strip()))
        except Exception:
            dpi = 300
            
        ray = 1 if self.parent_window.chk_ray.isChecked() else 0
        want_trans = self.parent_window.bg_trans.isChecked()
        want_white = self.parent_window.bg_white.isChecked() or not want_trans

        old_bg = cmd.get("bg_rgb")
        old_ray_bg = cmd.get("ray_opaque_background")
        if want_trans:
            cmd.bg_color("white"); cmd.set("ray_opaque_background", 0)
        elif want_white:
            cmd.bg_color("white"); cmd.set("ray_opaque_background", 1)

        cmd.png(fn, width=W, height=H, dpi=dpi, ray=ray)
        self.log(f"PNG export done: {os.path.basename(fn)} | {W}x{H}px @ {dpi} dpi | ray={ray} | bg={'transparent' if want_trans else 'white'}")

        # Restore background
        try:
            if isinstance(old_bg, (list, tuple)) and len(old_bg) == 3:
                r, g, b = [max(0.0, min(1.0, float(c))) for c in old_bg]
                cmd.set("bg_rgb", [r, g, b])
            cmd.set("ray_opaque_background", int(old_ray_bg))
        except Exception:
            pass

    def export_png(self):
        """Explicit PNG export from APBS tab."""
        try:
            obj = self.parent_window.obj_combo_apbs.currentText().strip()
            if not obj or obj == t("no_object"):
                QMessageBox.warning(self, t("title"), t("no_object"))
                return
            self._export_png_for_object(obj)
        except Exception as e:
            self.on_error(str(e))

    def export_dx(self):
        """Export DX potential map (placeholder for now, or implementation if available)."""
        # Note: The original code had a button connected to export_dx but I didn't see the implementation in the snippet.
        # I will add a placeholder log.
        self.log("Export DX feature requires APBS tools output or manual map saving.")
        # If we used apbs_run_quick, we have a map object in PyMOL, we could save it.
        try:
            obj = self.parent_window.obj_combo_apbs.currentText().strip()
            if obj in self._esp_maps:
                map_name = self._esp_maps[obj][0]
                fn, _ = QFileDialog.getSaveFileName(self, t("btn_export_dx"), f"{map_name}.dx", "DX Map (*.dx)")
                if fn:
                    from pymol import cmd
                    cmd.save(fn, map_name)
                    self.log(f"Saved map {map_name} to {fn}")
            else:
                self.log("No generated map found for this object. Run Quick ESP first.")
        except Exception as e:
            self.on_error(str(e))
