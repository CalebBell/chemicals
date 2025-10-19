# -*- mode: python ; coding: utf-8 -*-
import os
import glob

# Collect all data files from chemicals package
datas = []
data_dirs = ["Critical Properties", "Density", "Electrolytes", "Environment",
             "Heat Capacity", "Identifiers", "Law", "Misc", "Phase Change",
             "Reactions", "Safety", "Solubility", "Interface", "Triple Properties",
             "Thermal Conductivity", "Vapor Pressure", "Viscosity"]

for data_dir in data_dirs:
    dir_path = os.path.join('..', 'chemicals', data_dir)
    if os.path.exists(dir_path):
        datas.append((dir_path, os.path.join('chemicals', data_dir)))

a = Analysis(
    ['basic_standalone_chemicals_check.py'],
    pathex=[],
    binaries=[],
    datas=datas,
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='basic_standalone_chemicals_check',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
