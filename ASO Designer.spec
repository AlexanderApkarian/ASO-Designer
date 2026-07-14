# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['aso_designer.py'],
    pathex=[],
    binaries=[],
    datas=[],
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
    [],
    exclude_binaries=True,
    name='ASO Designer - by Alexander Apkarian',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch='universal2',
    codesign_identity=None,
    entitlements_file=None,
    icon=['assets/aso_designer_icon.icns'],
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='ASO Designer - by Alexander Apkarian',
)
app = BUNDLE(
    coll,
    name='ASO Designer - by Alexander Apkarian.app',
    icon='assets/aso_designer_icon.icns',
    bundle_identifier='com.alexanderapkarian.asodesigner',
)
