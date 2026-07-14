# ASO Designer - by Alexander Apkarian for Windows

This folder contains the Windows build kit for ASO Designer - by Alexander Apkarian.

## Build A Double-Clickable Windows App

On a Windows computer:

1. Install Python 3.10 or newer from python.org.
2. Unzip this folder.
3. Double-click `build_windows.bat`.
4. When it finishes, use:

```text
dist\ASO Designer - by Alexander Apkarian\ASO Designer - by Alexander Apkarian.exe
```

That `.exe` can be double-clicked on Windows.

## Run Without Building

If you only want to run it from source on Windows, double-click:

```text
run_from_source_windows.bat
```

## Why The EXE Is Not Already Included

The Mac app can be built on this Mac, but a real Windows `.exe` needs to be built on Windows. PyInstaller does not reliably cross-build Windows apps from macOS.
