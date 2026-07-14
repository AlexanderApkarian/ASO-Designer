@echo off
setlocal

cd /d "%~dp0"

echo Creating Windows build environment...
py -3 -m venv .venv_windows
if errorlevel 1 goto :error

call .venv_windows\Scripts\activate.bat
if errorlevel 1 goto :error

echo Installing dependencies...
python -m pip install --upgrade pip
if errorlevel 1 goto :error
python -m pip install -r requirements.txt pyinstaller
if errorlevel 1 goto :error

echo Building ASO Designer - by Alexander Apkarian.exe...
pyinstaller --clean --noconfirm --windowed --name "ASO Designer - by Alexander Apkarian" --icon "assets\aso_designer_icon.ico" aso_designer.py
if errorlevel 1 goto :error

echo.
echo Done. Your Windows app is here:
echo dist\ASO Designer - by Alexander Apkarian\ASO Designer - by Alexander Apkarian.exe
echo.
pause
exit /b 0

:error
echo.
echo Build failed. Make sure Python 3.10 or newer is installed from python.org,
echo and that the "py" launcher is available.
echo.
pause
exit /b 1
