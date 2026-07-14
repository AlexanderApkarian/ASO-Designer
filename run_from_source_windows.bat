@echo off
setlocal

cd /d "%~dp0"

if not exist .venv_windows (
    echo Creating Windows run environment...
    py -3 -m venv .venv_windows
    if errorlevel 1 goto :error
)

call .venv_windows\Scripts\activate.bat
if errorlevel 1 goto :error

python -m pip install -r requirements.txt
if errorlevel 1 goto :error

python aso_designer.py
exit /b 0

:error
echo.
echo Could not run ASO Designer. Make sure Python 3.10 or newer is installed.
echo.
pause
exit /b 1
