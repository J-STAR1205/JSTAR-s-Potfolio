@echo off
cd /d "%~dp0.."
call C:\Users\norma\miniconda3\Scripts\activate.bat
call conda activate qtcad
python -m notebook
