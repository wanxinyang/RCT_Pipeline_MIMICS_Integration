@echo off
REM Windows batch file to run MIMICS
REM Generated automatically

cd /d "D:\mimics\model\code"
echo Running MIMICS from %CD%
".\mimics1.5"
echo MIMICS execution completed with exit code %ERRORLEVEL%
