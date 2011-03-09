
CALL "c:\Program Files (x86)\Microsoft Visual C++ Toolkit 2003\vcvars32.bat"
cl.exe /Oa /Ot /Og /O2  /I ..\include ..\src\wripaca.c
wripaca
pause
exit
