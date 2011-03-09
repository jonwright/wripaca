
CALL "c:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\bin\vcvars32.bat

cl.exe /nologo /MD /W3 /GS- /Ox /Og /Ot /fp:fast /openmp /I ..\include ..\src\wripaca.c
mt.exe -manifest wripaca.exe.manifest -outputresource:wripaca.exe;1
wripaca
pause
REM exit

