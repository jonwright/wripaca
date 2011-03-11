
src: python/wripaca.py wripaca_wrap.c ../include/wripaca.h

wripaca_wrap.c python/wripaca.py: wripaca.i ../include/wripaca.h
	swig -python -outdir python wripaca.i
 

