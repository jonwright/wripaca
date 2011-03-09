

include make.inc

	
all : bin/testaffine$(EXE) bin/testwripaca$(EXE)

bin/testaffine$(EXE) : test/testaffine.c src/affine.c
	$(CC) $(CFLAGS) $^ $(MO)$@	

bin/testwripaca$(EXE) : test/testwripaca.c src/wripaca.c
	$(CC) $(CFLAGS) $^ $(MO)$@	

clean:
	cd bin && $(RM) testwripaca$(EXE) $(RM) testaffine$(EXE) && cd ..


