
COMPILER = msvc90_64
PLATFORM = win


INCDIR = include

include platform.inc.$(PLATFORM)
include make.inc.$(COMPILER)

# Root directory for building
R = builds$(S)$(PLATFORM)-$(COMPILER)

$(R):
	echo $(PLATFORM) $(COMPILER) $(R)
	$(MKDIR) $(R)$(S)bin
	$(MKDIR) $(R)$(S)test
	$(MKDIR) $(R)$(S)src

$(R)$(S)src$(S)affine$(OBJ) : src$(S)affine.c include$(S)affine.h
	$(CC) $(CFLAGS) $(MC) $< $(MOBJ)$@	

$(R)$(S)test$(S)testaffine$(OBJ) : test$(S)testaffine.c include$(S)affine.h
	$(CC) $(CFLAGS) $(MC) $< $(MOBJ)$@	

$(R)$(S)bin$(S)testaffine$(EXE) : $(R)$(S)src$(S)affine$(OBJ) $(R)$(S)test$(S)testaffine$(OBJ)
	$(CC) $(CFLAGS) $^ $(MEXE)$@	

$(R)$(S)src$(S)wripaca$(OBJ) : src$(S)wripaca.c include$(S)wripaca.h
	$(CC) $(CFLAGS) $(MC) $< $(MOBJ)$@

$(R)$(S)test$(S)testwripaca$(OBJ) : test$(S)testwripaca.c include$(S)wripaca.h
	$(CC) $(CFLAGS) $(MC) $< $(MOBJ)$@	

$(R)$(S)bin$(S)testwripaca$(EXE) : $(R)$(S)src$(S)wripaca$(OBJ) $(R)$(S)test$(S)testwripaca$(OBJ)
	$(CC) $(CFLAGS) $^ $(MEXE)$@	




all : $(R) $(R)$(S)bin$(S)testaffine$(EXE) $(R)$(S)bin$(S)testwripaca$(EXE)


