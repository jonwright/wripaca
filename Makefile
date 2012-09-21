
#COMPILER = msvc90_64
COMPILER = gcc4.5
PLATFORM = win


include platform.inc.$(PLATFORM)
include make.inc.$(COMPILER)
INCDIR = include


# Root directory for building
R = builds$(S)$(PLATFORM)-$(COMPILER)

$(R):
	echo $(PLATFORM) $(COMPILER) $(R)
	$(MKDIR) $(R)$(S)bin
	$(MKDIR) $(R)$(S)test
	$(MKDIR) $(R)$(S)src

$(R)$(S)src$(S)affine$(OBJ) : src$(S)affine.c include$(S)affine.h
	$(CC) $(CFLAGS) -I $(INCDIR) $(MC) $< $(MOBJ)$@	

$(R)$(S)test$(S)testaffine$(OBJ) : test$(S)testaffine.c include$(S)affine.h
	$(CC) $(CFLAGS) -I $(INCDIR) $(MC) $< $(MOBJ)$@	

$(R)$(S)bin$(S)testaffine$(EXE) : $(R)$(S)src$(S)affine$(OBJ) $(R)$(S)test$(S)testaffine$(OBJ)
	$(CC) $(CFLAGS) $^ $(MEXE)$@	

$(R)$(S)src$(S)wripaca$(OBJ) : src$(S)wripaca.c include$(S)wripaca.h
	$(CC) $(CFLAGS) -I $(INCDIR) $(MC) $< $(MOBJ)$@

$(R)$(S)test$(S)testwripaca$(OBJ) : test$(S)testwripaca.c include$(S)wripaca.h
	$(CC) $(CFLAGS) -I $(INCDIR) $(OPENMP) $(MC) $< $(MOBJ)$@	

$(R)$(S)bin$(S)testwripaca$(EXE) : $(R)$(S)src$(S)wripaca$(OBJ) $(R)$(S)test$(S)testwripaca$(OBJ)
	$(CC) $(CFLAGS) $(OPENMP) $^ $(MEXE)$@


$(R)$(S)src$(S)add$(OBJ) : src$(S)add.c 
	echo $(OCLINC)
	$(CC) $(CFLAGS) -I $(INCDIR) -I $(OCLINC) $(MC) $^ $(MOBJ)$@


$(R)$(S)bin$(S)add$(EXE) : $(R)$(S)src$(S)add$(OBJ) 
	$(CC) $(CFLAGS) $^ $(MEXE)$@ $(OCLLIB) $(ML)$(OCLLIBDIR) 



all : $(R) $(R)$(S)bin$(S)testaffine$(EXE) $(R)$(S)bin$(S)testwripaca$(EXE) $(R)$(S)bin$(S)add$(EXE)


clean :
	cd $(R)$(S)bin && $(RM) *.*
	cd $(R)$(S)test && $(RM) *.*
	cd $(R)$(S)src && $(RM) *.*


