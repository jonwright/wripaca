
TARGET = msvc90_64

include make.inc.$(TARGET)

	
all : $(TARGET)  $(TARGET)/bin \
       	$(TARGET)/bin/testaffine$(EXE) \
	$(TARGET)/bin/testwripaca$(EXE)

$(TARGET): 
	mkdir $(TARGET)

$(TARGET)/bin:
	cd  $(TARGET) && mkdir bin

$(TARGET)/bin/testaffine$(EXE) : test/testaffine.c src/affine.c
	$(CC) $(CFLAGS) $^ $(MO)$@	

$(TARGET)/bin/testwripaca$(EXE) : test/testwripaca.c src/wripaca.c
	$(CC) $(CFLAGS) $^ $(MO)$@	

clean:
	cd $(TARGET)/bin && $(RM) testwripaca$(EXE) $(RM) testaffine$(EXE) && cd .. 



