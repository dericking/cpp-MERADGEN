PROG     = meradgen10
EXE      = $(PROG).exe
TEST     = test.exe
SOURCES  = $(wildcard *.[f])
INC      = $(wildcard *.[inc])
OBJ      = fsir.o meradgen10.o run.o
OBJTEST  = fsir.o meradgen10.o test.o

COMPILER = f77

F77OPT  =-O2 -march=i486 -DLinux -fno-automatic -ffixed-line-length-none \
-fdollar-ok -fno-second-underscore -Iinclude

%.o:    %.f $(INC)
	$(COMPILER) $(F77OPT) -c $< -o $@
	
all:    $(EXE)

$(EXE):  $(OBJ) $(INC)
	$(COMPILER) $(F77OPT)  -o $(EXE) $(OBJ)

test:   $(OBJTEST) $(INC)
	$(COMPILER) $(F77OPT) -o $(TEST) $(OBJTEST)
clean:
	rm -f  *.o core

delete:
	rm -f $(EXE)

