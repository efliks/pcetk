PDYNAMO_PCORE = /home/mikolaj/local/opt/pDynamo-1.8.0/pCore-1.8.0
PY_INCLUDE    = /usr/include/python2.7
CFLAGS        = -fPIC -c -I$(PDYNAMO_PCORE)/extensions/cinclude -I$(PY_INCLUDE) -I../cinclude 
CC            = gcc


default: compile link
	@echo "\n*** Use 'make all' if you want to recompile Pyrex sources ***\n"

all: pyrex compile link


pyrex: ContinuumElectrostatics.StateVector.pyx ContinuumElectrostatics.StateVector.pxd
	python pyrex_compile.py $(PDYNAMO_PCORE)

compile: ContinuumElectrostatics.StateVector.c ../cinclude/StateVector.h
	$(CC) $(CFLAGS) ContinuumElectrostatics.StateVector.c -o ContinuumElectrostatics.StateVector.o

link: ContinuumElectrostatics.StateVector.o ../csource/StateVector.o
	$(CC) -shared ContinuumElectrostatics.StateVector.o ../csource/StateVector.o -o StateVector.so

install: StateVector.so
	mv StateVector.so ../../ContinuumElectrostatics/

clean:
	rm ContinuumElectrostatics.StateVector.o
