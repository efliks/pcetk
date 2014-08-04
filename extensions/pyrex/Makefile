PDYNAMO_PCORE = /home/mikolaj/local/opt/pDynamo-1.8.0/pCore-1.8.0
PY_INCLUDE    = /usr/include/python2.7
C_INCLUDE     = ../cinclude
CFLAGS        = -fPIC -c -I$(PDYNAMO_PCORE)/extensions/cinclude -I$(C_INCLUDE) -I$(PY_INCLUDE)
CC            = gcc

all: pyrex compile link install

pyrex:
	python pyrex_compile.py $(PDYNAMO_PCORE)

compile: ContinuumElectrostatics.StateVector.c
	$(CC) $(CFLAGS) ContinuumElectrostatics.StateVector.c -o ContinuumElectrostatics.StateVector.o

link: ContinuumElectrostatics.StateVector.o ../csource/StateVector.o
	$(CC) -shared ContinuumElectrostatics.StateVector.o ../csource/StateVector.o -o StateVector.so

install: StateVector.so
	mv StateVector.so ../../ContinuumElectrostatics/

clean:
	rm ContinuumElectrostatics.StateVector.c ContinuumElectrostatics.StateVector.o
