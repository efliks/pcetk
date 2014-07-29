PDYNAMO_ROOT = /home/mikolaj/local/opt/pDynamo-1.8.0

INC1   = $(PDYNAMO_ROOT)/pCore-1.8.0/extensions/cinclude 
INC2   = /usr/include/python2.7
INC3   = ../cinclude

CFLAGS = -fPIC -c -I$(INC1) -I$(INC2) -I$(INC3)
CC     = gcc

all: pyrex compile link install

pyrex:
	python pyrex_compile.py

compile: ContinuumElectrostatics.StateVector.c
	$(CC) $(CFLAGS) ContinuumElectrostatics.StateVector.c -o ContinuumElectrostatics.StateVector.o

link: ContinuumElectrostatics.StateVector.o ../csource/StateVector.o
	$(CC) -shared ContinuumElectrostatics.StateVector.o ../csource/StateVector.o -o StateVector.so

install: StateVector.so
	mv StateVector.so ../../ContinuumElectrostatics/

clean:
	rm ContinuumElectrostatics.StateVector.c ContinuumElectrostatics.StateVector.o
