#  pyrexc ContinuumElectrostatics.StateVector.pyx 
#  
#  gcc -c -fPIC -I../cinclude/ -I/usr/include/python2.7/ -I/home/mikolaj/local/opt/pDynamo-1.8.0/pCore-1.8.0/extensions/cinclude/  ContinuumElectrostatics.StateVector.c 
#  
#  gcc -shared ContinuumElectrostatics.StateVector.o ../csource/StateVector.o -o StateVector.so

INC1   = /usr/include/python2.7
INC2   = /home/mikolaj/local/opt/pDynamo-1.8.0/pCore-1.8.0/extensions/cinclude 
INC3   = ../cinclude

CFLAGS = -fPIC -c -I$(INC1) -I$(INC2) -I$(INC3)
CC     = gcc

all: pyrex object link install

pyrex:
	pyrexc ContinuumElectrostatics.StateVector.pyx 

object: ContinuumElectrostatics.StateVector.c
	$(CC) $(CFLAGS) ContinuumElectrostatics.StateVector.c -o ContinuumElectrostatics.StateVector.o

link: ContinuumElectrostatics.StateVector.o ../csource/StateVector.o
	$(CC) -shared ContinuumElectrostatics.StateVector.o ../csource/StateVector.o -o StateVector.so

install: StateVector.so
	mv StateVector.so ../../ContinuumElectrostatics/

clean:
	rm ContinuumElectrostatics.StateVector.c ContinuumElectrostatics.StateVector.o
