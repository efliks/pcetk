PDYNAMO_ROOT = /home/mikolaj/local/opt/pDynamo-1.8.0-vanilla

INC1         = $(PDYNAMO_ROOT)/pCore-1.8.0/extensions/cinclude 
INC2         = ../cinclude
CFLAGS       = -fPIC -c -W -Wall -pedantic -I$(INC1) -I$(INC2)
CC           = gcc

StateVector:
	$(CC) $(CFLAGS) StateVector.c -o StateVector.o

clean:
	rm StateVector.o
