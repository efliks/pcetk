INC1   = /home/mikolaj/local/opt/pDynamo-1.8.0/pCore-1.8.0/extensions/cinclude 
INC2   = ../cinclude

CFLAGS = -fPIC -c -I$(INC1) -I$(INC2)
CC     = gcc

StateVector:
	$(CC) $(CFLAGS) StateVector.c -o StateVector.o

clean:
	rm StateVector.o
