PDYNAMO_PCORE = /home/mikolaj/local/opt/pDynamo-1.8.0/pCore-1.8.0
CFLAGS        = -fPIC -c -W -Wall -pedantic -I$(PDYNAMO_PCORE)/extensions/cinclude -I../cinclude
CC            = gcc

StateVector: StateVector.c ../cinclude/StateVector.h
	$(CC) $(CFLAGS) StateVector.c -o StateVector.o

clean:
	rm StateVector.o
