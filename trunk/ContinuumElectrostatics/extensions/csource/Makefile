PDYNAMO_PCORE = /home/mikolaj/local/opt/pDynamo-1.8.0/pCore-1.8.0
C_INCLUDE     = ../cinclude
CFLAGS        = -fPIC -c -W -Wall -pedantic -I$(PDYNAMO_PCORE)/extensions/cinclude -I$(C_INCLUDE)
CC            = gcc

StateVector:
	$(CC) $(CFLAGS) StateVector.c -o StateVector.o

clean:
	rm StateVector.o
