# PDYNAMO_PCORE comes from the Makefile in the pyrex directory

CFLAGS        = -fPIC -c -W -Wall -pedantic -I$(PDYNAMO_PCORE)/extensions/cinclude -I../cinclude
CC            = gcc

StateVector.o: StateVector.c ../cinclude/StateVector.h
	$(CC) $(CFLAGS) StateVector.c -o StateVector.o

clean:
	if [ -e StateVector.o ] ; then rm StateVector.o ; fi
