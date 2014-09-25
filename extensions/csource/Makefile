# PDYNAMO_PCORE comes from the Makefile in the cython directory

# The options -fdata-sections and -ffunction-sections are for not linking the unused code
CFLAGS        = -fPIC -c -W -Wall -pedantic -I$(PDYNAMO_PCORE)/extensions/cinclude -I../cinclude
CC            = gcc

StateVector.a: StateVector.o lib/libpcore.a
	ar rcs StateVector.a StateVector.o lib/libpcore.a

StateVector.o: StateVector.c ../cinclude/StateVector.h
	$(CC) $(CFLAGS) StateVector.c -o StateVector.o

lib/libpcore.a:
	+$(MAKE) -C lib

clean:
	if [ -e StateVector.o ] ; then rm StateVector.o ; fi
	if [ -e StateVector.a ] ; then rm StateVector.a ; fi

clean_all: clean
	+$(MAKE) -C lib clean
