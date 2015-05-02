# PDYNAMO_PCORE comes from the Makefile in the cython directory

# The options -fdata-sections and -ffunction-sections are for not linking the unused code
CFLAGS        = -O2 -fPIC -c -W -Wall -pedantic -I$(PDYNAMO_PCORE)/extensions/cinclude -I../cinclude
CC            = gcc

default: MCModelDefault.o EnergyModel.o StateVector.o lib/libpcore.a

MCModelDefault.o: MCModelDefault.c ../cinclude/MCModelDefault.h
	$(CC) $(CFLAGS) MCModelDefault.c -o MCModelDefault.o

EnergyModel.o: EnergyModel.c ../cinclude/EnergyModel.h
	$(CC) $(CFLAGS) EnergyModel.c -o EnergyModel.o

StateVector.o: StateVector.c ../cinclude/StateVector.h
	$(CC) $(CFLAGS) StateVector.c -o StateVector.o

lib/libpcore.a:
	+$(MAKE) -C lib

clean:
	if [ -e MCModelDefault.o ] ; then rm MCModelDefault.o ; fi
	if [ -e EnergyModel.o    ] ; then rm EnergyModel.o    ; fi
	if [ -e StateVector.o    ] ; then rm StateVector.o    ; fi

clean_all: clean
	+$(MAKE) -C lib clean
