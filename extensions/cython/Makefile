PDYNAMO_PCORE = /home/mikolaj/local/opt/pDynamo-1.8.0/pCore-1.8.0
PY_INCLUDE    = /usr/include/python2.7

# The options -fdata-sections and -ffunction-sections are for not linking the unused code
CFLAGS        = -fPIC -c -I$(PDYNAMO_PCORE)/extensions/cinclude -I$(PY_INCLUDE) -I../cinclude 
CC            = gcc


default: StateVector.so
	@echo "\n*** Use 'make clean_all' and then 'make' if you want to recompile Cython sources ***\n"

../csource/StateVector.a:
	+$(MAKE) -C ../csource

#-Wl,--gc-sections 
StateVector.so: ../csource/StateVector.a ContinuumElectrostatics.StateVector.o
	$(CC) -shared ContinuumElectrostatics.StateVector.o ../csource/StateVector.a -o StateVector.so

ContinuumElectrostatics.StateVector.o: ContinuumElectrostatics.StateVector.c
	$(CC) $(CFLAGS) ContinuumElectrostatics.StateVector.c -o ContinuumElectrostatics.StateVector.o

ContinuumElectrostatics.StateVector.c: ContinuumElectrostatics.StateVector.pyx ContinuumElectrostatics.StateVector.pxd
	python cython_compile.py $(PDYNAMO_PCORE)

install: StateVector.so
	mv StateVector.so ../../ContinuumElectrostatics/

clean:
	if [ -e ContinuumElectrostatics.StateVector.o ]; then rm ContinuumElectrostatics.StateVector.o ; fi
	if [ -e StateVector.so                        ]; then rm StateVector.so                        ; fi
	+$(MAKE) -C ../csource clean

clean_all: clean
	if [ -e ContinuumElectrostatics.StateVector.c ]; then rm ContinuumElectrostatics.StateVector.c ; fi
