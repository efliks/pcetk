# PDYNAMO_PCORE comes from the upper-level make

# The options -fdata-sections and -ffunction-sections are for not linking the unused code
CFLAGS            = -O2 -fPIC -c -W -Wall -pedantic -I$(PDYNAMO_PCORE)/extensions/cinclude
CC                = gcc

OBJS =  Real2DArray.o MersenneTwister.o Transformation3Container.o IndexedSelection.o StatusHandler.o Integer1DArray.o \
 Transformation3.o AntisymmetricMatrix.o Matrix33.o PairList.o RegularGridOccupancy.o Integer.o \
 PolygonalSurface.o MarchingCubes.o CGLinearEquationSolver.o LAPACKEigenvalueSolver.o RandomNumberDistribution.o Real1DArray.o \
 Vector3.o JDEigenvalueSolver.o CubicSpline.o Memory.o Slice.o DoubleSymmetricMatrix.o \
 MemoryBlock.o Status.o RandomNumberGenerator.o NMREigenvalueSolver.o RegularGridDimensionData.o RealNDArray.o \
 SparseSymmetricMatrix.o BicubicSpline.o BlockStorage.o ArrayView.o Selection.o SliceOld.o \
 IntegerNDArray.o Coordinates3.o RegularGrid.o PairListGenerator.o Correlation.o Integer2DArray.o \
 List.o SelectionContainer.o Real.o Boolean1DArray.o Stack.o SymmetricMatrix.o \
 cblas_drotmg.o cblas_dtrmv.o cblas_dger.o cblas_drotg.o cblas_dsymv.o cblas_dtpsv.o \
 cblas_dtpmv.o cblas_globals.o cblas_dtrsm.o cblas_dgemv.o cblas_daxpy.o cblas_dscal.o \
 cblas_dsyrk.o cblas_dsyr2k.o cblas_dsyr2.o cblas_dgemm.o cblas_dcopy.o cblas_idamax.o \
 cblas_dspr2.o cblas_dtrsv.o cblas_xerbla.o cblas_dnrm2.o cblas_dspr.o cblas_ddot.o \
 cblas_dsyr.o cblas_dsymm.o cblas_drot.o cblas_dtbmv.o cblas_dswap.o cblas_dtbsv.o \
 cblas_dsbmv.o cblas_dspmv.o cblas_dtrmm.o cblas_dasum.o cblas_dgbmv.o cblas_drotm.o \
 d_atan.o derf_.o d_sin.o s_copy.o i_nint.o d_acos.o \
 d_sinh.o i_dim.o d_nint.o d_atn2.o i_len.o s_cmp.o \
 i_sign.o d_tanh.o i_abs.o d_lg10.o d_abs.o pow_di.o \
 pow_dd.o d_int.o d_sqrt.o d_tan.o d_asin.o d_log.o \
 d_exp.o pow_ii.o d_mod.o i_mod.o i_indx.o i_dnnt.o \
 d_cos.o derfc_.o d_prod.o d_dim.o d_sign.o d_cosh.o \
 dnrm2.o dgemm.o dspmv.o dsyr.o dsymv.o dtrmm.o \
 dtpmv.o ddot.o dtbsv.o dtrsv.o lsame.o drotmg.o \
 dsymm.o dcopy.o dscal.o idamax.o dtrmv.o dsdot.o \
 dspr2.o dgbmv.o dger.o xerbla.o daxpy.o dsyr2.o \
 drot.o dspr.o dswap.o dsyr2k.o dsbmv.o dtpsv.o \
 dsyrk.o drotm.o dtbmv.o dasum.o dtrsm.o drotg.o \
 dgemv.o dlamch.o dlansp.o iladlr.o ilaenv.o dlaebz.o \
 dlamrg.o dlarfb.o dlarnv.o dlaed8.o dlaed2.o dlansy.o \
 dlanst.o dstebz.o dspevx.o dsterf.o ieeeck.o dlaswp.o \
 dlatrd.o dstein.o dlaed4.o dgtsv.o dlassq.o lsame.o \
 dlaeda.o dsytf2.o dlaset.o dsteqr.o dlaed1.o dlaed9.o \
 dsptrs.o dsptrf.o dlagts.o dlascl.o iladlc.o dlagtf.o \
 dsprfs.o dlaev2.o dlasr.o dlaed7.o dlacpy.o iparmq.o \
 dlartg.o dopgtr.o xerbla.o dspevd.o disnan.o dlarfg.o \
 dlarft.o dlasrt.o dsytd2.o dlasyf.o dorgtr.o dopmtr.o \
 dsptrd.o dstedc.o dsytrf.o dspsvx.o dgetrf.o dsytrd.o \
 dlaed6.o dorg2r.o dlaed5.o dorg2l.o dlaed0.o dlapy2.o \
 dsytrs.o dlae2.o dlacn2.o dlaed3.o dlarf.o dlaisnan.o \
 dlaruv.o dgetf2.o dorgqr.o dsyev.o dspcon.o dorgql.o

all: $(OBJS)
	ar rcs libpcore.a *.o

clean:
	if [ -e libpcore.a ] ; then rm libpcore.a ; fi
	rm *.o

Real2DArray.o: $(PDYNAMO_PCORE)/extensions/csource/Real2DArray.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/Real2DArray.c -o Real2DArray.o

MersenneTwister.o: $(PDYNAMO_PCORE)/extensions/csource/MersenneTwister.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/MersenneTwister.c -o MersenneTwister.o

Transformation3Container.o: $(PDYNAMO_PCORE)/extensions/csource/Transformation3Container.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/Transformation3Container.c -o Transformation3Container.o

IndexedSelection.o: $(PDYNAMO_PCORE)/extensions/csource/IndexedSelection.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/IndexedSelection.c -o IndexedSelection.o

StatusHandler.o: $(PDYNAMO_PCORE)/extensions/csource/StatusHandler.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/StatusHandler.c -o StatusHandler.o

Integer1DArray.o: $(PDYNAMO_PCORE)/extensions/csource/Integer1DArray.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/Integer1DArray.c -o Integer1DArray.o

Transformation3.o: $(PDYNAMO_PCORE)/extensions/csource/Transformation3.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/Transformation3.c -o Transformation3.o

AntisymmetricMatrix.o: $(PDYNAMO_PCORE)/extensions/csource/AntisymmetricMatrix.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/AntisymmetricMatrix.c -o AntisymmetricMatrix.o

Matrix33.o: $(PDYNAMO_PCORE)/extensions/csource/Matrix33.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/Matrix33.c -o Matrix33.o

PairList.o: $(PDYNAMO_PCORE)/extensions/csource/PairList.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/PairList.c -o PairList.o

RegularGridOccupancy.o: $(PDYNAMO_PCORE)/extensions/csource/RegularGridOccupancy.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/RegularGridOccupancy.c -o RegularGridOccupancy.o

Integer.o: $(PDYNAMO_PCORE)/extensions/csource/Integer.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/Integer.c -o Integer.o

PolygonalSurface.o: $(PDYNAMO_PCORE)/extensions/csource/PolygonalSurface.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/PolygonalSurface.c -o PolygonalSurface.o

MarchingCubes.o: $(PDYNAMO_PCORE)/extensions/csource/MarchingCubes.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/MarchingCubes.c -o MarchingCubes.o

CGLinearEquationSolver.o: $(PDYNAMO_PCORE)/extensions/csource/CGLinearEquationSolver.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/CGLinearEquationSolver.c -o CGLinearEquationSolver.o

LAPACKEigenvalueSolver.o: $(PDYNAMO_PCORE)/extensions/csource/LAPACKEigenvalueSolver.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/LAPACKEigenvalueSolver.c -o LAPACKEigenvalueSolver.o

RandomNumberDistribution.o: $(PDYNAMO_PCORE)/extensions/csource/RandomNumberDistribution.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/RandomNumberDistribution.c -o RandomNumberDistribution.o

Real1DArray.o: $(PDYNAMO_PCORE)/extensions/csource/Real1DArray.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/Real1DArray.c -o Real1DArray.o

Vector3.o: $(PDYNAMO_PCORE)/extensions/csource/Vector3.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/Vector3.c -o Vector3.o

JDEigenvalueSolver.o: $(PDYNAMO_PCORE)/extensions/csource/JDEigenvalueSolver.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/JDEigenvalueSolver.c -o JDEigenvalueSolver.o

CubicSpline.o: $(PDYNAMO_PCORE)/extensions/csource/CubicSpline.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/CubicSpline.c -o CubicSpline.o

Memory.o: $(PDYNAMO_PCORE)/extensions/csource/Memory.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/Memory.c -o Memory.o

Slice.o: $(PDYNAMO_PCORE)/extensions/csource/Slice.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/Slice.c -o Slice.o

DoubleSymmetricMatrix.o: $(PDYNAMO_PCORE)/extensions/csource/DoubleSymmetricMatrix.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/DoubleSymmetricMatrix.c -o DoubleSymmetricMatrix.o

MemoryBlock.o: $(PDYNAMO_PCORE)/extensions/csource/MemoryBlock.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/MemoryBlock.c -o MemoryBlock.o

Status.o: $(PDYNAMO_PCORE)/extensions/csource/Status.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/Status.c -o Status.o

RandomNumberGenerator.o: $(PDYNAMO_PCORE)/extensions/csource/RandomNumberGenerator.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/RandomNumberGenerator.c -o RandomNumberGenerator.o

NMREigenvalueSolver.o: $(PDYNAMO_PCORE)/extensions/csource/NMREigenvalueSolver.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/NMREigenvalueSolver.c -o NMREigenvalueSolver.o

RegularGridDimensionData.o: $(PDYNAMO_PCORE)/extensions/csource/RegularGridDimensionData.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/RegularGridDimensionData.c -o RegularGridDimensionData.o

RealNDArray.o: $(PDYNAMO_PCORE)/extensions/csource/RealNDArray.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/RealNDArray.c -o RealNDArray.o

SparseSymmetricMatrix.o: $(PDYNAMO_PCORE)/extensions/csource/SparseSymmetricMatrix.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/SparseSymmetricMatrix.c -o SparseSymmetricMatrix.o

BicubicSpline.o: $(PDYNAMO_PCORE)/extensions/csource/BicubicSpline.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/BicubicSpline.c -o BicubicSpline.o

BlockStorage.o: $(PDYNAMO_PCORE)/extensions/csource/BlockStorage.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/BlockStorage.c -o BlockStorage.o

ArrayView.o: $(PDYNAMO_PCORE)/extensions/csource/ArrayView.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/ArrayView.c -o ArrayView.o

Selection.o: $(PDYNAMO_PCORE)/extensions/csource/Selection.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/Selection.c -o Selection.o

SliceOld.o: $(PDYNAMO_PCORE)/extensions/csource/SliceOld.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/SliceOld.c -o SliceOld.o

IntegerNDArray.o: $(PDYNAMO_PCORE)/extensions/csource/IntegerNDArray.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/IntegerNDArray.c -o IntegerNDArray.o

Coordinates3.o: $(PDYNAMO_PCORE)/extensions/csource/Coordinates3.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/Coordinates3.c -o Coordinates3.o

RegularGrid.o: $(PDYNAMO_PCORE)/extensions/csource/RegularGrid.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/RegularGrid.c -o RegularGrid.o

PairListGenerator.o: $(PDYNAMO_PCORE)/extensions/csource/PairListGenerator.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/PairListGenerator.c -o PairListGenerator.o

Correlation.o: $(PDYNAMO_PCORE)/extensions/csource/Correlation.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/Correlation.c -o Correlation.o

Integer2DArray.o: $(PDYNAMO_PCORE)/extensions/csource/Integer2DArray.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/Integer2DArray.c -o Integer2DArray.o

List.o: $(PDYNAMO_PCORE)/extensions/csource/List.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/List.c -o List.o

SelectionContainer.o: $(PDYNAMO_PCORE)/extensions/csource/SelectionContainer.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/SelectionContainer.c -o SelectionContainer.o

Real.o: $(PDYNAMO_PCORE)/extensions/csource/Real.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/Real.c -o Real.o

Boolean1DArray.o: $(PDYNAMO_PCORE)/extensions/csource/Boolean1DArray.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/Boolean1DArray.c -o Boolean1DArray.o

Stack.o: $(PDYNAMO_PCORE)/extensions/csource/Stack.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/Stack.c -o Stack.o

SymmetricMatrix.o: $(PDYNAMO_PCORE)/extensions/csource/SymmetricMatrix.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/csource/SymmetricMatrix.c -o SymmetricMatrix.o

cblas_drotmg.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_drotmg.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_drotmg.c -o cblas_drotmg.o

cblas_dtrmv.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dtrmv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dtrmv.c -o cblas_dtrmv.o

cblas_dger.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dger.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dger.c -o cblas_dger.o

cblas_drotg.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_drotg.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_drotg.c -o cblas_drotg.o

cblas_dsymv.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dsymv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dsymv.c -o cblas_dsymv.o

cblas_dtpsv.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dtpsv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dtpsv.c -o cblas_dtpsv.o

cblas_dtpmv.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dtpmv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dtpmv.c -o cblas_dtpmv.o

cblas_globals.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_globals.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_globals.c -o cblas_globals.o

cblas_dtrsm.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dtrsm.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dtrsm.c -o cblas_dtrsm.o

cblas_dgemv.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dgemv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dgemv.c -o cblas_dgemv.o

cblas_daxpy.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_daxpy.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_daxpy.c -o cblas_daxpy.o

cblas_dscal.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dscal.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dscal.c -o cblas_dscal.o

cblas_dsyrk.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dsyrk.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dsyrk.c -o cblas_dsyrk.o

cblas_dsyr2k.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dsyr2k.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dsyr2k.c -o cblas_dsyr2k.o

cblas_dsyr2.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dsyr2.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dsyr2.c -o cblas_dsyr2.o

cblas_dgemm.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dgemm.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dgemm.c -o cblas_dgemm.o

cblas_dcopy.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dcopy.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dcopy.c -o cblas_dcopy.o

cblas_idamax.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_idamax.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_idamax.c -o cblas_idamax.o

cblas_dspr2.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dspr2.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dspr2.c -o cblas_dspr2.o

cblas_dtrsv.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dtrsv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dtrsv.c -o cblas_dtrsv.o

cblas_xerbla.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_xerbla.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_xerbla.c -o cblas_xerbla.o

cblas_dnrm2.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dnrm2.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dnrm2.c -o cblas_dnrm2.o

cblas_dspr.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dspr.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dspr.c -o cblas_dspr.o

cblas_ddot.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_ddot.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_ddot.c -o cblas_ddot.o

cblas_dsyr.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dsyr.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dsyr.c -o cblas_dsyr.o

cblas_dsymm.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dsymm.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dsymm.c -o cblas_dsymm.o

cblas_drot.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_drot.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_drot.c -o cblas_drot.o

cblas_dtbmv.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dtbmv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dtbmv.c -o cblas_dtbmv.o

cblas_dswap.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dswap.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dswap.c -o cblas_dswap.o

cblas_dtbsv.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dtbsv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dtbsv.c -o cblas_dtbsv.o

cblas_dsbmv.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dsbmv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dsbmv.c -o cblas_dsbmv.o

cblas_dspmv.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dspmv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dspmv.c -o cblas_dspmv.o

cblas_dtrmm.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dtrmm.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dtrmm.c -o cblas_dtrmm.o

cblas_dasum.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dasum.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dasum.c -o cblas_dasum.o

cblas_dgbmv.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_dgbmv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_dgbmv.c -o cblas_dgbmv.o

cblas_drotm.o: $(PDYNAMO_PCORE)/extensions/cblas/cblas_drotm.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/cblas/cblas_drotm.c -o cblas_drotm.o

d_atan.o: $(PDYNAMO_PCORE)/extensions/f2c/d_atan.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/d_atan.c -o d_atan.o

derf_.o: $(PDYNAMO_PCORE)/extensions/f2c/derf_.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/derf_.c -o derf_.o

d_sin.o: $(PDYNAMO_PCORE)/extensions/f2c/d_sin.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/d_sin.c -o d_sin.o

s_copy.o: $(PDYNAMO_PCORE)/extensions/f2c/s_copy.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/s_copy.c -o s_copy.o

i_nint.o: $(PDYNAMO_PCORE)/extensions/f2c/i_nint.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/i_nint.c -o i_nint.o

d_acos.o: $(PDYNAMO_PCORE)/extensions/f2c/d_acos.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/d_acos.c -o d_acos.o

d_sinh.o: $(PDYNAMO_PCORE)/extensions/f2c/d_sinh.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/d_sinh.c -o d_sinh.o

i_dim.o: $(PDYNAMO_PCORE)/extensions/f2c/i_dim.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/i_dim.c -o i_dim.o

d_nint.o: $(PDYNAMO_PCORE)/extensions/f2c/d_nint.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/d_nint.c -o d_nint.o

d_atn2.o: $(PDYNAMO_PCORE)/extensions/f2c/d_atn2.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/d_atn2.c -o d_atn2.o

i_len.o: $(PDYNAMO_PCORE)/extensions/f2c/i_len.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/i_len.c -o i_len.o

s_cmp.o: $(PDYNAMO_PCORE)/extensions/f2c/s_cmp.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/s_cmp.c -o s_cmp.o

i_sign.o: $(PDYNAMO_PCORE)/extensions/f2c/i_sign.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/i_sign.c -o i_sign.o

d_tanh.o: $(PDYNAMO_PCORE)/extensions/f2c/d_tanh.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/d_tanh.c -o d_tanh.o

i_abs.o: $(PDYNAMO_PCORE)/extensions/f2c/i_abs.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/i_abs.c -o i_abs.o

d_lg10.o: $(PDYNAMO_PCORE)/extensions/f2c/d_lg10.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/d_lg10.c -o d_lg10.o

d_abs.o: $(PDYNAMO_PCORE)/extensions/f2c/d_abs.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/d_abs.c -o d_abs.o

pow_di.o: $(PDYNAMO_PCORE)/extensions/f2c/pow_di.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/pow_di.c -o pow_di.o

pow_dd.o: $(PDYNAMO_PCORE)/extensions/f2c/pow_dd.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/pow_dd.c -o pow_dd.o

d_int.o: $(PDYNAMO_PCORE)/extensions/f2c/d_int.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/d_int.c -o d_int.o

d_sqrt.o: $(PDYNAMO_PCORE)/extensions/f2c/d_sqrt.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/d_sqrt.c -o d_sqrt.o

d_tan.o: $(PDYNAMO_PCORE)/extensions/f2c/d_tan.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/d_tan.c -o d_tan.o

d_asin.o: $(PDYNAMO_PCORE)/extensions/f2c/d_asin.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/d_asin.c -o d_asin.o

d_log.o: $(PDYNAMO_PCORE)/extensions/f2c/d_log.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/d_log.c -o d_log.o

d_exp.o: $(PDYNAMO_PCORE)/extensions/f2c/d_exp.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/d_exp.c -o d_exp.o

pow_ii.o: $(PDYNAMO_PCORE)/extensions/f2c/pow_ii.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/pow_ii.c -o pow_ii.o

d_mod.o: $(PDYNAMO_PCORE)/extensions/f2c/d_mod.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/d_mod.c -o d_mod.o

i_mod.o: $(PDYNAMO_PCORE)/extensions/f2c/i_mod.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/i_mod.c -o i_mod.o

i_indx.o: $(PDYNAMO_PCORE)/extensions/f2c/i_indx.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/i_indx.c -o i_indx.o

i_dnnt.o: $(PDYNAMO_PCORE)/extensions/f2c/i_dnnt.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/i_dnnt.c -o i_dnnt.o

d_cos.o: $(PDYNAMO_PCORE)/extensions/f2c/d_cos.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/d_cos.c -o d_cos.o

derfc_.o: $(PDYNAMO_PCORE)/extensions/f2c/derfc_.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/derfc_.c -o derfc_.o

d_prod.o: $(PDYNAMO_PCORE)/extensions/f2c/d_prod.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/d_prod.c -o d_prod.o

d_dim.o: $(PDYNAMO_PCORE)/extensions/f2c/d_dim.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/d_dim.c -o d_dim.o

d_sign.o: $(PDYNAMO_PCORE)/extensions/f2c/d_sign.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/d_sign.c -o d_sign.o

d_cosh.o: $(PDYNAMO_PCORE)/extensions/f2c/d_cosh.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2c/d_cosh.c -o d_cosh.o

dnrm2.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dnrm2.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dnrm2.c -o dnrm2.o

dgemm.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dgemm.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dgemm.c -o dgemm.o

dspmv.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dspmv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dspmv.c -o dspmv.o

dsyr.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dsyr.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dsyr.c -o dsyr.o

dsymv.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dsymv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dsymv.c -o dsymv.o

dtrmm.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dtrmm.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dtrmm.c -o dtrmm.o

dtpmv.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dtpmv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dtpmv.c -o dtpmv.o

ddot.o: $(PDYNAMO_PCORE)/extensions/f2cblas/ddot.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/ddot.c -o ddot.o

dtbsv.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dtbsv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dtbsv.c -o dtbsv.o

dtrsv.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dtrsv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dtrsv.c -o dtrsv.o

lsame.o: $(PDYNAMO_PCORE)/extensions/f2cblas/lsame.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/lsame.c -o lsame.o

drotmg.o: $(PDYNAMO_PCORE)/extensions/f2cblas/drotmg.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/drotmg.c -o drotmg.o

dsymm.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dsymm.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dsymm.c -o dsymm.o

dcopy.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dcopy.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dcopy.c -o dcopy.o

dscal.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dscal.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dscal.c -o dscal.o

idamax.o: $(PDYNAMO_PCORE)/extensions/f2cblas/idamax.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/idamax.c -o idamax.o

dtrmv.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dtrmv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dtrmv.c -o dtrmv.o

dsdot.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dsdot.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dsdot.c -o dsdot.o

dspr2.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dspr2.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dspr2.c -o dspr2.o

dgbmv.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dgbmv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dgbmv.c -o dgbmv.o

dger.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dger.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dger.c -o dger.o

xerbla.o: $(PDYNAMO_PCORE)/extensions/f2cblas/xerbla.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/xerbla.c -o xerbla.o

daxpy.o: $(PDYNAMO_PCORE)/extensions/f2cblas/daxpy.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/daxpy.c -o daxpy.o

dsyr2.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dsyr2.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dsyr2.c -o dsyr2.o

drot.o: $(PDYNAMO_PCORE)/extensions/f2cblas/drot.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/drot.c -o drot.o

dspr.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dspr.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dspr.c -o dspr.o

dswap.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dswap.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dswap.c -o dswap.o

dsyr2k.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dsyr2k.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dsyr2k.c -o dsyr2k.o

dsbmv.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dsbmv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dsbmv.c -o dsbmv.o

dtpsv.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dtpsv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dtpsv.c -o dtpsv.o

dsyrk.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dsyrk.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dsyrk.c -o dsyrk.o

drotm.o: $(PDYNAMO_PCORE)/extensions/f2cblas/drotm.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/drotm.c -o drotm.o

dtbmv.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dtbmv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dtbmv.c -o dtbmv.o

dasum.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dasum.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dasum.c -o dasum.o

dtrsm.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dtrsm.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dtrsm.c -o dtrsm.o

drotg.o: $(PDYNAMO_PCORE)/extensions/f2cblas/drotg.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/drotg.c -o drotg.o

dgemv.o: $(PDYNAMO_PCORE)/extensions/f2cblas/dgemv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cblas/dgemv.c -o dgemv.o

dlamch.o: $(PDYNAMO_PCORE)/extensions/f2cdlamch/dlamch.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2cdlamch/dlamch.c -o dlamch.o

dlansp.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlansp.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlansp.c -o dlansp.o

iladlr.o: $(PDYNAMO_PCORE)/extensions/f2clapack/iladlr.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/iladlr.c -o iladlr.o

ilaenv.o: $(PDYNAMO_PCORE)/extensions/f2clapack/ilaenv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/ilaenv.c -o ilaenv.o

dlaebz.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlaebz.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlaebz.c -o dlaebz.o

dlamrg.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlamrg.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlamrg.c -o dlamrg.o

dlarfb.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlarfb.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlarfb.c -o dlarfb.o

dlarnv.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlarnv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlarnv.c -o dlarnv.o

dlaed8.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlaed8.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlaed8.c -o dlaed8.o

dlaed2.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlaed2.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlaed2.c -o dlaed2.o

dlansy.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlansy.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlansy.c -o dlansy.o

dlanst.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlanst.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlanst.c -o dlanst.o

dstebz.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dstebz.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dstebz.c -o dstebz.o

dspevx.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dspevx.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dspevx.c -o dspevx.o

dsterf.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dsterf.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dsterf.c -o dsterf.o

ieeeck.o: $(PDYNAMO_PCORE)/extensions/f2clapack/ieeeck.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/ieeeck.c -o ieeeck.o

dlaswp.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlaswp.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlaswp.c -o dlaswp.o

dlatrd.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlatrd.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlatrd.c -o dlatrd.o

dstein.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dstein.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dstein.c -o dstein.o

dlaed4.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlaed4.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlaed4.c -o dlaed4.o

dgtsv.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dgtsv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dgtsv.c -o dgtsv.o

dlassq.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlassq.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlassq.c -o dlassq.o

lsame.o: $(PDYNAMO_PCORE)/extensions/f2clapack/lsame.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/lsame.c -o lsame.o

dlaeda.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlaeda.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlaeda.c -o dlaeda.o

dsytf2.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dsytf2.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dsytf2.c -o dsytf2.o

dlaset.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlaset.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlaset.c -o dlaset.o

dsteqr.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dsteqr.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dsteqr.c -o dsteqr.o

dlaed1.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlaed1.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlaed1.c -o dlaed1.o

dlaed9.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlaed9.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlaed9.c -o dlaed9.o

dsptrs.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dsptrs.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dsptrs.c -o dsptrs.o

dsptrf.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dsptrf.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dsptrf.c -o dsptrf.o

dlagts.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlagts.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlagts.c -o dlagts.o

dlascl.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlascl.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlascl.c -o dlascl.o

iladlc.o: $(PDYNAMO_PCORE)/extensions/f2clapack/iladlc.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/iladlc.c -o iladlc.o

dlagtf.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlagtf.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlagtf.c -o dlagtf.o

dsprfs.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dsprfs.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dsprfs.c -o dsprfs.o

dlaev2.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlaev2.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlaev2.c -o dlaev2.o

dlasr.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlasr.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlasr.c -o dlasr.o

dlaed7.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlaed7.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlaed7.c -o dlaed7.o

dlacpy.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlacpy.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlacpy.c -o dlacpy.o

iparmq.o: $(PDYNAMO_PCORE)/extensions/f2clapack/iparmq.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/iparmq.c -o iparmq.o

dlartg.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlartg.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlartg.c -o dlartg.o

dopgtr.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dopgtr.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dopgtr.c -o dopgtr.o

xerbla.o: $(PDYNAMO_PCORE)/extensions/f2clapack/xerbla.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/xerbla.c -o xerbla.o

dspevd.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dspevd.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dspevd.c -o dspevd.o

disnan.o: $(PDYNAMO_PCORE)/extensions/f2clapack/disnan.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/disnan.c -o disnan.o

dlarfg.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlarfg.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlarfg.c -o dlarfg.o

dlarft.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlarft.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlarft.c -o dlarft.o

dlasrt.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlasrt.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlasrt.c -o dlasrt.o

dsytd2.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dsytd2.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dsytd2.c -o dsytd2.o

dlasyf.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlasyf.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlasyf.c -o dlasyf.o

dorgtr.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dorgtr.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dorgtr.c -o dorgtr.o

dopmtr.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dopmtr.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dopmtr.c -o dopmtr.o

dsptrd.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dsptrd.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dsptrd.c -o dsptrd.o

dstedc.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dstedc.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dstedc.c -o dstedc.o

dsytrf.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dsytrf.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dsytrf.c -o dsytrf.o

dspsvx.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dspsvx.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dspsvx.c -o dspsvx.o

dgetrf.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dgetrf.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dgetrf.c -o dgetrf.o

dsytrd.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dsytrd.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dsytrd.c -o dsytrd.o

dlaed6.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlaed6.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlaed6.c -o dlaed6.o

dorg2r.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dorg2r.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dorg2r.c -o dorg2r.o

dlaed5.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlaed5.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlaed5.c -o dlaed5.o

dorg2l.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dorg2l.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dorg2l.c -o dorg2l.o

dlaed0.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlaed0.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlaed0.c -o dlaed0.o

dlapy2.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlapy2.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlapy2.c -o dlapy2.o

dsytrs.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dsytrs.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dsytrs.c -o dsytrs.o

dlae2.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlae2.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlae2.c -o dlae2.o

dlacn2.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlacn2.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlacn2.c -o dlacn2.o

dlaed3.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlaed3.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlaed3.c -o dlaed3.o

dlarf.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlarf.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlarf.c -o dlarf.o

dlaisnan.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlaisnan.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlaisnan.c -o dlaisnan.o

dlaruv.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dlaruv.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dlaruv.c -o dlaruv.o

dgetf2.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dgetf2.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dgetf2.c -o dgetf2.o

dorgqr.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dorgqr.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dorgqr.c -o dorgqr.o

dsyev.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dsyev.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dsyev.c -o dsyev.o

dspcon.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dspcon.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dspcon.c -o dspcon.o

dorgql.o: $(PDYNAMO_PCORE)/extensions/f2clapack/dorgql.c
	$(CC) $(CFLAGS) $(PDYNAMO_PCORE)/extensions/f2clapack/dorgql.c -o dorgql.o

