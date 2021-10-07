#-------------------------------------------------------------------------------
# defaults
#-------------------------------------------------------------------------------
MAKE = make
AR = ar


all:  mstm-scan

mstm-scan: mpidefs-serial-scan.o intrinsics-scan.o modules-scan.o main-scan.o
	gfortran mpidefs-serial-scan.o intrinsics-scan.o modules-scan.o main-scan.o -o mstm-scan.exe

mpidefs-serial-scan.o:
	gfortran -c mpidefs-serial-scan.f90

intrinsics-scan.o:
	gfortran -c intrinsics-scan.f90


modules-scan.o:
	gfortran -c modules-scan.f90


main-scan.o:
	gfortran -c main-scan.f90


clean:
	rm -f *.o *.mod mstm-scan.exe
