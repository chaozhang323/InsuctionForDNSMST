#FC=mpif90
FC=ifort
#FFLAG=-O3 -cpp -r8 -g
#FFLAG=-r8 -check all,noarg_temp_created -align -fpe0 -g -O0 -traceback -heap-arrays -debug all -ftrapuv
#FFLAG=-O3 -assume buffered_io -axAVX -xSSE4.1 -r8 -zero -align -traceback -heap-arrays -fpe0 -ftrapuv -mcmodel medium
#FFLAG=-O3 -axAVX -xSSE4.1 -ip -r8 -zero -traceback -heap-arrays
FFLAG=-O3  -r8 -zero -traceback -heap-arrays -mcmodel medium

HDF5=/u/US_HOME2/czb58/Package/HDF5

#tec360=/home/czb58/tecplot360
#fftpack=/share/duan/czb58/dnsmst_utilities/src/serial/fftpack5.1

tec360=$(DIR)/../tecio
fftpack=$(DIR)/../fftpack5.1


INCLUDE=-I $(HDF5)/include -I $(COMMDIR) -I $(tec360)/include
LIB=-L$(HDF5)/lib -L $(tec360)/lib  -L$(fftpack)/lib  -lfftpack
LIBLINK= $(HDF5)/lib/libhdf5_fortran.a $(tec360)/lib/libtecio.a  -lhdf5   -lz  -lm -lstdc++



#FC=gfortran
#FFLAG=-O3 -fdefault-real-8
#LINKFLAG=
#LINKFLAG=/nasa/ae/tecplot360/lib/tecio64.a
