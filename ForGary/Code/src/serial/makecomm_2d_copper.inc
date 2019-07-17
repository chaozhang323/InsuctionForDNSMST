FC=ftn
#FC=mpif90
#FC=ifort
#FFLAG=-O -cpp -r8 -zero -traceback -heap-arrays -g
FFLAG=-O3 -cpp -r8 -zero -traceback -heap-arrays -g
#FFLAG=-r8 -check all,noarg_temp_created -align -fpe0 -g -O0 -traceback -heap-arrays -debug all -ftrapuv
#FFLAG=-O3 -cpp -axAVX -xSSE4.1 -r8 -zero -align -traceback -heap-arrays -fpe0 -ftrapuv -mcmodel medium
#FFLAG=-O3 -axAVX -xSSE4.1 -ip -r8 -zero -traceback -heap-arrays
#FFLAG=-O3 -Mpreprocess -r8

HDF5=/opt/cray/hdf5-parallel/1.8.13/INTEL/140
tec360=${DIR}/../tecio
fftpack=${DIR}/../fftpack5.1
2decomp_fft=/u/czb58/Package/2decomp_fft

INCLUDE=-I $(HDF5)/include  -I $(tec360)/include -I $(2decomp_fft)/include 
LIB=-L$(HDF5)/lib -L $(tec360)/lib  -L$(2decomp_fft)/lib      -L$(fftpack)/lib  -lfftpack
LIBLINK= $(HDF5)/lib/libhdf5_fortran.a $(tec360)/lib/libtecio.a  $(2decomp_fft)/lib/lib2decomp_fft.a  -lhdf5   -lz  -lm -lstdc++


#FC=gfortran
#FFLAG=-O3 -fdefault-real-8
LINKFLAG=
#LINKFLAG=/nasa/ae/tecplot360/lib/tecio64.a
