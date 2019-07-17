FC=mpif90
#FC=ifort
#FFLAG=-O -cpp -r8 -zero -traceback -heap-arrays -g
#FFLAG=-O3 -cpp -r8 -zero -traceback -heap-arrays -g
#FFLAG=-r8 -check all,noarg_temp_created -align -fpe0 -g -O0 -traceback -heap-arrays -debug all -ftrapuv
FFLAG=-O3 -cpp -axAVX -xSSE4.1 -r8 -zero -align -traceback -heap-arrays -fpe0 -ftrapuv -mcmodel medium
#FFLAG=-O3 -axAVX -xSSE4.1 -ip -r8 -zero -traceback -heap-arrays

# For Lib at nic-cluster ###
HDF5=/share/duan/czb58/Package/HDF5_Szip
2decomp_fft=/share/duan/czb58/Package/2decomp_fft_new
T3PIO=/share/duan/czb58/Package/T3PIO

# For Lib at forge ####
#HDF5=/opt/hdf5/intel2015.2.164/openmpi_ib
#2decomp_fft=/home/czb58/Package/2decomp_fft
#T3PIO=/home/czb58/Package/t3pio

# For Lib at AirForce US ###
#HDF5=/u/US_HOME2/czb58/Package/PHDF5
#2decomp_fft=/u/home/czb58/Package/2decomp_fft

tec360=$(DIR)/../tecio
fftpack=$(DIR)/../fftpack5.1

INCLUDE=-I $(HDF5)/include -I $(COMMDIR) -I $(tec360)/include -I $(2decomp_fft)/include 
LIB=-L$(HDF5)/lib -L $(tec360)/lib  -L$(2decomp_fft)/lib      -L$(fftpack)/lib  -lfftpack
LIBLINK= $(HDF5)/lib/libhdf5_fortran.a $(tec360)/lib/libtecio.a  $(2decomp_fft)/lib/lib2decomp_fft.a  -lhdf5   -lz  -lm -lstdc++


#FC=gfortran
#FFLAG=-O3 -fdefault-real-8
LINKFLAG=
#LINKFLAG=/nasa/ae/tecplot360/lib/tecio64.a
