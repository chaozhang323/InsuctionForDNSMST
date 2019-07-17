# Set to appropriate C++ compiler
           CPP=g++
        MPICPP=mpic++

#
# If linking against the mpi tecio built using the source code care package then
# the library is static instead of dynamic.  Use the following instead:
# TECIOMPILIB=full-path-to-care-package/libteciompi.a
#
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	TECIOMPILIB=../../../../MacOS/libteciompi.dylib
else
	TECIOMPILIB=../../../../bin/libteciompi.so
endif

#
# If linking against the tecio built using the source code care package then
# the library is static instead of dynamic.  Use the following instead:
# TECIOLIB=full-path-to-care-package/libtecio.a
#
ifeq ($(UNAME_S),Darwin)
	TECIOLIB=../../../../MacOS/libtecio.dylib
else
	TECIOLIB=../../../../bin/libtecio.so
endif

# Set to appropriate Fortran compiler
            FC=gfortran
         MPIFC=mpif90
        FFLAGS=-fcray-pointer

ifeq ($(UNAME_S),Darwin)
    EXTRALIBS=-lstdc++
else
    FOUND_INSTALLED_LIBSTDCXX_S := $(shell test -f ../../../../bin/sys/libstdc++.so.6 && echo found || echo missing)
    #
    # Note:
    #     On Linux the examples must link against the libstdc++.so.6 that libtecio.so was build
    #     with. For customer installations this is located up and over in the bin directory of the
    #     installation. For internal developer builds, libstdc++.so.6 isn't in the CMake binary
    #     directory and is located with the compiler.
    #
    ifeq ($(FOUND_INSTALLED_LIBSTDCXX_S),found)
        EXTRALIBS=../../../../bin/sys/libstdc++.so.6
    else
        EXTRALIBS=-lstdc++
    endif
endif

   EXTRAINCLUDES=-I../../../../include
#
# (If not needed reset to empty in actual Makefile)
# link libraries
      LINKLIBS=-lpthread

   PATHTOEXECUTABLE=.

all: $(TARGETS)

clean:
	rm -f $(PATHTOEXECUTABLE)/$(EXECUTABLE)

cppbuild:
	$(CPP) $(CPPFILES) $(EXTRAINCLUDES) $(TECIOLIB) $(EXTRALIBS) $(LINKLIBS) -o $(PATHTOEXECUTABLE)/$(EXECUTABLE)

mpicppbuild:
	$(MPICPP) -DTECIOMPI $(CPPFILES) $(EXTRAINCLUDES) $(TECIOMPILIB) $(EXTRALIBS) $(LINKLIBS) -o $(PATHTOEXECUTABLE)/$(EXECUTABLE)-mpi

fbuild:
	$(FC) $(FFILES) $(FFLAGS) $(EXTRAINCLUDES) $(TECIOLIB) $(EXTRALIBS) $(LINKLIBS) -o $(PATHTOEXECUTABLE)/$(EXECUTABLE)-f

f90build:
	$(FC) $(F90FILES) $(FFLAGS) $(EXTRAINCLUDES) $(TECIOLIB) $(EXTRALIBS) $(LINKLIBS) -o $(PATHTOEXECUTABLE)/$(EXECUTABLE)-f90

mpif90build:
	$(MPIFC) -DTECIOMPI $(F90FILES) $(FFLAGS) $(EXTRAINCLUDES) $(TECIOMPILIB) $(EXTRALIBS) $(LINKLIBS) -lmpi_cxx -o $(PATHTOEXECUTABLE)/$(EXECUTABLE)-mpif90
