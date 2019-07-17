fftutil.f90
contains subroutines to rearrange the fourier series output by subroutines
in modFFT.f90 so that the zero-frequency component is at the center of spectrum.

interp1d.f90
contains one-dimensional interpolation subroutines and related subroutines,
including bi-section search, linear and spline interpolation, derivatives
and integral computed through spline interpolation.

modBuffer.f90
contains module to generate a buffer zone using tanh shape.

modFFT.f90
contains modules that wrap the fftpack5 subroutines upto 3D FFT.

modFileIO_Compatible.f90
contains subroutines to read/write 2D/3D Plot3D unformatted files and unformatted
files for the cartesian grid.

modTecplotIO.f90
contains subroutines to read/write Tecplot ASCII format files.

nonlinearsolver.f90
contains Newton-Rapson nonlinear solver.

random.f90
contains subroutines to generate uniform and Gaussian random numbers.

spectprocess.f90
contains subroutines to do spectral interpolation.

tridiagonal.f90
contains subroutines to solve tri-diagonal systems or inverse tri-diagonal matrices.
