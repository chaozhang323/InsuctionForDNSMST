!subroutine to rearrange spectra output by 1D FFT so that
!the zero-frequency component is at the center of the series.
!
!np (input): number of points in the spectrum
!spectin (input): the spectrum output by 1D FFT
!spectout (output): the rearranged spectrum
!wavenum (output): the corresponding wave number of the output spectrum
subroutine RearrangeFFT1D(np,spectin,spectout,wavenum)
  integer, intent(in) :: np
  complex, dimension(np), intent(inout) :: spectin
  complex, dimension(np), intent(out) :: spectout
  integer, dimension(np), intent(out) :: wavenum
  integer :: ihalf, i
  ihalf = np/2+1
  do i=1,np
    wavenum(i) = ihalf-1+i-np
    if (i.le.np-ihalf) then
      spectout(i) = spectin(ihalf+i)
    else
      spectout(i) = spectin(i-np+ihalf)
    end if
  end do
end subroutine RearrangeFFT1D

!subroutine to rearrange spectra output by 1D FFT so that
!the zero-frequency component is at the center of the series.
!
!np1 (input): number of points in 1st dimension of the spectrum
!np2 (input): number of points in 2nd dimension of the spectrum
!spectin (input): the spectrum output by 2D FFT
!spectout (output): the rearranged spectrum
!wavenum1 (output): the corresponding wave number in the 1st dimension of the output spectrum
!wavenum2 (output): the corresponding wave number in the 2nd dimension of the output spectrum
subroutine RearrangeFFT2D(np1,np2,spectin,spectout,wavenum1,wavenum2)
  integer, intent(in) :: np1, np2
  complex, dimension(np1,np2), intent(inout) :: spectin
  complex, dimension(np1,np2), intent(out) :: spectout
  complex, dimension(np2) :: tmp
  integer, intent(out) :: wavenum1(np1),wavenum2(np2)
  integer :: ihalf, jhalf, i, j
  ihalf = np1/2+1; jhalf = np2/2+1
  do i=1,np1
    wavenum1(i) = ihalf-1+i-np1
    if (i.le.np1-ihalf) then
      spectout(i,:) = spectin(ihalf+i,:)
    else
      spectout(i,:) = spectin(i-np1+ihalf,:)
    end if
  end do
  do j=1,np2
    wavenum2(j) = jhalf-1+j-np2
  end do
  do i=1,np1
    tmp = spectout(i,:)
    do j=1,np2
    if (j.le.np2-jhalf) then
      spectout(i,j) = tmp(jhalf+j)
    else
      spectout(i,j) = tmp(j-np2+jhalf)
    end if
    end do
  end do
end subroutine RearrangeFFT2D
