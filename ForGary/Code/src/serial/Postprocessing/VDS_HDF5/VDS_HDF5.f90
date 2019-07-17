!>
!> \file VDS
!> @brief
!> @details

PROGRAM VDS
  use modRWHDF5
  use MFileIO
  implicit none

  integer ::   i, j, k, kk, k2, kp, n, nn, imax, jmax, kmax, ires,jres, icurrent,jcurrent, kcurrent
  integer, dimension(:), allocatable :: isize, jsize
  real(8), parameter :: Univ = 8314.3D0
  integer :: iConvert_grid, iConvert_flow
  integer :: iFormat_grid, iFormat_flow
  character(300) :: gridfile_rd, gridfile_wt, flowpath_rd, flowpath_wt
  integer :: file_be, file_end, file_skip, num_file,num_links
  character(400) :: fname, VDSname
  character(8) :: fnum
  character(4) :: fnum4_1, fnum4_2
  type(tp_rdwt_hdf5) :: grd, fsol
  type(tp_rdwt_hdf5), dimension(:), allocatable :: FILES
  integer :: kioplane, inode, jnode
  character(5) :: gridname(3) = (/'x','y','z'/)
  character(5) :: flowdataname(5) = (/'u','v','w','p','T'/)
  integer(HID_T) :: file_id, group_id, dset_id, dspace_id, memspace_id, prp_id, src_space_id
  integer :: hdferr

  call InitHDF5()
  call Initial()

  
  if(iConvert_flow.eq.1) then ! Convert flow file
    call InitFlowHDF5(fsol)
    num_file = (file_end - file_be)/file_skip + 1
    
    
    
    if(iFormat_flow.eq.1) then !kioplane selection
      
      !finding amount of files
      k=kmax/kioplane
      ires=kmax-k*kioplane
      num_links = k
      if(ires.gt.0) num_links=k+1
      
      allocate(FILES(num_links))
      
      

      do n=1, num_file
        write(unit=fnum,fmt='(I08.8)') (file_be + (n-1)*file_skip)
        
        kcurrent=0
        do kk=1, kmax, kioplane !defining info of each file to be linked
          kcurrent=kcurrent+1    
          call InitFlowHDF5(FILES(kcurrent))
           
          k2 = min(kk+kioplane-1,kmax)
          kp = k2-kk+1
          
          FILES(kcurrent)%dimsf = (/kp,imax,jmax/)
          FILES(kcurrent)%offset = (/kk-1,0,0/)
          
          !generating expected file name
          write(unit=fnum4_1,fmt='(I04.4)') kk
          write(unit=fnum4_2,fmt='(I04.4)') k2
          FILES(kcurrent)%fname = 'flowdata_'//fnum//'_kplane'//fnum4_1//'_'//fnum4_2//'.h5'
    
        enddo!end kk loop
        
        CALL MakeVDS()
        
        
      enddo ! end n loop
          
    deallocate(FILES)
    
    endif!i_formatt =1

    
    if(iFormat_flow.eq.2) then !inode
      

      num_links=inode
      allocate(FILES(num_links))
      
     
      allocate(isize(0:inode-1))
      
      !defining size of each file
      i=imax/inode
      ires = imax - i*inode
      do n=0, inode-1
        isize(n) = i
        if(n.lt.ires) isize(n) = i+1
      enddo

      do n=1, num_file !loop each file group
        write(unit=fnum,fmt='(I08.8)') (file_be + (n-1)*file_skip)
        

        icurrent = 0
        do nn=0,inode-1 !loop each file to generate information for VDS
          call InitFlowHDF5(FILES(nn+1)) 
          
          FILES(nn+1)%offset=(/0,icurrent,0/)
          
          icurrent = icurrent + isize(nn)
          FILES(nn+1)%dimsf = (/kmax,isize(nn),jmax/)
          write(unit=fnum4_1,fmt='(I04.4)') icurrent-isize(nn)+1
          write(unit=fnum4_2,fmt='(I04.4)') icurrent
          FILES(nn+1)%fname = 'flowdata_i'//fnum4_1//'_'//fnum4_2//'_'//fnum//'.h5'
          

        enddo!end nn loop defining each file data

        CALL MakeVDS()
      enddo ! end n loop
      
      deallocate(isize)        
      deallocate(FILES)
    
    endif!i_formatt =2
    
    if(iFormat_flow.eq.3) then !jnode

      num_links=jnode
      allocate(FILES(num_links))
      
     
      allocate(jsize(0:jnode-1))
      
      !defining size of each file
      j=jmax/jnode
      jres = jmax - j*jnode
      do n=0, jnode-1
        jsize(n) = j
        if(n.lt.jres) jsize(n) = j+1
      enddo

      do n=1, num_file
        write(unit=fnum,fmt='(I08.8)') (file_be + (n-1)*file_skip)
        print *, jsize
        
      
        jcurrent = 0
        do nn=0,jnode-1 !loop each file to generate information for VDS
          call InitFlowHDF5(FILES(nn+1))
          
          FILES(nn+1)%offset=(/0,0,jcurrent/)
          
          jcurrent = jcurrent + jsize(nn)
          FILES(nn+1)%dimsf = (/kmax,imax,jsize(nn)/)
          write(unit=fnum4_1,fmt='(I04.4)') jcurrent-jsize(nn)+1
          write(unit=fnum4_2,fmt='(I04.4)') jcurrent
          FILES(nn+1)%fname = 'flowdata_j'//fnum4_1//'_'//fnum4_2//'_'//fnum//'.h5'
          

        enddo!end nn loop defining each file data

        CALL MakeVDS()
      enddo ! end n loop
      
      deallocate(jsize)        
      deallocate(FILES)
    
    endif!i_formatt =3


  endif  ! end if iConvert_flow

 
  call FinalizeHDF5()

contains

  subroutine Initial()
     implicit none
     read(*,*)
     read(*,*) imax, jmax, kmax, kioplane, inode, jnode
     read(*,*)
     read(*,*) iConvert_grid, iFormat_grid
     read(*,*)
     read(*,*) iConvert_flow, iFormat_flow
     read(*,*)
     read(*,*)
     read(*,*)
     read(*,'(a)') gridfile_rd
     read(*,*)
     read(*,'(a)') gridfile_wt
     read(*,*)
     read(*,*)
     read(*,*)
     read(*,'(a)') flowpath_rd
     read(*,*)
     read(*,'(a)') flowpath_wt
     read(*,*)
     read(*,*)
     read(*,*)
     read(*,*) file_be, file_end, file_skip

   if(iConvert_grid.eq.1.or.iConvert_flow.eq.1) then
       print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
       if(iFormat_grid.eq.1.or.iFormat_flow.eq.1) then
         print *, 'Make VDS file from multiple HDF5 files using kioplane'
       elseif(iFormat_grid.eq.2.or.iFormat_flow.eq.2) then
         print *, 'Make VDS file from multiple HDF5 files using inode'
       elseif(iFormat_grid.eq.3.or.iFormat_flow.eq.3) then
         print *, 'Make VDS file from multiple HDF5 files using jnode'
       else
         print *, 'unknown grid format. Stop!'
         stop
       endif
       print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
     endif

     print *, 'File dimension: imax = ', imax, 'jmax = ', jmax, 'kmax = ', kmax

  end subroutine Initial
  
  subroutine MakeVDS()
  
  
  
  
   fsol%dimsf=(/kmax,imax,jmax/)
   VDSname= 'VDS_flowdata_'//fnum//'.h5'
      
        
   CALL h5fcreate_f(VDSname, H5F_ACC_TRUNC_F, file_id, hdferr)
   CALL h5screate_simple_f(fsol%rank, fsol%dimsf, dspace_id, hdferr)
   CALL H5Pcreate_f(H5P_DATASET_CREATE_F,prp_id,hdferr)
   CALL H5Pset_fill_value_f(prp_id, H5T_NATIVE_DOUBLE, -1.,hdferr)
        
      do i=1,fsol%dnum  
      
      print *, 'linking dataset...', fsol%dname(i)
      icurrent = 0
        do nn=1,num_links
          
          fsol%offset=FILES(nn)%offset
          
          
          fsol%dimsf = FILES(nn)%dimsf

          fsol%fname = FILES(nn)%fname
          print *, 'linking file ... ', trim(fsol%fname)

        



    CALL h5screate_simple_f(fsol%rank, fsol%dimsf, src_space_id, hdferr) ! create source space
    CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, fsol%offset, &    
                               fsol%count, hdferr, fsol%stride, fsol%dimsf) !select location to be linked
    CALL h5pset_virtual_f(prp_id,dspace_id,fsol%fname,fsol%dname(i), src_space_id, hdferr) !set link
    
    
    CALL h5sclose_f(src_space_id, hdferr)
       
enddo!end nn loop
    CALL h5dcreate_f(file_id, fsol%dname(i), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr,prp_id ) !create data set    
    CALL h5dclose_f(dset_id, hdferr)
    
           enddo !end i loop
           
    CALL h5sclose_f(dspace_id, hdferr)
    Call h5pclose_f(prp_id, hdferr) 
    
    
  print *, 'linking scalar data ',fsol%sname, ' using file...' ,fsol%fname
  
  CALL H5Pcreate_f(H5P_DATASET_CREATE_F,prp_id,hdferr)
  CALL H5Pset_fill_value_f(prp_id, H5T_NATIVE_DOUBLE, -1.,hdferr)
  
  !adding link to scalar         
  CALL h5screate_f(H5S_SCALAR_F, dspace_id, hdferr) !Making scalar data space
  CALL h5screate_f(H5S_SCALAR_F, src_space_id, hdferr) !Source file
  CALL h5pset_virtual_f(prp_id,dspace_id,fsol%fname,fsol%sname, src_space_id, hdferr) !making link
  CALL h5dcreate_f(file_id, trim(fsol%sname),H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr,prp_id)
  !closing scalar 
  CALL h5dclose_f(dset_id, hdferr)
  CALL h5sclose_f(dspace_id, hdferr)
  CALL h5sclose_f(src_space_id, hdferr)
  Call h5pclose_f(prp_id, hdferr) 
  
  
  !closing file  
  CALL h5fclose_f(file_id, hdferr)

      
      
  end subroutine MakeVDS
  
  

end program VDS




