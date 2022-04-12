program vector2tile_perturbation

  use netcdf
  implicit none

  integer, parameter :: max_n_var_lndp = 6

  type vector_type
    double precision, allocatable :: var                (:)
  end type vector_type

  type tile_type
    double precision, allocatable :: var                (:,:,:)
    real,             allocatable :: land_frac          (:,:,:)
  end type tile_type

  type namelist_type
    character*256      :: namelist_name = ""
    character*15       :: direction = ""
    character*3        :: layout
    character*256      :: tile_path = ""
    integer            :: tile_size
    character*256      :: input_file = ""
    character*256      :: output_file = ""
    character(len=3)   :: lndp_var_list(max_n_var_lndp)
    integer            :: n_var_lndp
  end type namelist_type
  type(vector_type)   :: vector
  type(tile_type)     :: tile
  type(namelist_type) :: namelist
  character*256       :: vector_filename
  character*256       :: tile_filename
  character*256       :: input_filename
  character*256       :: output_filename
  character*2         :: tile1, tile2
  real, allocatable   :: tmp2d(:,:)
  integer             :: filename_length
  integer             :: vector_length = 0
  integer             :: layout_x, layout_y, nx, ny
  integer             :: itile, ix, iy, iloc, ivar
  integer             :: i, j, m, n, i1, i2, j1, j2, t2
  integer             :: ncid, dimid, varid, status
  integer             :: dim_id_xdim, dim_id_ydim, dim_id_time
  integer             :: ncid_landp, ncid_vec, ncid_tile(6)
  logical             :: file_exists

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get namelist file name from command line
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call get_command_argument(1, namelist%namelist_name)
  if(namelist%namelist_name == "") then
        print *,  "add namelist to the command line: "
        stop 10
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read namelist information and create date string
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call ReadNamelist(namelist)

  print*, "conversion direction: ",namelist%direction

  if(namelist%direction /= "layout2vector" .and. namelist%direction /= "layout2tile" ) then
    print*, "conversion direction: ",namelist%direction, " not recognized"
    stop 10
  end if

  if(trim(namelist%layout) == '1x4') then
    layout_x = 1
    layout_y = 4
  else if(trim(namelist%layout) == '2x2') then
    layout_x = 2
    layout_y = 2
  else if(trim(namelist%layout) == '4x1') then
    layout_x = 4
    layout_y = 1
  else
    print*, "layout: ",namelist%layout, " not recognized"
    stop 10
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocate tile variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(tile%var(      namelist%tile_size,namelist%tile_size,6))
  allocate(tile%land_frac(namelist%tile_size,namelist%tile_size,6))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read FV3 tile information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do itile = 1, 6

    if(namelist%tile_size < 100) then
      write(tile_filename,'(a5,i2,a11,i1,a3)') "oro_C", namelist%tile_size, ".mx100.tile", itile, ".nc"
    elseif(namelist%tile_size < 1000) then
      write(tile_filename,'(a5,i3,a11,i1,a3)') "oro_C", namelist%tile_size, ".mx100.tile", itile, ".nc"
    elseif(namelist%tile_size < 10000) then
      write(tile_filename,'(a5,i4,a11,i1,a3)') "oro_C", namelist%tile_size, ".mx100.tile", itile, ".nc"
    else
      print *, "unknown tile size"
      stop 10
    end if

    tile_filename = trim(namelist%tile_path)//trim(tile_filename)

    inquire(file=tile_filename, exist=file_exists)

    if(.not.file_exists) then
      print*, trim(tile_filename), " does not exist"
      print*, "Check paths and file name"
      stop 10
    end if

    status = nf90_open(tile_filename, NF90_NOWRITE, ncid)
      if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_varid(ncid, "land_frac", varid)
    status = nf90_get_var(ncid, varid , tile%land_frac(:,:,itile))

    status = nf90_close(ncid)

    vector_length = vector_length + count(tile%land_frac(:,:,itile) > 0)

  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocate vector variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(vector%var(vector_length))

  nx = namelist%tile_size/layout_x
  ny = namelist%tile_size/layout_y

  if(namelist%n_var_lndp > 0) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define the output file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if(trim(namelist%direction) == "layout2vector") then
        output_filename = namelist%output_file
        status = nf90_create(output_filename, NF90_CLOBBER, ncid_vec)
           if (status /= nf90_noerr) call handle_err(status)

! Define dimensions in the file.

        status = nf90_def_dim(ncid_vec, "location", vector_length, dim_id_xdim)
           if (status /= nf90_noerr) call handle_err(status)
        status = nf90_def_dim(ncid_vec, "Time", NF90_UNLIMITED, dim_id_time)
           if (status /= nf90_noerr) call handle_err(status)

! Define variables in the file. 
        do ivar = 1, namelist%n_var_lndp
           status = nf90_def_var(ncid_vec, namelist%lndp_var_list(ivar), &
                    NF90_FLOAT, (/dim_id_xdim,dim_id_time/), varid)
              if (status /= nf90_noerr) call handle_err(status)
        enddo
        status = nf90_enddef(ncid_vec)
     else if(trim(namelist%direction) == "layout2tile" ) then
        filename_length = len_trim(namelist%output_file)
        do itile = 1, 6
           write(tile1,fmt='(I2.2)') itile
           output_filename = trim(namelist%output_file(1:filename_length-5))//trim(tile1)//'.nc'
           status = nf90_create(output_filename, NF90_CLOBBER, ncid_tile(itile))
              if (status /= nf90_noerr) call handle_err(status)

! Define dimensions in the file.

           status = nf90_def_dim(ncid_tile(itile), "xaxis_1", namelist%tile_size , dim_id_xdim)
              if (status /= nf90_noerr) call handle_err(status)
           status = nf90_def_dim(ncid_tile(itile), "yaxis_1", namelist%tile_size , dim_id_ydim)
              if (status /= nf90_noerr) call handle_err(status)
           status = nf90_def_dim(ncid_tile(itile), "Time", NF90_UNLIMITED, dim_id_time)
              if (status /= nf90_noerr) call handle_err(status)

! Define variables in the file. 
           do ivar = 1, namelist%n_var_lndp
              status = nf90_def_var(ncid_tile(itile), namelist%lndp_var_list(ivar), &
                       NF90_FLOAT, (/dim_id_xdim,dim_id_ydim,dim_id_time/), varid)
                 if (status /= nf90_noerr) call handle_err(status)
           enddo
           status = nf90_enddef(ncid_tile(itile))
        enddo
     endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocate temporary variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     allocate(tmp2d(nx,ny))

     do ivar = 1, namelist%n_var_lndp
        t2 = 1
        filename_length = len_trim(namelist%input_file)
        do itile = 1, 6
           i1=1
           i2=i1+nx-1
           j1=1
           j2=j1+ny-1
           do j=1,layout_y
              do i=1,layout_x
                 write(tile2,fmt='(I2.2)') t2
                 if(t2 > 1) then
                    i1=i1+nx
                    i2=i2+nx
                    if (i2 .GT. namelist%tile_size) then
                       i1=1
                       i2=i1+nx-1
                    endif
                 endif
                 input_filename = trim(namelist%input_file(1:filename_length-5))//trim(tile2)//'.nc'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read the perturbation pattern
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                 status = nf90_open(input_filename, NF90_NOWRITE, ncid_landp)
                 status = nf90_inq_varid(ncid_landp, namelist%lndp_var_list(ivar), varid)
                 status = nf90_get_var(ncid_landp, varid, tmp2d, start = (/1,1,1/), count = (/nx, ny, 1/))

                 do m = i1, i2
                    do n = j1, j2
                       tile%var(m,n,itile) = tmp2d(m-i1+1,n-j1+1)
                    enddo
                 enddo
                 t2 = t2+1
              enddo
              j1=j1+ny
              j2=j2+ny
              if (j2 .GT. namelist%tile_size) then
                 j1=1
                 j2=j1+ny-1
              endif
           enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write the perturbation pattern for the tile files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           if(trim(namelist%direction) == "layout2tile") then
              status = nf90_inq_varid(ncid_tile(itile), namelist%lndp_var_list(ivar), varid)
              status = nf90_put_var(ncid_tile(itile), varid , tile%var(:,:,itile), &
                                    start = (/1,1,1/), count = (/namelist%tile_size, namelist%tile_size, 1/))
           endif
        enddo    ! for each tile
        iloc = 0
        do itile = 1, 6
           do j = 1, namelist%tile_size
              do i = 1, namelist%tile_size
                 if(tile%land_frac(i,j,itile) > 0.0) then
                    iloc = iloc + 1
                    vector%var(iloc) = tile%var(i,i,itile)
                 endif
              enddo
           enddo
        enddo    ! for each tile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write the perturbation pattern for the vector file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if(trim(namelist%direction) == "layout2vector") then
           status = nf90_inq_varid(ncid_vec, namelist%lndp_var_list(ivar), varid)
           status = nf90_put_var(ncid_vec, varid , vector%var(:), &
                              start = (/1,1/), count = (/vector_length, 1/))
        endif
     enddo       ! for each variable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Close the netcdf file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if(trim(namelist%direction) == "layout2vector") then
        status = nf90_close(ncid_vec)
     else if(trim(namelist%direction) == "layout2tile") then
        do itile = 1, 6
           status = nf90_close(ncid_tile(itile))
        enddo
     endif
  endif

contains

  subroutine ReadNamelist(namelist)

    integer, parameter  :: max_n_var_lndp = 6
    type(namelist_type) :: namelist
    character*15        :: direction
    character*3         :: layout
    character*256       :: tile_path
    integer             :: tile_size
    character*256       :: input_file
    character*256       :: output_file
    character(len=3)    :: lndp_var_list(max_n_var_lndp)
    integer             :: n_var_lndp
    integer             :: k

    namelist / run_setup  / direction, layout, tile_path, tile_size, input_file, output_file, lndp_var_list, n_var_lndp

    lndp_var_list = 'XXX'

    open(30, file=namelist%namelist_name, form="formatted")
     read(30, run_setup)
    close(30)

    namelist%direction   = direction
    namelist%layout      = layout
    namelist%tile_path   = tile_path
    namelist%tile_size   = tile_size
    namelist%input_file  = input_file
    namelist%output_file = output_file
    n_var_lndp= 0
    do k =1,size(lndp_var_list)
       if (lndp_var_list(k) .EQ. 'XXX') then
          cycle
       else
          n_var_lndp=n_var_lndp+1
          namelist%lndp_var_list(n_var_lndp) = lndp_var_list(k)
       endif
    enddo
    namelist%n_var_lndp = n_var_lndp
    if (n_var_lndp > max_n_var_lndp) then
       print*, 'ERROR: land perturbation requested for too many parameters', &
                       'increase max_n_var_lndp'
       stop 10
    endif

  end subroutine ReadNamelist

  subroutine handle_err(status)
    use netcdf
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 10
    end if
  end subroutine handle_err

end program vector2tile_perturbation

!integer, parameter   :: location=18360
!integer, parameter   :: CRES=96
!integer, parameter   :: ntiles=6
!double precision     :: glon(cres,cres,ntiles), glat(cres,cres,ntiles)
!real                 :: rlon(cres,cres,ntiles), rlat(cres,cres,ntiles)
!real                 :: vgf(cres,cres,ntiles), smc(cres,cres,ntiles)
!real                 :: vvgf(location), vsmc(location)
!real, allocatable    :: xlon(:,:), xlat(:,:), xvgf(:,:), xsmc(:,:)
!real                 :: err_lon, err_lat
!integer              :: cube_tile(location), cube_i(location), cube_j(location)
!character*128        :: dir, ifname, ofname, sfname, vfname
!character*10         :: compress
!character*3          :: layout
!character*2          :: tile1,tile2
!integer              :: layout_x, layout_y
!integer              :: nx, ny, incid, oncid, sncid, vncid
!integer              :: i, j, i1, i2, j1, j2, t, t2, m, n
!layout = '1x4'
!compress = 'false'
!dir = '/scratch1/NCEPDEV/stmp4/Zhichang.Guo/Work/Test/vp/phil/stochastic_physics/unit_tests/stochy_out/'
!vfname = '/scratch1/NCEPDEV/stmp2/Michael.Barlage/forcing/C96/static/ufs-land_C96_static_fields.nc'
!sfname = '/scratch1/NCEPDEV/stmp4/Zhichang.Guo/Work/tmp/fv3-jedi-tools/src/raster/fv3grid/fv3grid_c0096.nc4'
!ofname = './perturbation_pattern.nc'
!call nf90i_open_file('R', vfname, vncid)
!call nf90i_read1d_int(vncid, location, cube_tile(:), 'cube_tile')
!call nf90i_read1d_int(vncid, location, cube_i(:),    'cube_i')
!call nf90i_read1d_int(vncid, location, cube_j(:),    'cube_j')
!call nf90i_close_file(vncid)
!call nf90i_open_file('R', sfname, sncid)
!if(trim(layout) == '1x4') then
!   layout_x = 1
!   layout_y = 4
!else if(trim(layout) == '2x2') then
!   layout_x = 2
!   layout_y = 2
!else if(trim(layout) == '4x1') then
!   layout_x = 4
!   layout_y = 1
!endif
!allocate(xlon(nx,ny))
!allocate(xlat(nx,ny))
!allocate(xvgf(nx,ny))
!allocate(xsmc(nx,ny))
!t2=1
!do t = 1, ntiles
!   call nf90i_read2d_double_time(sncid, cres, cres, t, glon(:,:,t), 'flons')
!   call nf90i_read2d_double_time(sncid, cres, cres, t, glat(:,:,t), 'flats')
!   write(tile1,fmt='(I2.2)') t
!!   i1=1
!   i2=i1+nx-1
!   j1=1
!!   j2=j1+ny-1
!   do j=1,layout_y
!      do i=1,layout_x
!         write(tile2,fmt='(I2.2)') t2
!         if(t2 > 1) then
!            i1=i1+nx
!            i2=i2+nx
!            if (i2.GT.CRES) then
!               i1=1
!               i2=i1+nx-1
!            endif
!         endif
!!         ifname = trim(dir)//'workg_T162_984x488.tile'//trim(tile2)//'.nc'
!         call nf90i_open_file('R', ifname, incid)
!         call nf90i_read2d_real_time(incid, nx, ny, 1, xlon(:,:), 'grid_lon')
!         call nf90i_read2d_real_time(incid, nx, ny, 1, xlat(:,:), 'grid_lat')
!         call nf90i_read2d_real_time(incid, nx, ny, 1, xvgf(:,:), 'vgf')
!         call nf90i_read2d_real_time(incid, nx, ny, 1, xsmc(:,:), 'smc')
!         call nf90i_close_file(incid)
!         do m = i1, i2
!            do n = j1, j2
!               rlon(m,n,t) = xlon(m-i1+1,n-j1+1)
!               rlat(m,n,t) = xlat(m-i1+1,n-j1+1)
!               vgf(m,n,t)  = xvgf(m-i1+1,n-j1+1)
!               smc(m,n,t)  = xsmc(m-i1+1,n-j1+1)
!            enddo
!         enddo
!         write(*,'(8i5)') t, t2, i1, i2, j1, j2
!!        write(*,'(A128)') ifname
!         t2=t2+1
!      enddo
!      j1=j1+ny
!      j2=j2+ny
!      if (j2.GT.CRES) then
!         j1=1
!!         j2=j1+ny-1
!      endif
!   enddo
!enddo
!call nf90i_close_file(sncid)
!do t = 1, ntiles
!   do i = 1, cres
!      do j = 1, cres
!         err_lon = (glon(i,j,t)-rlon(i,j,t))/glon(i,j,t)*1.E6
!         err_lat = (glat(i,j,t)-rlat(i,j,t))/glat(i,j,t)*1.E6
!         if(abs(err_lon) > 1. .or. err_lat > 1.) then
!            write(*,'(A9,3i4,6f10.5)') 'warning: ', t, i, j, glon(i,j,t), rlon(i,j,t), glat(i,j,t), rlat(i,j,t), err_lon, err_lat
!         else
!            write(*,'(3i4,6f10.5)') t, i, j, glon(i,j,t), rlon(i,j,t), glat(i,j,t), rlat(i,j,t), err_lon, err_lat
!         endif
!      enddo
!   enddo
!enddo
!do i = 1, location
!   vvgf(i) = vgf(cube_i(i),cube_j(i),cube_tile(i))
!   vsmc(i) = smc(cube_i(i),cube_j(i),cube_tile(i))
!enddo
!call nf90i_open_file('C', ofname, oncid)
!call nf90i_def_dim(oncid, location, 'location')
!call nf90i_def_var_type(oncid, 'int',   'location', 'cube_i',    compress=compress)
!call nf90i_def_var_type(oncid, 'int',   'location', 'cube_j',    compress=compress)
!call nf90i_def_var_type(oncid, 'int',   'location', 'cube_tile', compress=compress)
!call nf90i_def_var_type(oncid, 'float', 'location', 'vgf',       compress=compress)
!call nf90i_def_var_type(oncid, 'float', 'location', 'smc',       compress=compress)
!call nf90i_def_end(oncid)
!call nf90i_out1d_int( oncid, location, cube_i(:),    'cube_i')
!call nf90i_out1d_int( oncid, location, cube_j(:),    'cube_j')
!call nf90i_out1d_int( oncid, location, cube_tile(:), 'cube_tile')
!call nf90i_out1d_real(oncid, location, vvgf(:),      'vgf')
!call nf90i_out1d_real(oncid, location, vsmc(:),      'smc')
!call nf90i_close_file(oncid)
!deallocate(xlon)
!deallocate(xlat)
!deallocate(xvgf)
!deallocate(xsmc)
!print*, 'Program ended normally!'
!stop
