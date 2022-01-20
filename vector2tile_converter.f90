program vector2tile_converter

  use netcdf
  implicit none

  type vector_type
    double precision, allocatable :: swe                (:)
    double precision, allocatable :: snow_depth         (:)
    double precision, allocatable :: active_snow_layers (:)
    double precision, allocatable :: swe_previous       (:)
    double precision, allocatable :: snow_soil_interface(:,:)
    double precision, allocatable :: temperature_snow   (:,:)
    double precision, allocatable :: snow_ice_layer     (:,:)
    double precision, allocatable :: snow_liq_layer     (:,:)
    double precision, allocatable :: temperature_soil   (:,:)
! needed for IMSaggregate_mod
    double precision, allocatable :: vegetation_type(:)
! needed by JEDI to mask out land-ice
    double precision, allocatable :: soil_moisture     (:)
  end type vector_type    

  type tile_type
    double precision, allocatable :: swe                (:,:,:)
    double precision, allocatable :: snow_depth         (:,:,:)
    double precision, allocatable :: active_snow_layers (:,:,:)
    double precision, allocatable :: swe_previous       (:,:,:)
    double precision, allocatable :: snow_soil_interface(:,:,:,:)
    double precision, allocatable :: temperature_snow   (:,:,:,:)
    double precision, allocatable :: snow_ice_layer     (:,:,:,:)
    double precision, allocatable :: snow_liq_layer     (:,:,:,:)
    double precision, allocatable :: temperature_soil   (:,:,:,:)
    real,             allocatable :: land_frac          (:,:,:)
    double precision, allocatable :: soil_moisture      (:, :, :)
    double precision, allocatable :: vegetation_type(:,:,:)
! needed by add increments
    double precision, allocatable :: slmsk              (:, :, :)
  end type tile_type    

  type namelist_type
    character*256      :: namelist_name = ""
    character*11       :: direction = ""
    character*256      :: tile_path = ""
    integer            :: tile_size
    character*19       :: restart_date = ""
    character*256      :: vector_restart_path = ""
    character*256      :: tile_restart_path = ""
    character*256      :: output_path = ""
    character*256      :: static_filename = ""
  end type namelist_type    

  type(vector_type)   :: vector
  type(tile_type)     :: tile
  type(namelist_type) :: namelist
  character*256       :: vector_filename
  character*256       :: tile_filename
  character*19        :: date
  integer             :: vector_length = 0
  integer             :: yyyy,mm,dd,hh,nn,ss
  integer             :: itile, ix, iy, iloc
  integer             :: ncid, dimid, varid, status
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
  
  if(namelist%direction /= "tile2vector" .and. namelist%direction /= "vector2tile" ) then
    print*, "conversion direction: ",namelist%direction, " not recognized"
    stop 10 
  end if

  read(namelist%restart_date( 1: 4),'(i4.4)') yyyy
  read(namelist%restart_date( 6: 7),'(i2.2)') mm
  read(namelist%restart_date( 9:10),'(i2.2)') dd
  read(namelist%restart_date(12:13),'(i2.2)') hh
  read(namelist%restart_date(15:16),'(i2.2)') nn
  read(namelist%restart_date(18:19),'(i2.2)') ss

  write(date,'(i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)') &
      yyyy, "-", mm, "-", dd, "_", hh, "-", nn, "-", ss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocate tile variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(tile%swe                (namelist%tile_size,namelist%tile_size,6))
  allocate(tile%snow_depth         (namelist%tile_size,namelist%tile_size,6))
  allocate(tile%active_snow_layers (namelist%tile_size,namelist%tile_size,6))
  allocate(tile%swe_previous       (namelist%tile_size,namelist%tile_size,6))
  allocate(tile%snow_soil_interface(namelist%tile_size,namelist%tile_size,7,6))
  allocate(tile%temperature_snow   (namelist%tile_size,namelist%tile_size,3,6))
  allocate(tile%snow_ice_layer     (namelist%tile_size,namelist%tile_size,3,6))
  allocate(tile%snow_liq_layer     (namelist%tile_size,namelist%tile_size,3,6))
  allocate(tile%temperature_soil   (namelist%tile_size,namelist%tile_size,4,6))
  allocate(tile%soil_moisture      (namelist%tile_size,namelist%tile_size,6))
  allocate(tile%land_frac          (namelist%tile_size,namelist%tile_size,6))
  allocate(tile%slmsk              (namelist%tile_size,namelist%tile_size,6))
  allocate(tile%vegetation_type(namelist%tile_size,namelist%tile_size,6))

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
  
  print*, "The FV3 tiles report ",vector_length, "land grids"
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocate vector variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(vector%swe                (vector_length))
  allocate(vector%snow_depth         (vector_length))
  allocate(vector%active_snow_layers (vector_length))
  allocate(vector%swe_previous       (vector_length))
  allocate(vector%snow_soil_interface(vector_length,7))
  allocate(vector%temperature_snow   (vector_length,3))
  allocate(vector%snow_ice_layer     (vector_length,3))
  allocate(vector%snow_liq_layer     (vector_length,3))
  allocate(vector%temperature_soil   (vector_length,4))
  allocate(vector%soil_moisture      (vector_length))
  allocate(vector%vegetation_type    (vector_length))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Direction of transfer branch
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(namelist%direction == "vector2tile") then
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read vector restart file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call ReadVectorRestart(namelist, date, vector, vector_length)
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Transfer vector to tiles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! explicitly initialize to 0.
    tile%slmsk=0. 

    iloc = 0
    do itile = 1, 6
    do iy = 1, namelist%tile_size
    do ix = 1, namelist%tile_size
      if(tile%land_frac(ix,iy,itile) > 0.0) then
        iloc = iloc + 1
        tile%swe(ix,iy,itile)                   = vector%swe(iloc)
        tile%vegetation_type(ix,iy,itile)       = vector%vegetation_type(iloc)
        tile%snow_depth(ix,iy,itile)            = vector%snow_depth(iloc)
        tile%active_snow_layers(ix,iy,itile)    = vector%active_snow_layers(iloc)
        tile%swe_previous(ix,iy,itile)          = vector%swe_previous(iloc)
        tile%snow_soil_interface(ix,iy,:,itile) = vector%snow_soil_interface(iloc,:)
        tile%temperature_snow(ix,iy,:,itile)    = vector%temperature_snow(iloc,:)
        tile%snow_ice_layer(ix,iy,:,itile)      = vector%snow_ice_layer(iloc,:)
        tile%snow_liq_layer(ix,iy,:,itile)      = vector%snow_liq_layer(iloc,:)
        tile%temperature_soil(ix,iy,:,itile)    = vector%temperature_soil(iloc,:)
        tile%soil_moisture(ix,iy,itile)         = vector%soil_moisture(iloc)
        tile%slmsk(ix,iy,itile)         = 1.
      end if
      
    end do
    end do
    end do
      
   print*, "Transferred ",iloc, "land grids"  
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write FV3 tile file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call WriteTileRestart(namelist, date, tile)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! tile2vector branch
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  elseif(namelist%direction == "tile2vector") then
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read tile restart files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call ReadTileRestart(namelist, date, tile)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Transfer tile to vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    iloc = 0
    do itile = 1, 6
    do iy = 1, namelist%tile_size
    do ix = 1, namelist%tile_size
    
      if(tile%land_frac(ix,iy,itile) > 0.0) then
        iloc = iloc + 1
        vector%swe(iloc)                   = tile%swe(ix,iy,itile)
        vector%snow_depth(iloc)            = tile%snow_depth(ix,iy,itile)
        vector%active_snow_layers(iloc)    = tile%active_snow_layers(ix,iy,itile)
        vector%swe_previous(iloc)          = tile%swe_previous(ix,iy,itile)
        vector%snow_soil_interface(iloc,:) = tile%snow_soil_interface(ix,iy,:,itile)
        vector%temperature_snow(iloc,:)    = tile%temperature_snow(ix,iy,:,itile)
        vector%snow_ice_layer(iloc,:)      = tile%snow_ice_layer(ix,iy,:,itile)
        vector%snow_liq_layer(iloc,:)      = tile%snow_liq_layer(ix,iy,:,itile)
        vector%temperature_soil(iloc,:)    = tile%temperature_soil(ix,iy,:,itile)
      end if
      
    end do
    end do
    end do
      
   print*, "Transferred ",iloc, "land grids"  
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write FV3 tile file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call WriteVectorRestart(namelist, date, vector, vector_length)
  
  end if ! "vector2tile" or "tile2vector" branch
     
contains   
  
  subroutine ReadVectorRestart(namelist, date, vector, vector_length)
  
  use netcdf

  type(namelist_type) :: namelist
  type(vector_type)   :: vector
  character*19        :: date
  integer             :: vector_length
  character*256       :: vector_filename, filename
  integer             :: ncid, dimid, varid, status
  logical             :: file_exists
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create vector file name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(vector_filename,'(a17,a19,a3)') "ufs_land_restart.", date, ".nc"

  filename = trim(namelist%vector_restart_path)//trim(vector_filename)
  
  inquire(file=filename, exist=file_exists)
  
  if(.not.file_exists) then 
    print*, trim(filename), " does not exist"
    print*, "Check paths and file name"
    stop 10 
  end if
    
  print*, "Reading vector file: ", trim(filename)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check the vector length, fail if not consistent with tile-calculated length
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call ReadVectorLength(filename, vector_length)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read the vector fields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  status = nf90_open(filename, NF90_NOWRITE, ncid)

  status = nf90_inq_varid(ncid, "weasd", varid)
  status = nf90_get_var(ncid, varid , vector%swe   , &
      start = (/1,1/), count = (/vector_length, 1/))

  status = nf90_inq_varid(ncid, "snwdph", varid)
  status = nf90_get_var(ncid, varid , vector%snow_depth  , &
      start = (/1,1/), count = (/vector_length, 1/))

  status = nf90_inq_varid(ncid, "snowxy", varid)
  status = nf90_get_var(ncid, varid , vector%active_snow_layers  , &
      start = (/1,1/), count = (/vector_length, 1/))

  status = nf90_inq_varid(ncid, "sneqvoxy", varid)
  status = nf90_get_var(ncid, varid , vector%swe_previous, &
      start = (/1,1/), count = (/vector_length, 1/))

  status = nf90_inq_varid(ncid, "tsnoxy", varid)
  status = nf90_get_var(ncid, varid , vector%temperature_snow  , &
      start = (/1            , 1, 1/)                , &
      count = (/vector_length, 3, 1/))

  status = nf90_inq_varid(ncid, "zsnsoxy", varid)
  status = nf90_get_var(ncid, varid , vector%snow_soil_interface , &
      start = (/1            , 1, 1/)                , &
      count = (/vector_length, 7, 1/))

  status = nf90_inq_varid(ncid, "snicexy", varid)
  status = nf90_get_var(ncid, varid , vector%snow_ice_layer , &
      start = (/1            , 1, 1/)                , &
      count = (/vector_length, 3, 1/))

  status = nf90_inq_varid(ncid, "snliqxy", varid)
  status = nf90_get_var(ncid, varid , vector%snow_liq_layer , &
      start = (/1            , 1, 1/)                , &
      count = (/vector_length, 3, 1/))
 
  status = nf90_inq_varid(ncid, "stc", varid)
  status = nf90_get_var(ncid, varid , vector%temperature_soil , &
      start = (/1            , 1, 1/)                , &
      count = (/vector_length, 4, 1/))

  status = nf90_inq_varid(ncid, "smc", varid)
  status = nf90_get_var(ncid, varid , vector%soil_moisture , &
      start = (/1            , 1, 1/)                , &
      count = (/vector_length, 1, 1/))

  status = nf90_close(ncid)

  ! read vegetation from static file 

  filename = trim(namelist%static_filename)
  
  inquire(file=filename, exist=file_exists)
  
  if(.not.file_exists) then 
    print*, trim(filename), " does not exist"
    print*, "Check paths and file name"
    stop 10 
  end if
  
  status = nf90_open(filename, NF90_NOWRITE, ncid)

  status = nf90_inq_varid(ncid, "vegetation_category", varid)
  status = nf90_get_var(ncid, varid , vector%vegetation_type, &
      start = (/1,1/), count = (/vector_length, 1/))
      
  status = nf90_close(ncid)
  end subroutine ReadVectorRestart

  subroutine ReadTileRestart(namelist, date, tile)
  
  use netcdf

  type(namelist_type) :: namelist
  type(tile_type)     :: tile
  character*19        :: date
  character*256       :: tile_filename
  integer             :: ncid, dimid, varid, status
  logical             :: file_exists
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create tile file name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do itile = 1, 6
    write(tile_filename,'(a4,a2,a2,a1,a2,a18,i1,a3)')  & 
        date(1:4), date(6:7), date(9:10),".",date(12:13), "0000.sfc_data.tile",itile,".nc"

    tile_filename = trim(namelist%tile_restart_path)//trim(tile_filename)
    
    inquire(file=tile_filename, exist=file_exists)
  
    if(.not.file_exists) then 
      print*, trim(tile_filename), " does not exist"
      print*, "Check paths and file name"
      stop 10 
    end if
    
    print*, "Reading tile file: ", trim(tile_filename)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read the tile fields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    status = nf90_open(tile_filename, NF90_NOWRITE, ncid)
      if (status /= nf90_noerr) call handle_err(status)

! Start reading restart file
  
    status = nf90_inq_varid(ncid, "sheleg", varid)
    status = nf90_get_var(ncid, varid , tile%swe(:,:,itile)   , &
      start = (/1,1,1/), count = (/namelist%tile_size, namelist%tile_size, 1/))

    status = nf90_inq_varid(ncid, "snwdph", varid)
    status = nf90_get_var(ncid, varid , tile%snow_depth(:,:,itile)   , &
      start = (/1,1,1/), count = (/namelist%tile_size, namelist%tile_size, 1/))

    status = nf90_inq_varid(ncid, "snowxy", varid)
    status = nf90_get_var(ncid, varid , tile%active_snow_layers(:,:,itile)   , &
      start = (/1,1,1/), count = (/namelist%tile_size, namelist%tile_size, 1/))

    status = nf90_inq_varid(ncid, "sneqvoxy", varid)
    status = nf90_get_var(ncid, varid , tile%swe_previous(:,:,itile)   , &
      start = (/1,1,1/), count = (/namelist%tile_size, namelist%tile_size, 1/))

    status = nf90_inq_varid(ncid, "zsnsoxy", varid)
    status = nf90_get_var(ncid, varid , tile%snow_soil_interface(:,:,:,itile) , &
      start = (/1                , 1                , 1, 1/), &
      count = (/namelist%tile_size, namelist%tile_size, 7, 1/))

    status = nf90_inq_varid(ncid, "tsnoxy", varid)
    status = nf90_get_var(ncid, varid , tile%temperature_snow(:,:,:,itile)  , &
      start = (/1                , 1                , 1, 1/), &
      count = (/namelist%tile_size, namelist%tile_size, 3, 1/))

    status = nf90_inq_varid(ncid, "snicexy", varid)
    status = nf90_get_var(ncid, varid , tile%snow_ice_layer(:,:,:,itile) , &
      start = (/1                , 1                , 1, 1/), &
      count = (/namelist%tile_size, namelist%tile_size, 3, 1/))

    status = nf90_inq_varid(ncid, "snliqxy", varid)
    status = nf90_get_var(ncid, varid , tile%snow_liq_layer(:,:,:,itile) , &
      start = (/1                , 1                , 1, 1/), &
      count = (/namelist%tile_size, namelist%tile_size, 3, 1/))

    status = nf90_inq_varid(ncid, "stc", varid)
    status = nf90_get_var(ncid, varid , tile%temperature_soil(:,:,:,itile)   , &
      start = (/1                , 1                , 1, 1/), &
      count = (/namelist%tile_size, namelist%tile_size, 4, 1/))

    status = nf90_close(ncid)

  end do
  
  end subroutine ReadTileRestart

  subroutine WriteVectorRestart(namelist, date, vector, vector_length)
  
  use netcdf

  type(namelist_type) :: namelist
  type(vector_type)   :: vector
  character*19        :: date
  integer             :: vector_length
  character*256       :: vector_filename
  integer             :: ncid, dimid, varid, status
  logical             :: file_exists
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create vector file name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(vector_filename,'(a17,a19,a3)') "ufs_land_restart.", date, ".nc"

  vector_filename = trim(namelist%output_path)//trim(vector_filename)
  
  inquire(file=vector_filename, exist=file_exists)
  
  if(.not.file_exists) then 
    print*, trim(vector_filename), " does not exist"
    print*, "Check paths and file name"
    stop 10
  end if
    
  print*, "Writing vector file: ", trim(vector_filename)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check the vector length, fail if not consistent with tile-calculated length
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call ReadVectorLength(vector_filename, vector_length)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write the vector fields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  status = nf90_open(vector_filename, NF90_WRITE, ncid)
      if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_varid(ncid, "weasd", varid)
  status = nf90_put_var(ncid, varid , vector%swe   , &
      start = (/1,1/), count = (/vector_length, 1/))

  status = nf90_inq_varid(ncid, "snwdph", varid)
  status = nf90_put_var(ncid, varid , vector%snow_depth  , &
      start = (/1,1/), count = (/vector_length, 1/))

  status = nf90_inq_varid(ncid, "snowxy", varid)
  status = nf90_put_var(ncid, varid , vector%active_snow_layers  , &
      start = (/1,1/), count = (/vector_length, 1/))

  status = nf90_inq_varid(ncid, "sneqvoxy", varid)
  status = nf90_put_var(ncid, varid , vector%swe_previous, &
      start = (/1,1/), count = (/vector_length, 1/))

  status = nf90_inq_varid(ncid, "tsnoxy", varid)
  status = nf90_put_var(ncid, varid , vector%temperature_snow  , &
      start = (/1            , 1, 1/)                , &
      count = (/vector_length, 3, 1/))

  status = nf90_inq_varid(ncid, "zsnsoxy", varid)
  status = nf90_put_var(ncid, varid , vector%snow_soil_interface , &
      start = (/1            , 1, 1/)                , &
      count = (/vector_length, 7, 1/))

  status = nf90_inq_varid(ncid, "snicexy", varid)
  status = nf90_put_var(ncid, varid , vector%snow_ice_layer , &
      start = (/1            , 1, 1/)                , &
      count = (/vector_length, 3, 1/))

  status = nf90_inq_varid(ncid, "snliqxy", varid)
  status = nf90_put_var(ncid, varid , vector%snow_liq_layer , &
      start = (/1            , 1, 1/)                , &
      count = (/vector_length, 3, 1/))
 
  status = nf90_inq_varid(ncid, "stc", varid)
  status = nf90_put_var(ncid, varid , vector%temperature_soil , &
      start = (/1            , 1, 1/)                , &
      count = (/vector_length, 4/))

  status = nf90_close(ncid)

  end subroutine WriteVectorRestart

  subroutine WriteTileRestart(namelist, date, tile)
  
  use netcdf

  type(namelist_type) :: namelist
  type(tile_type)     :: tile
  character*19        :: date
  character*256       :: tile_filename
  integer             :: itile
  integer             :: ncid, varid, status, i
  integer             :: dim_id_xdim, dim_id_ydim, dim_id_soil, dim_id_snow, dim_id_snso, dim_id_time
  
  do itile = 1, 6

    !write(tile_filename,'(a17,a19,a5,i1,a3)') "ufs_land_restart.", date, ".tile", itile, ".nc"
    write(tile_filename,'(a4,a2,a2,a1,a2,a18,i1,a3)')  & 
        date(1:4), date(6:7), date(9:10),".",date(12:13), "0000.sfc_data.tile",itile,".nc"

    tile_filename = trim(namelist%output_path)//trim(tile_filename)
    
    print*, "Writing tile file: ", trim(tile_filename)

    status = nf90_create(tile_filename, NF90_CLOBBER, ncid)
      if (status /= nf90_noerr) call handle_err(status)

! Define dimensions in the file.

    status = nf90_def_dim(ncid, "xaxis_1"          , namelist%tile_size , dim_id_xdim)
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_dim(ncid, "yaxis_1"          , namelist%tile_size , dim_id_ydim)
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_dim(ncid, "zaxis_2"   , 4                  , dim_id_soil)
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_dim(ncid, "zaxis_3"   , 3                  , dim_id_snow)
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_dim(ncid, "zaxis_4"   , 7                  , dim_id_snso)
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_dim(ncid, "Time"          , NF90_UNLIMITED     , dim_id_time)
      if (status /= nf90_noerr) call handle_err(status)

! define dimension variables (for JEDI) 

    status = nf90_def_var(ncid, "Time", NF90_DOUBLE,    &
      (/dim_id_time/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_def_var(ncid, "xaxis_1", NF90_DOUBLE,    &
      (/dim_id_xdim/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_def_var(ncid, "yaxis_1", NF90_DOUBLE,    &
      (/dim_id_ydim/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_def_var(ncid, "zaxis_2", NF90_DOUBLE,    &
      (/dim_id_soil/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_def_var(ncid, "zaxis_3", NF90_DOUBLE,    &
      (/dim_id_snow/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_def_var(ncid, "zaxis_4", NF90_DOUBLE,    &
      (/dim_id_snso/), varid)
    if (status /= nf90_noerr) call handle_err(status)

  
! Define variables in the file.

    status = nf90_def_var(ncid, "sheleg", NF90_DOUBLE,    & ! note: this is weasd in vector file.
      (/dim_id_xdim,dim_id_ydim,dim_id_time/), varid)
      if (status /= nf90_noerr) call handle_err(status)

    status = nf90_def_var(ncid, "snwdph", NF90_DOUBLE,   &
      (/dim_id_xdim,dim_id_ydim,dim_id_time/), varid)
      if (status /= nf90_noerr) call handle_err(status)

    status = nf90_def_var(ncid, "snowxy", NF90_DOUBLE,   &
      (/dim_id_xdim,dim_id_ydim,dim_id_time/), varid)
      if (status /= nf90_noerr) call handle_err(status)

    status = nf90_def_var(ncid, "sneqvoxy", NF90_DOUBLE, &
      (/dim_id_xdim,dim_id_ydim,dim_id_time/), varid)
      if (status /= nf90_noerr) call handle_err(status)

    status = nf90_def_var(ncid, "zsnsoxy", NF90_DOUBLE,  &
      (/dim_id_xdim,dim_id_ydim,dim_id_snso,dim_id_time/), varid)
      if (status /= nf90_noerr) call handle_err(status)

    status = nf90_def_var(ncid, "tsnoxy", NF90_DOUBLE,   &
      (/dim_id_xdim,dim_id_ydim,dim_id_snow,dim_id_time/), varid)
      if (status /= nf90_noerr) call handle_err(status)

    status = nf90_def_var(ncid, "snicexy", NF90_DOUBLE,  &
      (/dim_id_xdim,dim_id_ydim,dim_id_snow,dim_id_time/), varid)
      if (status /= nf90_noerr) call handle_err(status)

    status = nf90_def_var(ncid, "snliqxy", NF90_DOUBLE,  &
      (/dim_id_xdim,dim_id_ydim,dim_id_snow,dim_id_time/), varid)
      if (status /= nf90_noerr) call handle_err(status)

    status = nf90_def_var(ncid, "stc", NF90_DOUBLE,      &
      (/dim_id_xdim,dim_id_ydim,dim_id_soil,dim_id_time/), varid)
      if (status /= nf90_noerr) call handle_err(status)

    status = nf90_def_var(ncid, "smc", NF90_DOUBLE,   &
      (/dim_id_xdim,dim_id_ydim,dim_id_time/), varid)
      if (status /= nf90_noerr) call handle_err(status)

    status = nf90_def_var(ncid, "slmsk", NF90_DOUBLE,   &
      (/dim_id_xdim,dim_id_ydim,dim_id_time/), varid)
      if (status /= nf90_noerr) call handle_err(status)
      
  status = nf90_def_var(ncid, "vtype", NF90_DOUBLE,   &
      (/dim_id_xdim,dim_id_ydim,dim_id_time/), varid)
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_enddef(ncid)

! fill dimension variables 

    status = nf90_inq_varid(ncid, "Time", varid)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, varid ,(/1/) )

    status = nf90_inq_varid(ncid, "xaxis_1", varid)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, varid ,(/(i, i=1, namelist%tile_size)/) )

    status = nf90_inq_varid(ncid, "yaxis_1", varid)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, varid ,(/(i, i=1, namelist%tile_size)/) )

    status = nf90_inq_varid(ncid, "zaxis_2", varid)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, varid ,(/(i, i=1, 4)/) )

    status = nf90_inq_varid(ncid, "zaxis_3", varid)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, varid ,(/(i, i=1, 3)/) )

    status = nf90_inq_varid(ncid, "zaxis_4", varid)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(ncid, varid ,(/(i, i=1, 7)/) )

! Start writing restart file
  
    status = nf90_inq_varid(ncid, "sheleg", varid)
    status = nf90_put_var(ncid, varid , tile%swe(:,:,itile)   , &
      start = (/1,1,1/), count = (/namelist%tile_size, namelist%tile_size, 1/))

    status = nf90_inq_varid(ncid, "snwdph", varid)
    status = nf90_put_var(ncid, varid , tile%snow_depth(:,:,itile)   , &
      start = (/1,1,1/), count = (/namelist%tile_size, namelist%tile_size, 1/))

    status = nf90_inq_varid(ncid, "snowxy", varid)
    status = nf90_put_var(ncid, varid , tile%active_snow_layers(:,:,itile)   , &
      start = (/1,1,1/), count = (/namelist%tile_size, namelist%tile_size, 1/))

    status = nf90_inq_varid(ncid, "sneqvoxy", varid)
    status = nf90_put_var(ncid, varid , tile%swe_previous(:,:,itile)   , &
      start = (/1,1,1/), count = (/namelist%tile_size, namelist%tile_size, 1/))

    status = nf90_inq_varid(ncid, "zsnsoxy", varid)
    status = nf90_put_var(ncid, varid , tile%snow_soil_interface(:,:,:,itile) , &
      start = (/1                , 1                , 1, 1/), &
      count = (/namelist%tile_size, namelist%tile_size, 7, 1/))

    status = nf90_inq_varid(ncid, "tsnoxy", varid)
    status = nf90_put_var(ncid, varid , tile%temperature_snow(:,:,:,itile)  , &
      start = (/1                , 1                , 1, 1/), &
      count = (/namelist%tile_size, namelist%tile_size, 3, 1/))

    status = nf90_inq_varid(ncid, "snicexy", varid)
    status = nf90_put_var(ncid, varid , tile%snow_ice_layer(:,:,:,itile) , &
      start = (/1                , 1                , 1, 1/), &
      count = (/namelist%tile_size, namelist%tile_size, 3, 1/))

    status = nf90_inq_varid(ncid, "snliqxy", varid)
    status = nf90_put_var(ncid, varid , tile%snow_liq_layer(:,:,:,itile) , &
      start = (/1                , 1                , 1, 1/), &
      count = (/namelist%tile_size, namelist%tile_size, 3, 1/))

    status = nf90_inq_varid(ncid, "stc", varid)
    status = nf90_put_var(ncid, varid , tile%temperature_soil(:,:,:,itile)   , &
      start = (/1,1,1/), count = (/namelist%tile_size, namelist%tile_size, 4, 1/))

    status = nf90_inq_varid(ncid, "smc", varid)
    status = nf90_put_var(ncid, varid , tile%soil_moisture(:,:,itile)   , &
      start = (/1,1,1/), count = (/namelist%tile_size, namelist%tile_size, 1/))

! include in output, so can be used to id which tile grid cells are being simulated
    status = nf90_inq_varid(ncid, "slmsk", varid)
    status = nf90_put_var(ncid, varid , tile%slmsk(:,:,itile)   , &
      start = (/1,1,1/), count = (/namelist%tile_size, namelist%tile_size, 1/))

    status = nf90_inq_varid(ncid, "vtype", varid)
    status = nf90_put_var(ncid, varid , tile%vegetation_type(:,:,itile)   , &
      start = (/1,1,1/), count = (/namelist%tile_size, namelist%tile_size, 1/))
      
  status = nf90_close(ncid)

  end do
  
  end subroutine WriteTileRestart
  
  subroutine ReadVectorLength(filename, vector_length)
  
  use netcdf

  character*256     :: filename
  integer           :: vector_length
  integer           :: length_from_file
  integer           :: ncid, dimid, varid, status
  
  status = nf90_open(filename, NF90_NOWRITE, ncid)

  status = nf90_inq_dimid(ncid, "location", dimid)
  status = nf90_inquire_dimension(ncid, dimid, len = length_from_file)
  
  status = nf90_close(ncid)
  
  if(vector_length /= length_from_file) then
    print*, "number of land points in tiles not consistent with land model vector length"
    stop 10 
  else
    print*, "number of land points in tiles consistent with land model vector length"
  end if

  end subroutine ReadVectorLength

  subroutine ReadNamelist(namelist)
  
    type(namelist_type) :: namelist    
    character*11        :: direction
    character*256       :: tile_path
    integer             :: tile_size
    character*19        :: restart_date
    character*256       :: vector_restart_path
    character*256       :: tile_restart_path
    character*256       :: output_path
    character*256       :: static_filename
  
    namelist / run_setup  / direction, tile_path, tile_size, restart_date, vector_restart_path, tile_restart_path, output_path, static_filename

    open(30, file=namelist%namelist_name, form="formatted")
     read(30, run_setup)
    close(30)
    
    namelist%direction           = direction
    namelist%tile_path           = tile_path
    namelist%tile_size           = tile_size
    namelist%restart_date        = restart_date
    namelist%vector_restart_path = vector_restart_path
    namelist%tile_restart_path   = tile_restart_path
    namelist%output_path         = output_path
    namelist%static_filename     = static_filename
  end subroutine ReadNamelist

  subroutine handle_err(status)
    use netcdf
    integer, intent ( in) :: status
 
    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 10
    end if
  end subroutine handle_err

end program vector2tile_converter
