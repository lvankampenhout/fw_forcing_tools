program def_blanking_mask

   ! ---------------------------------------------------------------------------|
   !
   !   Function:   Finds in a NetCDF file the land cells on, and ocean cells 
   !               around a target landmass.
   !               The landmass is defined by the variables IX_START and IY_START
   !               which are hardcoded. 
   !               A width is defined for the ocean part of the mask. Any ocean cell
   !               that has distance <= WIDTH from a masked land point will be
   !               masked as well.
   !
   !               Input is a file with dimensions nlon=320 x nlat=384
   !               and variable TOPO that defines the topography.
   !  
   !               Output is the same file, with additional variable MASK (integer)
   !               MASK is
   !                       +1 :  on masked cells (land or ocean)
   !                        0 :  everywhere else
   !
   !               In case MASK already exists it will be overwritten.
   !
   !   Usage:      ifort -fpp -lnetcdff def_blanking_mask.f90
   !
   !   Date:       19-Mar-2014
   !   Author:     Leo van Kampenhout
   !
   !   Date:       31-Mar-2014
   !   Author:     Leo van Kampenhout
   !   Updates:    Add starting coordinates for Antartica, update comments.
   !
   !   Date:       04-Jun-2014
   !   Author:     Leo van Kampenhout
   !   Updates:    Adapted from find_bounding_cells.f90
   !
   !   Date:       05-Jun-2014
   !   Author:     Leo van Kampenhout
   !   Updates:    Remove canadian arctic from mask.
   !
   !   Date:       08-Aug-2014
   !   Author:     Leo van Kampenhout
   !   Updates:    Remove magic number 320.
   !
   ! ---------------------------------------------------------------------------|
 
   use netcdf
   implicit none
 
   ! *** PARAMETERS

   ! This is the name of the data file we will read. 
   character (len = *), parameter :: TOPO_NAME = "TOPO"
   character (len = *), parameter :: MASK_NAME = "BLK_MASK"
   character (len = *), parameter :: NLON_S = "nlon"
   character (len = *), parameter :: NLAT_S = "nlat"

   ! How many cells wide is the masked (ocean) region
   integer, parameter :: WIDTH = 10

   character (len = *), parameter :: FILE_NAME = "blanking_mask/regions.nc"
   integer, parameter :: IX_START = 311, IY_START = 369 ! Starting point inside landmass (Greenland)
   integer, parameter :: IX_START2 = 40, IY_START2 = 1 ! Starting point inside landmass (Antartica)
 
   ! *** LOCALS

   ! dimension ids of main dimensions
   integer :: idx
   integer :: idy
   integer :: idt

   ! We are reading 3D data: x, y, t
   integer :: nlon, nlat, nt
   integer, allocatable :: data_in(:,:) ! Only X/Y
   integer, allocatable :: data_out(:,:) ! Only X/Y
 
   ! This will be the netCDF ID for the file and data variable.
   integer :: ncid, varid, dimid
 
   ! Loop indexes, and error handling.
   integer :: ix, iy, it, ix2, iy2
   integer :: dx, dy
   
   integer :: count(2), start(2), dimids(2)
   real :: runmin, runmax
 
   ! ---------------------------------------------------------------------------|
   ! *** Start of executable code *** 
   ! ---------------------------------------------------------------------------|


   ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
   ! the file.
   call check( nf90_open(FILE_NAME, NF90_WRITE, ncid) )
 
   call check( nf90_inq_dimid(ncid, NLON_S, dimid) )
   call check( nf90_inquire_dimension(ncid, dimid, len=nlon))
   idx = dimid

   call check( nf90_inq_dimid(ncid, NLAT_S, dimid) )
   call check( nf90_inquire_dimension(ncid, dimid, len=nlat))
   idy = dimid

   print *,'FOUND DIMENSIONS ', NLON_S,', ', NLAT_S, ' = ', nlon, nlat
   
   allocate(data_in(nlon,nlat))
   allocate(data_out(nlon,nlat))
   count = (/ nlon, nlat /)
   start = (/  1,  1 /)
 
   ! Get the varid of the data variable, based on its name.
   call check( nf90_inq_varid(ncid, TOPO_NAME, varid) )
 
   ! Check the data.
   call check( nf90_get_var(ncid, varid, data_in, start=start, count=count ))

   ! Now data_in contains the field TOPO
   print *,'LANDMASS STARTING POINT: ', IX_START, IY_START
   print *, 'VALUE OF TOPO AT THIS POINT = ', data_in(IX_START,IY_START)
   print *,'LANDMASS STARTING POINT: ', IX_START2, IY_START2
   print *, 'VALUE OF TOPO AT THIS POINT = ', data_in(IX_START2,IY_START2)

   if (data_in(IX_START,IY_START) /= 1) then
      stop 'ERROR: starting point not in landmass'
   endif

   if (data_in(IX_START2,IY_START2) /= 1) then
      stop 'ERROR: starting point not in landmass'
   endif


   if (nf90_inq_varid(ncid, MASK_NAME, varid) /= nf90_noerr ) then

      print *, 'VARIABLE ', MASK_NAME,' DOES NOT EXIST, WILL CREATE'
      call check( nf90_redef(ncid) )
      dimids = (/ idx, idy /)
      call check( nf90_def_var(ncid, MASK_NAME, NF90_DOUBLE, dimids, varid) )
      call check( nf90_enddef(ncid) )

   else

      print *, 'VARIABLE ', MASK_NAME,' EXISTS, WILL OVERWRITE'
   endif


   call check( nf90_redef(ncid) )
   call check( nf90_put_att(ncid, varid, "history", "written by program: def_blanking_mask"))
   call check( nf90_put_att(ncid, varid, "long_name", "blanking mask"))
   call check( nf90_enddef(ncid) )


   ! *** START of search
   data_out = 0

   ! Figure out land points
   call explore_point_lnd(IX_START, IY_START)
   call explore_point_lnd(IX_START2, IY_START2)

   ! Then, figure out coastal ocean points around those
   do iy = 1, nlat
   do ix = 1, nlon
      if (data_in(ix,iy) == 0) then ! ocean
         ! Loop mask land cells
         do iy2 = 1, nlat
         do ix2 = 1, nlon
            if (data_out(ix2,iy2) == 2) then
               ! Test distance to land mass
               dy = abs(iy2 - iy)

               ! also wrap around in longitude dimension
               dx = min( abs(ix2 - ix), nlon - abs(ix2 - ix))
               if ( dx + dy <= WIDTH) then
                  ! Add ocean point to mask
                  data_out(ix,iy) = 1
                  exit
               endif
            endif
         enddo
         if (data_out(ix,iy) == 1) exit
         enddo
      endif
   enddo
   enddo

   ! Manually Exclude Canadian arctic from mask... (This is kind of a hack to
   ! prevent the Canadian runoff to be wiped away since it is not replaced like
   ! greenland and artartica)
   data_out(230:236,364:373) = 0

   ! *** END of search

   ! Insert data values in file
   call check( nf90_put_var(ncid, varid, data_out) )

   ! Close the file, freeing all resources.
   call check( nf90_close(ncid) )
 
   print *,"*** SUCCESS, RESULT WRITTEN TO ", FILE_NAME
 
contains

   recursive subroutine explore_point_lnd(ix, iy)
      ! For a given point, add land points to MASK
      ! Also, call nearby land points that haven't been called already
      integer, intent(in) :: ix, iy 

      integer :: ixx, iyy, np

      ! ocean cell: do nothing for now
      if (data_in(ix,iy) == 0) &
         return

      ! already masked: do nothing
      if (data_out(ix,iy) == 2) &
         return

      ! Add current cell to mask
      data_out(ix,iy) = 2
      
      neighbours: do np = 1, 4

         select case(np)
            case (1)
               ixx = ix + 1
               iyy = iy
            case (2)
               ixx = ix - 1
               iyy = iy
            case (3)
               ixx = ix
               iyy = iy + 1
            case (4)
               ixx = ix
               iyy = iy - 1
            case default
               stop 'case select np programming error'
         end select 
         
         ! Only process points that are on the grid
         if (ixx > 0 .and. ixx < nlon+1 .and.  &
             iyy > 0 .and. iyy < nlat+1)       & 
            call explore_point_lnd(ixx,iyy)

      enddo neighbours

   end subroutine explore_point_lnd


   subroutine check(status)
      integer, intent ( in) :: status
     
      if(status /= nf90_noerr) then 
         print *, "ERROR: ", trim(nf90_strerror(status))
         stop "Stopped"
      end if
   end subroutine check  
end program def_blanking_mask

