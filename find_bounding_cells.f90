program find_bounding_cells

   ! ---------------------------------------------------------------------------|
   !
   !   Function:   Finds in a NetCDF file the coastal points around a target landmass.
   !               This landmass is determined by the variables IX_START and IY_START
   !               which are hardcoded.
   !
   !               Input is a file with dimensions nlon=320 x nlat=384
   !               and variable TOPO that defines the topography.
   !  
   !               Output is the same file, with additional variable COAST (integer)
   !               COAST is
   !                       -1 :  outside the target landmass, ocean or land
   !                        0 :  inside the landmass,
   !                       +1 :  on ocean cells directly adjacent the landmass
   !
   !               In case COAST is already existing it will be overwritten.
   !
   !   Usage:      ifort -fpp -lnetcdff find_bounding_cells.f90
   !
   !   Date:       19-Mar-2014
   !   Author:     Leo van Kampenhout
   !
   !   Date:       31-Mar-2014
   !   Author:     Leo van Kampenhout
   !   Updates:    Add starting coordinates for Antartica, update comments.
   !
   ! ---------------------------------------------------------------------------|
 
   use netcdf
   implicit none
 
   ! *** PARAMETERS

   ! This is the name of the data file we will read. 
   character (len = *), parameter :: TOPO_NAME = "TOPO"
   character (len = *), parameter :: COAST_NAME = "COAST"
   character (len = *), parameter :: XDIM = "nlon"
   character (len = *), parameter :: YDIM = "nlat"

#define ANTARTICA

#ifdef GREENLAND
   character (len = *), parameter :: FILE_NAME = "antartica/regions.nc"
   integer, parameter :: IX_START = 311, IY_START = 369 ! Starting point inside landmass (Greenland)
#endif

#ifdef ANTARTICA
   character (len = *), parameter :: FILE_NAME = "antartica/regions.nc"
   integer, parameter :: IX_START = 40, IY_START = 1 ! Starting point inside landmass (Antartica)
#endif
 
   ! *** LOCALS

   ! dimension ids of main dimensions
   integer :: idx
   integer :: idy
   integer :: idt

   ! We are reading 3D data: x, y, t
   integer :: nx, ny, nt
   integer, allocatable :: data_in(:,:) ! Only X/Y
   integer, allocatable :: data_out(:,:) ! Only X/Y
 
   ! This will be the netCDF ID for the file and data variable.
   integer :: ncid, varid, dimid
 
   ! Loop indexes, and error handling.
   integer :: ix, iy, it
   
   integer :: count(2), start(2), dimids(2)
   real :: runmin, runmax
 
   ! ---------------------------------------------------------------------------|
   ! *** Start of executable code *** 
   ! ---------------------------------------------------------------------------|


   ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
   ! the file.
   call check( nf90_open(FILE_NAME, NF90_WRITE, ncid) )
 
   call check( nf90_inq_dimid(ncid, XDIM, dimid) )
   call check( nf90_inquire_dimension(ncid, dimid, len=nx))
   idx = dimid

   call check( nf90_inq_dimid(ncid, YDIM, dimid) )
   call check( nf90_inquire_dimension(ncid, dimid, len=ny))
   idy = dimid

   print *,'FOUND DIMENSIONS ',XDIM,', ', YDIM, ' = ', nx, ny
   
   allocate(data_in(nx,ny))
   allocate(data_out(nx,ny))
   count = (/ nx, ny /)
   start = (/  1,  1 /)
 
   ! Get the varid of the data variable, based on its name.
   call check( nf90_inq_varid(ncid, TOPO_NAME, varid) )
 
   ! Check the data.
   call check( nf90_get_var(ncid, varid, data_in, start=start, count=count ))

   ! Now data_in contains the field TOPO
   print *,'LANDMASS STARTING POINT: ', IX_START, IY_START
   print *, 'VALUE OF TOPO AT THIS POINT = ', data_in(IX_START,IY_START)

   if (data_in(IX_START,IY_START) /= 1) then
      stop 'ERROR: starting point not in landmass'
   endif


   if (nf90_inq_varid(ncid, COAST_NAME, varid) /= nf90_noerr ) then

      print *, 'VARIABLE ', COAST_NAME,' DOES NOT EXIST, WILL CREATE'
      call check( nf90_redef(ncid) )
      dimids = (/ idx, idy /)
      call check( nf90_def_var(ncid, COAST_NAME, NF90_DOUBLE, dimids, varid) )
      call check( nf90_enddef(ncid) )

   else

      print *, 'VARIABLE ', COAST_NAME,' EXISTS, WILL OVERWRITE'
   endif


   call check( nf90_redef(ncid) )
   call check( nf90_put_att(ncid, varid, "history", "written by program: find_bounding_cells"))
   call check( nf90_put_att(ncid, varid, "long_name", "Greenland coastal cells"))
   call check( nf90_enddef(ncid) )


   ! *** START of search
   data_out = -1

   call explore_point(IX_START, IY_START)

   ! *** END of search

   ! Insert data values in file
   call check( nf90_put_var(ncid, varid, data_out) )

   ! Close the file, freeing all resources.
   call check( nf90_close(ncid) )
 
   print *,"*** SUCCESS, RESULT WRITTEN TO ", FILE_NAME
 
contains

   recursive subroutine explore_point(ix, iy)
      ! For a given point, add coastal points to COAST
      ! Also, call nearby land points that haven't been called already
      integer, intent(in) :: ix, iy 

      integer :: ixx, iyy, np

      if (data_in(ix,iy) == 0) then
         ! coastal point
         data_out(ix,iy) = 1
         return
      else
         ! land point
         if (data_out(ix,iy) == 0) then
            ! point was already handled 
            return
         endif
      endif

      ! Mark current cell as handled by algorithm
      data_out(ix,iy) = 0
      
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
         if (ixx > 0 .and. ixx < nx+1 .and.  &
             iyy > 0 .and. iyy < ny+1)       & 
            call explore_point(ixx,iyy)

      enddo neighbours

   end subroutine explore_point


   subroutine check(status)
      integer, intent ( in) :: status
     
      if(status /= nf90_noerr) then 
         print *, "ERROR: ", trim(nf90_strerror(status))
         stop "Stopped"
      end if
   end subroutine check  
end program find_bounding_cells

