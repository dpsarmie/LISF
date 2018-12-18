!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: read_petlvt
! \label{read_petlvt}
!
! !REVISION HISTORY:
!  06Sep2018: D. Sarmiento; Initial code, modeled after PET USGS code
!
! !INTERFACE:
subroutine read_petlvt (n, kk, findex, pet_filename, ferror_petlvt )

! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use LIS_logMod
  use LIS_metforcingMod, only : LIS_forc
  use petlvt_forcingMod, only : petlvt_struc
  use netcdf

  implicit none

! !ARGUMENTS:   
  integer, intent(in) :: n 
  integer, intent(in) :: kk
  integer, intent(in) :: findex
  character(99), intent(in) :: pet_filename  
  integer,intent(out) :: ferror_petlvt
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  PET LVT data and interpolates to the LIS domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the forcing dataset
!  \item[kk]
!    index of the forecast member
!  \item[pet\_filename]
!    name of the PET LVT file
!  \item[ferror\_petlvt]
!    flag to indicate success of the call (=0 indicates success)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_petlvt](\ref{interp_petlvt}) \newline
!    spatially interpolates the PET LVT data
!  \end{description}
!EOP

!==== Local Variables ==============================

!  integer, parameter :: xd=361, yd=576      ! Dimension of NOAA coarse ET global product
!  integer, parameter :: xd=359, yd=575       ! Dimension of PET LVT coarse 0.5 x 0.625 deg data
  integer, parameter :: xd=1433, yd=2871    ! Dimension of PET LVT fine 0.125 x 0.125 deg data
  real               :: realpet(yd,xd)
  integer            :: i,j
  integer            :: ftn_pet, ETId
  integer            :: index1, index2
  real, pointer      :: pet_regrid(:,:)      ! Interpolated PET array

  logical            :: file_exists
  real               :: ETos(yd,xd)

!====================  End Variable Definition  =======================

   allocate (pet_regrid(LIS_rc%lnc(n),LIS_rc%lnr(n)))
   pet_regrid = -1.0
   realpet    = -1.0
   ferror_petlvt = 1  ! Default file available
 
!----------------------------------------------------------------------
! Read PET LVT data from NetCDF files 
!----------------------------------------------------------------------

!- Inquire if LVT file exists:
   inquire( file=trim(pet_filename), exist=file_exists )
   if ( file_exists .eqv. .false. ) then
      write(LIS_logunit,*) "Missing LVT PET File, ", trim(pet_filename)
      write(LIS_logunit,*)  " Calling end run."
      call LIS_endrun()
      ferror_petlvt = 0
      return
   endif
   
   call LIS_verify(nf90_open(path=trim(pet_filename), mode=NF90_NOWRITE, &
          ncid=ftn_pet), 'nf90_open failed for PET_LVT file in read_petlvt')

!- Code to read in Ref ET product from Mike Hobbins at NOAA
!   call LIS_verify(nf90_inq_varid(ftn_pet,'ETos', ETId), &
!          'nf90_inq_varid failed for ETos in read_petlvt')

!   call LIS_verify(nf90_get_var(ftn_pet,ETId, ETos), &
!          'nf90_get_var failed for ETos in read_petlvt')

!- Code to read in Ref ET product that is produced in LVT
   call LIS_verify(nf90_inq_varid(ftn_pet,'RefET_from_RefET_v_RefET_ds1', ETId), &
          'nf90_inq_varid failed for RefET in read_petlvt')

   call LIS_verify(nf90_get_var(ftn_pet,ETId, ETos), &
          'nf90_get_var failed for RefET in read_petlvt')

!=== End of data reading

!- Data read in successfully:
   if( ferror_petlvt == 1 ) then

   !- Set any undefined points to -1 and everything else to realpet:
      index1 = 0
      do i = 1, xd         
         do j = 1, yd
            if( ETos(j,i) /= ETos(j,i) ) then !Check for NaN values
              realpet(j,i) = -1
            elseif( ETos(j,i) > 0 ) then
              realpet(j,i) = ETos(j,i)
            else
              realpet(j,i) = -1
            endif
         enddo
      enddo
  !-- Spatially interpolate original file domain to LIS grid domain:

      call interp_petlvt( n, findex, xd, yd, realpet, LIS_rc%gridDesc(n,:), &
                           LIS_rc%lnc(n), LIS_rc%lnr(n), pet_regrid )

      do j = 1,LIS_rc%lnr(n)
         do i = 1,LIS_rc%lnc(n)
            if( pet_regrid(i,j) .ne. -1.0 ) then
               index2 = LIS_domain(n)%gindex(i,j)
               if( index2 .ne. -1 ) then
                  petlvt_struc(n)%metdata2(kk,1,index2) = pet_regrid(i,j)   ! mm/day
               endif
            endif
         enddo
      enddo

      write(LIS_logunit,*) "Obtained PET LVT data:: ", pet_filename

   elseif( ferror_petlvt == 0 ) then
      write(LIS_logunit,*) "Missing PET LVT data ", pet_filename

   endif

   deallocate(pet_regrid)

end subroutine read_petlvt
