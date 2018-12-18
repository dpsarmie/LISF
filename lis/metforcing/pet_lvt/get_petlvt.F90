!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_petlvt
! \label{get_petlvt}
!
!
! !REVISION HISTORY:
! 06Sep2018: D. Sarmiento; Initial code, modeled after PET USGS code
!
! !INTERFACE:
subroutine get_petlvt(n, findex)

! !USES:
  use LIS_coreMod,    only : LIS_rc
  use LIS_timeMgrMod, only : LIS_tick, LIS_get_nstep
  use LIS_logMod,     only : LIS_logunit
  use petlvt_forcingMod, only : petlvt_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex

!  
! !DESCRIPTION:
!  Opens, reads, and interpolates daily LVT-generated Ref ET forcing. 
!  At the beginning of a simulation, the code reads the 
!  most recent ``past'' data file (current 1-day), and nearest 
!  ``future'' data file (next 1-day). These two datasets are used to 
!  temporally interpolate the data to the current model timestep. 
!.
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!     determines the PET LVT data times
!  \item[petlvtfile](\ref{petlvtfile}) \newline
!     Puts together appropriate file name daily intervals
!  \item[read\_petlvt](\ref{read_petlvt}) \newline
!     Interpolate PET LVT data to LIS grid
!  \end{description}
!EOP

  integer :: ferror_pet                  ! Error flags for PETdata sources
  integer :: endtime_petlvt              ! 1=get a new file
  real*8  :: ctime, ftime_petlvt         ! Current time and end boundary times for PET data sources 
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1   ! Time parameters for current LIS time
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2   ! Time parameters for file boundary end time
  real    :: ts1, ts2
  real    :: gmt1, gmt2
  integer :: kk                          ! Forecast index

  character(99) :: filename              ! Filename variables for PET data sources

!=== End Variable Definition =======================

!------------------------------------------------------------------------
! Determine required observed PET LVT data times 
! (current, accumulation end time)
! - Model current time
!------------------------------------------------------------------------
  yr1 = LIS_rc%yr    !current time
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = LIS_rc%mn
  ss1 = 0
  ts1 = 0
  call LIS_tick( ctime, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )

!------------------------------------------------------------------------
! LVT PET end time
!------------------------------------------------------------------------
  yr2 = LIS_rc%yr    !end accumulation time data
  mo2 = LIS_rc%mo
  da2 = LIS_rc%da
  hr2 = 0
  mn2 = 0
  ss2 = 0
  ts2 = ((23*60)+59)*60  ! 23 hrs + 59 mins --> Seconds
  call LIS_tick( ftime_petlvt, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, ts2 )

!------------------------------------------------------------------------
! Ensure that data is found during the first time step
!------------------------------------------------------------------------

   endtime_petlvt = 0

   if( LIS_get_nstep(LIS_rc,n) == 1 .or. LIS_rc%rstflag(n) == 1 ) then
      endtime_petlvt = 1
      LIS_rc%rstflag(n) = 0
   endif

!------------------------------------------------------------------------
! Check for and get PET LVT data file
!------------------------------------------------------------------------
    if( ctime > petlvt_struc(n)%pettime )  endtime_petlvt = 1
   
    if( endtime_petlvt == 1 ) then  ! get new time2 daily data file

      do kk= petlvt_struc(n)%st_iterid, petlvt_struc(n)%en_iterid

        ! Form PET LVT filename:
        call petlvtfile( n, kk, findex, filename, petlvt_struc(n)%pettype,&
                          petlvt_struc(n)%petdir, yr1, mo1, da1 )

        write(LIS_logunit,*) "Getting new PET LVT data (new 24 hour period): "
        write(LIS_logunit,*)  trim(filename)

        ! Read in and spatially interpolate latest PET LVT file obtained:
        ferror_pet = 0
        call read_petlvt( n, kk, findex, filename, ferror_pet )

        petlvt_struc(n)%pettime = ftime_petlvt

      end do  ! end forecast loop
    endif     ! need new time2

end subroutine get_petlvt

! --------------------------

!BOP
! !ROUTINE: petlvtfile
! \label{petlvtfile}
!
!
! !INTERFACE:
subroutine petlvtfile( n, kk, findex, filename, filetype, &
                        petlvtdir, yr, mo, da )

  use LIS_coreMod
  use LIS_forecastMod

  implicit none
! !ARGUMENTS: 
  integer         , intent(in)  :: n
  integer         , intent(in)  :: kk
  integer         , intent(in)  :: findex
  character(len=*), intent(out) :: filename
  character(len=*), intent(in)  :: petlvtdir
  character(len=*), intent(in)  :: filetype
  integer         , intent(in)  :: yr, mo, da

! !DESCRIPTION:
!   This subroutine puts together daily PET LVT filenames.
! 
!  The arguments are:
!  \begin{description}
!  \item[petlvtdir]
!    Name of the PET LVT directory
!  \item[yr]
!    year 
!  \item[mo]
!    month
!  \item[da]
!    day of month
!  \item[filename]
!    name of the timestamped PET LVT file
!  \end{description}
!
!EOP

   character(4) :: fyr4
   character(2) :: fyr2
   character(2) :: fmo, fda, fhr

!=== end variable definition =============================================

   if(LIS_rc%forecastMode.eq.0) then !hindcast run

     write(unit=fyr4,fmt='(i4.4)') yr
     write(unit=fyr2,fmt='(i2.2)') yr-2000
     write(unit=fmo, fmt='(i2.2)') mo
     write(unit=fda, fmt='(i2.2)') da

     if( filetype == "current" ) then 
!       filename = trim(petlvtdir)//"/"//fyr4//"/ETos_"//fyr4//fmo//fda//".nc"    !- Mike Hobbins ET data
       filename = trim(petlvtdir)//"/"//fyr4//"/MEAN_TS."//fyr4//fmo//fda//"0000.d01.nc" 
     elseif( filetype == "climatology" ) then
       filename = trim(petlvtdir)//"/LVT_MEAN_FINAL."//fyr4//fmo//fda//"0000.d01.nc"
     endif

   ! Forecast mode (e.g., ESP):
   else
     call LIS_sample_forecastDate(n,kk,findex,yr,mo,da)

     write(unit=fyr4,fmt='(i4.4)') yr
     write(unit=fyr2,fmt='(i2.2)') yr-2000
     write(unit=fmo, fmt='(i2.2)') mo
     write(unit=fda, fmt='(i2.2)') da

     if( filetype == "current" ) then
!       filename = trim(petlvtdir)//"/"//fyr4//"/ETos_"//fyr4//fmo//fda//".nc" !- Mike Hobbins ET data
       filename = trim(petlvtdir)//"/"//fyr4//"/MEAN_TS."//fyr4//fmo//fda//"0000.d01.nc"
     elseif( filetype == "climatology" ) then
       filename = trim(petlvtdir)//"/LVT_MEAN_FINAL."//fyr4//fmo//fda//"0000.d01.nc"
     endif
   endif

end subroutine petlvtfile
