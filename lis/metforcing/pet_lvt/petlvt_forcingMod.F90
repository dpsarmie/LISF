!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module petlvt_forcingMod

!BOP
! !MODULE: petlvt_forcingMod
! 
! !REVISION HISTORY:
!  05Sep2018: D. Sarmiento;  Added LVT-derived PET to LIS (using petusgs code as sample)
!
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the potential evapotranspiration (PET) data 
!  from LVT output. The inital iteration of this reader uses 
!  the National Oceanic and Atmospheric Administration (NOAA) method 
!  for estimating PET from MERRA2 files and fields.  The files are global 
!  and summed daily, reflecting a 00Z timestamp. 
!
!  The implementation in LIS has the derived data type {\tt petlvt\_struc}
!  that includes the variables to specify the runtime options, and 
!  the weights and neighbor information for spatial interpolation
! 
!  They are desribed below: 
! \begin{description}
!  \item[petdir]
!    Directory containing the input data
!  \item[pettype]
!    Option to select PET file type: climatology or current (retrospective)
!  \item[pettime]
!    The nearest (daily) instance of the incoming 
!    data (as a real time).
!  \item[griduptime1]
!    The time to switch the input resolution to T170
!  \item[mi]
!    Number of points in the input grid
!  \item[n112,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for conservative interpolation. 
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for conservative interpolation.
!  \end{description}
!
! !USES: 
  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_petlvt     !defines the native resolution of 
                                    !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: petlvt_struc
!EOP

  type, public :: petlvt_type_dec

     character*100     :: petdir      ! PET LVT Forcing Directory
     character*20      :: pettype     ! PET LVT File type (climatology|current)
     real*8            :: pettime
     real*8            :: griduptime1
     logical           :: gridchange1
     integer           :: mi
     integer, allocatable :: n111(:,:)
     integer, allocatable :: n121(:,:)
     integer, allocatable :: n211(:,:)
     integer, allocatable :: n221(:,:)
     real, allocatable    :: w111(:,:),w121(:,:)
     real, allocatable    :: w211(:,:),w221(:,:)

     integer           :: nIter, st_iterid,en_iterid  ! Forecast parameters

     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:) 

  end type petlvt_type_dec

  type(petlvt_type_dec), allocatable :: petlvt_struc(:)
contains
  
!BOP
!
! !ROUTINE: init_petlvt
!
! !REVISION HISTORY: 
! 05Sep2018: D. Sarmiento;  Initial Specification, similar code as PET_USGS
! 
! !INTERFACE:
  subroutine init_petlvt(findex)

! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
    use LIS_logMod,     only : LIS_logunit, LIS_endrun
    use LIS_forecastMod

    implicit none
    integer, intent(in) :: findex

    integer :: n 
    real    :: gridDesci(50)
    integer :: updoy,r1,mo1,da1,hr1,mn1,ss1
    real    :: upgmt

!-----------------------------
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for PET LVT
!  data. The grid description arrays are based on MERRA2 and
!  followed in the LIS interpolation schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readpetlvtcrd](\ref{readpetlvtcrd}) \newline
!     reads the runtime options specified for LVT-derived PET data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!
!EOP
    allocate(petlvt_struc(LIS_rc%nnest))
    call readpetlvtcrd()

    LIS_rc%met_nf(findex) = 1  ! number of met variables in PET LVT forcing

    do n=1,LIS_rc%nnest   

       ! Forecast mode:
       if(LIS_rc%forecastMode.eq.1) then
         if(mod(LIS_rc%nensem(n),&
              LIS_forecast_struc(1)%niterations).ne.0) then
            write(LIS_logunit,*) '[ERR] The number of ensembles must be a multiple'
            write(LIS_logunit,*) '[ERR] of the number of iterations '
            write(LIS_logunit,*) '[ERR] nensem = ',LIS_rc%nensem(n)
            write(LIS_logunit,*) '[ERR] niter = ',LIS_forecast_struc(1)%niterations
            call LIS_endrun()
         endif

         petlvt_struc(n)%st_iterid = LIS_forecast_struc(1)%st_iterId
         petlvt_struc(n)%en_iterId = LIS_forecast_struc(1)%niterations
         petlvt_struc(n)%nIter = LIS_forecast_struc(1)%niterations

         allocate(petlvt_struc(n)%metdata1(LIS_forecast_struc(1)%niterations,&
              LIS_rc%met_nf(findex),&
              LIS_rc%ngrid(n)))
         allocate(petlvt_struc(n)%metdata2(LIS_forecast_struc(1)%niterations,&
              LIS_rc%met_nf(findex),&
              LIS_rc%ngrid(n)))

       ! Regular retrospective or non-forecast mode:
       else 
         petlvt_struc(n)%st_iterid = 1
         petlvt_struc(n)%en_iterId = 1
         petlvt_struc(n)%nIter = 1

         allocate(petlvt_struc(n)%metdata1(1,&
              LIS_rc%met_nf(findex),&
              LIS_rc%ngrid(n)))
         allocate(petlvt_struc(n)%metdata2(1,&
              LIS_rc%met_nf(findex),&
              LIS_rc%ngrid(n)))
       endif
    
       petlvt_struc(n)%metdata1 = 0
       petlvt_struc(n)%metdata2 = 0

    !- MERRA2 scale parameters 
!       gridDesci     = 0
!       gridDesci(1)  = 0
!       gridDesci(2)  = 575 !number of columns
!       gridDesci(3)  = 359 !number of rows
!       gridDesci(4)  = -89.50 !lower-left latitude
!       gridDesci(5)  = -179.375 !lower-left longitude
!       gridDesci(6)  = 128
!       gridDesci(7)  = 89.50 !upper right latitude
!       gridDesci(8)  = 179.325 !upper right longitude
!       gridDesci(9)  = 0.625 !resolution E-W direction
!       gridDesci(10) = 0.50 !resolution N-S direction
!       gridDesci(20) = 0

    !- Fine-scale (0.125 x 0.125deg) scale parameters
       gridDesci     = 0
       gridDesci(1)  = 0
       gridDesci(2)  = 2871 !number of columns
       gridDesci(3)  = 1433 !number of rows
       gridDesci(4)  = -89.50 !lower-left latitude
       gridDesci(5)  = -179.375 !lower-left longitude
       gridDesci(6)  = 128
       gridDesci(7)  = 89.50 !upper right latitude
       gridDesci(8)  = 179.325 !upper right longitude
       gridDesci(9)  = 0.125 !resolution E-W direction
       gridDesci(10) = 0.125 !resolution N-S direction
       gridDesci(20) = 0

    !- Timestep call:
       call LIS_update_timestep( LIS_rc, n, 86400.0 )

    !- Setting up weights for spatial interpolation:

       allocate(petlvt_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(petlvt_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(petlvt_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(petlvt_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(petlvt_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(petlvt_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(petlvt_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(petlvt_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       
       if( LIS_rc%met_interp(findex) == "bilinear" ) then

          call bilinear_interp_input( n,gridDesci,&
             petlvt_struc(n)%n111,petlvt_struc(n)%n121, &
             petlvt_struc(n)%n211,petlvt_struc(n)%n221, &
             petlvt_struc(n)%w111,petlvt_struc(n)%w121, &
             petlvt_struc(n)%w211,petlvt_struc(n)%w221 )

       elseif( LIS_rc%met_interp(findex) == "budget-bilinear" ) then

          call conserv_interp_input( n,gridDesci,&
             petlvt_struc(n)%n111,petlvt_struc(n)%n121,petlvt_struc(n)%n211,&
             petlvt_struc(n)%n221,petlvt_struc(n)%w111,petlvt_struc(n)%w121,&
             petlvt_struc(n)%w211,petlvt_struc(n)%w221 )


       elseif( LIS_rc%met_interp(findex) == "neighbor" ) then

          call neighbor_interp_input(n, gridDesci,&
             petlvt_struc(n)%n111 )

       else
          write(LIS_logunit,*) 'This interp. option not supported for PET LVT'
          write(LIS_logunit,*) 'Program stopping ... '
          call LIS_endrun()

       endif 
    enddo
  end subroutine init_petlvt
end module petlvt_forcingMod
