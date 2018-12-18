!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: timeinterp_petlvt
! \label{timeinterp_petlvt}
! 
! !REVISION HISTORY: 
!  06Sep2017:  D. Sarmiento: Initial Code, modeled after PET USGS 
!
! !INTERFACE:
subroutine timeinterp_petlvt( n, findex )

! !USES:
  use ESMF
  use LIS_coreMod,        only : LIS_rc,LIS_domain
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod,  only : LIS_FORC_Base_State, LIS_forc
  use LIS_logMod,         only : LIS_verify, LIS_logunit
  use petlvt_forcingMod,  only : petlvt_struc
  use LIS_forecastMod,    only : LIS_get_iteration_index

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model 
!   timestep. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP

   integer          :: t, index1
   integer          :: mfactor, m, k, kk
   integer          :: status
   type(ESMF_Field) :: petfield
   real, pointer    :: pet(:)

! ----------------------------------------
    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_PET%varname(1)),petfield,&
          rc=status)
    call LIS_verify(status, 'Error: enable PET in forcing variables list')
    call ESMF_FieldGet(petfield,localDE=0,farrayptr=pet,rc=status)
    call LIS_verify(status)

!------------------------------------------------------------------------
   mfactor = LIS_rc%nensem(n)/petlvt_struc(n)%nIter

   do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
       t = m + (k-1)*mfactor
       index1 = LIS_domain(n)%tile(t)%index
       kk = LIS_get_iteration_index(n, k, index1, mfactor)
       if( petlvt_struc(n)%metdata2(kk,1,index1) .ne. -1.0 ) &
          pet(t) = petlvt_struc(n)%metdata2(kk,1,index1)     
       !- Convert mm/day to rate with seconds:
          pet(t) = pet(t) / 86400.     ! mm/s
      enddo
    enddo
end subroutine timeinterp_petlvt
  
