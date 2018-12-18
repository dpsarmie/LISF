!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readpetlvtcrd
! \label{readpetlvtcrd}
!
! !REVISION HISTORY:
! 06Sep2018; D. Sarmiento, Initial Code, modeled after PET USGS code
!
! !INTERFACE:    
subroutine readpetlvtcrd()

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_config, LIS_rc
  use LIS_logMod,  only : LIS_logunit, LIS_verify, LIS_endrun
  use petlvt_forcingMod, only : petlvt_struc
!
! !DESCRIPTION:
!  This routine reads the options specific to PET LVT forcing from 
!  the LIS configuration file. 
!EOP
  implicit none
  integer :: n, rc

  call ESMF_ConfigFindLabel(LIS_config,"LVT PET forcing directory:",rc=rc)
  call LIS_verify(rc, 'LVT PET forcing directory: not defined')
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,petlvt_struc(n)%petdir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"LVT PET forcing type:",rc=rc)
  call LIS_verify(rc, 'LVT PET forcing type: not defined')
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,petlvt_struc(n)%pettype,rc=rc)
     if( trim(petlvt_struc(n)%pettype) .ne. "current" .and. &
         trim(petlvt_struc(n)%pettype) .ne. "climatology" ) then
        write(LIS_logunit,*) "ERR: This LVT PET forcing type is not supported: ",&
                              petlvt_struc(n)%pettype
        write(LIS_logunit,*) " ... Please select: 'current' or 'climatology'."
        write(LIS_logunit,*) " LIS end run being called ... "
        call LIS_endrun
     endif
  enddo

  do n=1, LIS_rc%nnest
     
    write(LIS_logunit,*) "Using LVT PET forcing"
    write(LIS_logunit,*) "LVT PET forcing directory :: ",petlvt_struc(n)%petdir
    write(LIS_logunit,*) "LVT PET forcing file type :: ",petlvt_struc(n)%pettype
    write(LIS_logunit,*) " "

!------------------------------------------------------------------------
! Setting global observed PET times to zero to ensure 
! data is read in during first time step
!------------------------------------------------------------------------
    petlvt_struc(n)%pettime = 0.0

  enddo

end subroutine readpetlvtcrd
