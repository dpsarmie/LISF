!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: finalize_petlvt
! \label{finalize_petlvt}
! 
! !REVISION HISTORY: 
! 06Sep2018: D. Sarmiento; Initial Specification, framed after PET_USGS
! 
! !INTERFACE:
subroutine finalize_petlvt(findex)

! !USES:
  use petlvt_forcingMod, only : petlvt_struc

! !DESCRIPTION: 
!  Routine to cleanup PET LVT forcing related memory allocations.   
!
!EOP
  implicit none
  integer, intent(IN) :: findex  


  deallocate( petlvt_struc )


end subroutine finalize_petlvt
