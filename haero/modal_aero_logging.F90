module modal_aero_logging

!> @file
!> Responsible for managing logging for the aerosol model's output log
!> 
!> @author
!> Hui Wan, Pacific Northwest National Laboratory, Richland, WA
!> 

   implicit none

   public  iulog_main
   public  set_modal_aero_iulog_main

!> logical unit of the main output log

   integer :: iulog_main = 6     ! 6 = stdout (screen)

!> @todo There seems to be additional logical units declared and used 
!> in the aerosol model, for example in modal_aero_amicphys.F90
!> If they are needed, it would be better to collection the 
!> declarations and move them all to this module. 

   contains

!> @brief Sets iulog_main to a new value inew passed in from the calling routine. 
!> It is preferred that elsewhere in the aerosol model or in the GCM,
!> this subroutine is used instead of "iulog(_main) = inew",
!> so as to clarify which logical unit is being set.

   subroutine set_modal_aero_iulog_main( inew )

     integer, intent(in) :: inew

     iulog_main = inew

   end subroutine set_modal_aero_iulog_main

end module modal_aero_logging
