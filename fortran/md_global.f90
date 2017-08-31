!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declerations - program global variables and parameters
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module global
  
  Public 
   
  ! constants for real size defginitions
  INTEGER, PARAMETER :: sp = kind(1.0),&
                        dp = selected_real_kind(2*precision(1.0_sp)), &
                        qp = selected_real_kind(2*precision(1.0_dp))
  
  REAL(dp), PARAMETER :: kb = 1.38064852e-23 ! Boltzman constant
  
  ! to save data of one timestep
  TYPE timestep_saver
  REAL(dp), ALLOCATABLE :: pos_save(:,:) &! coordinate 
                         , vel_save(:,:) &! velocitie
                         , acc_save(:,:) &! acceleration
                         , frc_save(:,:)  ! forces
  END TYPE timestep_saver
						
  
end module global
