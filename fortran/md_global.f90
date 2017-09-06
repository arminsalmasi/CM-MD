!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declerations - program global variables and parameters
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE global
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PUBLIC 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! constants for real size defginitions
  INTEGER, PARAMETER :: sp = kind(1.0),&
                        dp = selected_real_kind(2*precision(1.0_sp)), &
                        qp = selected_real_kind(2*precision(1.0_dp))
  
  REAL(dp), PARAMETER :: kb = 1.38064852e-23 ! Boltzman constant
!************************************************************************************
  ! Type to save EAM raw data from file -- AND DERIVITIVES: intent mobility of data
  TYPE EAM_data
    REAL(dp), ALLOCATABLE :: F1   (:)     &   
                           , rho1 (:)     &
                           , F2   (:)     &
                           , rho2 (:)     &
                           , phi11(:)     &
                           , phi21(:)     &
                           , phi22(:)     

    REAL(dp) :: drho, dr, cutoff       
    
    INTEGER :: nR, nRho

  END TYPE EAM_data
!************************************************************************************ 
  ! to save data of one timestep -- intent :sampling
  TYPE timestep_saver

    REAL(dp), ALLOCATABLE :: xyz_save(:,:) &! coordinate 
                           , vel_save(:,:) &! velocitie
                           , acc_save(:,:) &! acceleration
                           , frc_save(:,:) &! forces
                           , EAM_save(:)    ! EAM potential energy
    REAL(dp) :: Kin_save &! Kinetic energy
              , E_save   &! Total Energy 
              , err_save  ! relative error
  
  END TYPE timestep_saver
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE global
