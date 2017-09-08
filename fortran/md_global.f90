!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declerations - program global variables and parameters
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE glbl
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PUBLIC 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! constants for real size defginitions
  INTEGER, PARAMETER :: sp = kind(1.0),&
                        dp = selected_real_kind(2*precision(1.0_sp)), &
                        qp = selected_real_kind(2*precision(1.0_dp))
  
  REAL(dp), PARAMETER :: kb = 1.38064852e-23 ! Boltzman constant
!************************************************************************************
  ! Type to save EAM NIST Ni_Co data from file .alloy file
  TYPE EAM_data
    REAL(dp), ALLOCATABLE :: F1   (:) &   
                           , rho1 (:) &
                           , F2   (:) &
                           , rho2 (:) &
                           , phi11(:) &
                           , phi21(:) &
                           , phi22(:)     

    REAL(dp) :: drho, dr, cutoff       
    
    INTEGER :: nR, nRho

  END TYPE EAM_data
!************************************************************************************ 
  ! to save data of one timestep -- intent :sampling
  TYPE timestep_sample

    REAL(dp), ALLOCATABLE :: xyz(:,:) &! coordinate 
                           , vel(:,:) &! velocitie
                           , acc(:,:) &! acceleration
                           , frc(:,:) &! forces
                           , EAM(:)    ! EAM potential energy
    REAL(dp) :: Kin &! Kinetic energy
              , E   &! Total Energy 
              , error & ! relative error
              , time
  END TYPE timestep_sample
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE glbl
