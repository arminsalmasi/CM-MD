! Program MD
! MD is a fortran code for a simple MD simulation.
! It is possible to use eather Lennard-Jones or &&
! Embedded Atom Method (EAM) potentials.
!
! case study is a monoatomic metal in this case Ti
! 
! The code is developed by Armin Salmasi at KTH, Stockholm, Sweden.
! First code commit (commit#2-2): 2017-07-11
! 
! Feel free to use any part of this code.

PROGRAM molecular_dynamics
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Anounces
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  USE monitor
  USE global
  USE utilities
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declerations - main local variables
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  IMPLICIT NONE
  REAL(dp) :: t, dt, temperature, volume              &! simulation, time timestep, volume
            , samp_interval         &! time between saving data (purpes = speed up)
            , temperature_temp      &! 
            , E0                    &! Total energy initial
            , Kin                   &! Kinetic energy
            , E                      ! Total Energy & , A(4,3) &

  REAL(dp), ALLOCATABLE :: AtMas(:)  &! atomic mass of species 
                         , xyz(:,:)          &! coordinate 
                         , vel(:,:)          &! velocitie
                         , acc(:,:)          &! acceleration
                         , frc(:,:)          &! forces
                         , EAM(:)         

  INTEGER :: nPart       &! total number of atoms in the system
           , nTstps      &! number of timesteps
           , nSpcis      &! number of species in the system
           , nSamps      &! number of samples
           , sampOffset  &! number of timesteps between two samples except for the last one
           , sampStp     &! tstp of sampling
           , sampK		 &! sample counter  
           , tstp ,ix, j         ! counters
  INTEGER, ALLOCATABLE :: AtNums(:)     &! atomic numbers of each species A1, A2, ... 
                        , atNumsVec(:)  &!  A1,A1,A1,....; A2,A2,A2... 
                        , n_atoms(:)         ! numbers of atoms of each species n1, n2 ,..
  TYPE(timestep_saver), DIMENSION(:), ALLOCATABLE :: timestep_samples
  TYPE(EAM_data) :: EAMdata, EAMdiff         

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! MAIN
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  WRITE(*,*) 'This is my simple MD code, implemneted with fortran 90'
  ! Read from input: all vectores atomic number ordered, ascending   
  ! Lenght unit: angestrom
  CALL get_input(t, dt, temperature, volume, nTstps, samp_interval, &
                 nPart, nSpcis, n_atoms, AtNums, atNumsVec, AtMas)
  ! allocate and initialize saver
  sampOffset = Floor(samp_interval / dt * 1.0)
  nSamps = 1 + FLOOR(nTstps / sampOffset * 1.0)
  ALLOCATE(timestep_samples(0:nSamps))                           
  DO ix = 0, nSamps
    ALLOCATE(timestep_samples(ix)%xyz_save(nPart,3)) ! N*3 matrix
    ALLOCATE(timestep_samples(ix)%vel_save(nPart,3)) ! N*3 matrix
    ALLOCATE(timestep_samples(ix)%acc_save(nPart,3)) ! N*3 matrix
    ALLOCATE(timestep_samples(ix)%frc_save(nPart,3)) ! N*3 matrix
    ALLOCATE(timestep_samples(ix)%EAM_save(nPart))   ! N*1 matrix
    timestep_samples(ix)%xyz_save(:,:) = 0 
    timestep_samples(ix)%vel_save(:,:) = 0       
    timestep_samples(ix)%acc_save(:,:) = 0
    timestep_samples(ix)%frc_save(:,:) = 0
    timestep_samples(ix)%EAM_save(:)   = 0
    timestep_samples(ix)%kin_save      = 0
    timestep_samples(ix)%E_save        = 0
    timestep_samples(ix)%err_save      = 0
  END DO
    
  !allocate and initialize timestep variables
  ALLOCATE(xyz(nPart,3)) ! N*3 matrix
  ALLOCATE(vel(nPart,3)) ! N*3 matrix
  ALLOCATE(acc(nPart,3)) ! N*3 matrix
  ALLOCATE(frc(nPart,3)) ! N*3 matrix
  ALLOCATE(EAM(nPart))   ! N*1 matrix
  xyz(:,:) = 0 
  vel(:,:) = 0       
  acc(:,:) = 0
  frc(:,:) = 0

!ALLOCATE(dotEAM(nPart))
!********************************************************************************
  ! read EAM potential data from Ni-Co eam.alloy file
  ! unites: Angestrom, ev
  CALL do_get_EAMPotData(EAMdata)
  ! calculate derivitives of EAM nested functions
  CALL do_calc_drEAMData(EAMdata, EAMdiff) 
 
!********************************************************************************
!time step loop
!********************************************************************************
  sampStp = 0
  sampK = 0
  DO tstp= 0, 2!nTstps
    IF (tstp==0) THEN
      CALL do_rand_xyz(xyz, volume)
      CALL do_rand_vel(vel, temperature, AtMas)          
      CALL do_fix_centerOfmass(vel, AtMas)	   
    ELSE
      CALL do_velverlet(xyz, vel, acc, frc, AtMas, dt)
    END IF
    ! Thermostat - Scale velocity with designated T
    CALL do_scale_vel(vel, temperature, AtMas)
    ! Perodic boundary - Recursive
	CALL do_FixXYZ(xyz(1:nPart,1:3), volume)
!REAL(dp) :: A(4,3)  
!A(1,:)=(/1.2,-2.3,4.2/)
!A(2,:)=(/2.2,1.2,-4.3/)
!A(3,:)=(/1.2,1.2,10.0/)
!A(4,:)=(/-1.5,-1.8,-4.1/)
!PRINT*, '*******************original', tstp
!PRINT*, A(1,:)
!PRINT*, A(2,:)
!PRINT*, A(3,:)
!PRINT*, A(4,:)

    ! Calculate force
    CALL do_calcFrc(frc, EAM, xyz, AtMas, AtNums, atNumsVec, n_atoms, EAMdata, EAMdiff)
    ! Calculate total kinetic energy

    kin = 0 
    DO j= 1, 3
      Kin = kin + 0.5 * sum(AtMas(:) * vel(:,j)**2)	  
    END DO  
!Print*, '(2)', kin		

!    kin =0
!    Kin =  sum(0.5 * AtMas(1:nPart) * sum(vel(1:nPart,1:3)**2)	  	) / npart
!Print*, '(3)', kin	

    ! sample initial values
    E = sum(EAM(:)) + Kin
    IF (tstp==0) THEN
      E0 = E0
    END IF
    
    IF ( tstp==sampStp ) THEN
	  timestep_samples(sampK)%xyz_save(:,:) = xyz(:,:)    
      timestep_samples(sampK)%vel_save(:,:) = vel(:,:)
      timestep_samples(sampK)%acc_save(:,:) = acc(:,:)
      timestep_samples(sampK)%frc_save(:,:) = frc(:,:)
      timestep_samples(sampK)%EAM_save(:) = EAM(:)
      timestep_samples(sampK)%kin_save = kin
      timestep_samples(sampK)%E_save = E
      timestep_samples(sampK)%err_save = (E-E0)/E0
      IF (sampStp+sampOffset <= nTstps) THEN
        sampStp = sampStp + sampOffset  
      ELSE
        sampStp = nTstps  
      END IF    
	  sampK = sampK + 1
    END IF
!Print*, tstp
  END DO
  !ToDo  : print to file
  
!********************************************************************************
END PROGRAM molecular_dynamics
!********************************************************************************
!write(*,*) 'number of atoms: ', n_atoms(:)
!write(*,*) 'atomic numbers: ', AtNums(:)
!write(*,*) 'AtMas: ',AtMas(1:size(AtMas))
!write(*,*) 't: ', t, 'dt: ', dt, 'T: ', temperature, 'V:', volume, &
!           'sampling offset: ', sampOffset,  'number of timesteps: ', n_timestps
!write(*,*) 'number of particles: ', nPart, 'number of species:', nSpcis           
