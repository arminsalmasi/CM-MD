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

  REAL(dp), ALLOCATABLE :: atomic_masses(:)  &! atomic mass of species 
                         , xyz(:,:)          &! coordinate 
                         , vel(:,:)          &! velocitie
                         , acc(:,:)          &! acceleration
                         , frc(:,:)          &! forces
, EAM(:)         

  INTEGER :: n_particles     &! total number of atoms in the system
           , n_timesteps     &! number of timesteps
           , n_species       &! number of species in the system
           , n_samples       &! number of samples
           , samp_offset     &! number of timesteps between two samples except for the last one
           , sampleStep      &! tstp of sampling
           , i ,j , k         ! counters
  INTEGER, ALLOCATABLE :: atomic_numbers(:)  &! atomic numbers of each species A1, A2, ... 
                        , n_atoms(:)         ! numbers of atoms of each species n1, n2 ,..
  TYPE(timestep_saver), DIMENSION(:), ALLOCATABLE :: timestep_samples
  TYPE(EAM_data) :: EAMdata, EAMdiff         

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! MAIN
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  WRITE(*,*) 'This is my simple MD code, implemneted with fortran 90'
  ! Read from input: all vectores atomic number ordered, ascending   
  ! Lenght unit: angestrom
  CALL get_input(t, dt, temperature, volume, n_timesteps, samp_interval, &
                 n_particles, n_species, n_atoms, atomic_numbers, atomic_masses)
!write(*,*) 'number of atoms: ', n_atoms(:)
!write(*,*) 'atomic numbers: ', atomic_numbers(:)
!write(*,*) 'atomic_masses: ',atomic_masses(1:size(atomic_masses))
!write(*,*) 't: ', t, 'dt: ', dt, 'T: ', temperature, 'V:', volume, &
!           'sampling offset: ', samp_offset,  'number of timesteps: ', n_timestps
!write(*,*) 'number of particles: ', n_particles, 'number of species:', n_species           

  ! allocate and initialize saver
  samp_offset = Floor(samp_interval / dt * 1.0)
!  samp_offset = FLOOR(samp_offset)  
  n_samples = 1 + FLOOR(n_timesteps / samp_offset * 1.0)
  ALLOCATE(timestep_samples(0:n_samples))                           
  DO i = 0, n_samples
    ALLOCATE(timestep_samples(i)%xyz_save(n_particles,3)) ! N*3 matrix
    ALLOCATE(timestep_samples(i)%vel_save(n_particles,3)) ! N*3 matrix
    ALLOCATE(timestep_samples(i)%acc_save(n_particles,3)) ! N*3 matrix
    ALLOCATE(timestep_samples(i)%frc_save(n_particles,3)) ! N*3 matrix
    ALLOCATE(timestep_samples(i)%EAM_save(n_particles))   ! N*1 matrix
    timestep_samples(i)%xyz_save(:,:) = 0 
    timestep_samples(i)%vel_save(:,:) = 0       
    timestep_samples(i)%acc_save(:,:) = 0
    timestep_samples(i)%frc_save(:,:) = 0
    timestep_samples(i)%EAM_save(:)   = 0
    timestep_samples(i)%kin_save      = 0
    timestep_samples(i)%E_save        = 0
    timestep_samples(i)%err_save      = 0
  END DO
    
  !allocate and initialize timestep variables
  ALLOCATE(xyz(n_particles,3)) ! N*3 matrix
  ALLOCATE(vel(n_particles,3)) ! N*3 matrix
  ALLOCATE(acc(n_particles,3)) ! N*3 matrix
  ALLOCATE(frc(n_particles,3)) ! N*3 matrix
  ALLOCATE(EAM(n_particles))   ! N*1 matrix
  xyz(:,:) = 0 
  vel(:,:) = 0       
  acc(:,:) = 0
  frc(:,:) = 0
  
!ALLOCATE(dotEAM(n_particles))
!********************************************************************************
  ! read EAM potential data from Ni-Co eam.alloy file
  ! unites: Angestrom, ev
  CALL do_get_EAMPotData(EAMdata)
  ! calculate derivitives of EAM nested functions
  CALL do_calc_drEAMData(EAMdata, EAMdiff) 
 
!********************************************************************************
!time step loop
!********************************************************************************
  sampleStep = 0
  DO i= 0, n_timesteps
Print*, i
    IF (i==0) THEN
      CALL do_rand_xyz(xyz, volume)
      CALL do_rand_vel(vel, temperature, atomic_masses)          
      CALL do_fix_centerOfmass(vel, atomic_masses)
    ELSE
      CALL do_velverlet(xyz, vel, acc, frc, atomic_masses, dt)
    END IF
    ! Thermostat - Scale velocity with designated T
    CALL do_scale_vel(vel, temperature, atomic_masses)
!A(1,:)=(/1.2,-2.3,4.2/)
!A(2,:)=(/2.2,1.2,-4.3/)
!A(3,:)=(/1.2,1.2,10.0/)
!A(4,:)=(/-1.5,-1.8,-4.1/)
    ! Perodic boundary - Recursive
    DO j = 1 , n_particles
      CALL do_FixXYZ(xyz(j,:), volume)
    END DO
!PRINT*, A(1,:)
!PRINT*, A(2,:)
!PRINT*, A(3,:)
!PRINT*, A(4,:)
    ! Calculate force
    CALL do_calcEAMfrc(frc, EAM, xyz, atomic_masses, n_atoms, EAMdata, EAMdiff)
    ! Calculate total kinetic energy
    Kin = 0
    DO j= 1, n_particles
      Kin = 0.5 * atomic_masses(j) * sum(vel(j,:)**2)
    END DO   
    ! sample initial values
    E = sum(EAM(:)) + Kin
    IF (i==0) THEN
      E0 = E0
    END IF
    
    IF ( i==sampleStep ) THEN
      timestep_samples(i)%xyz_save(:,:) = xyz(:,:)    
      timestep_samples(i)%vel_save(:,:) = vel(:,:)
      timestep_samples(i)%acc_save(:,:) = acc(:,:)
      timestep_samples(i)%frc_save(:,:) = frc(:,:)
      timestep_samples(i)%EAM_save(:) = EAM(:)
      timestep_samples(i)%kin_save = kin
      timestep_samples(i)%E_save = E
      timestep_samples(i)%err_save = (E-E0)/E0
      IF (sampleStep+samp_offset <= n_timesteps) THEN
        sampleStep = sampleStep + samp_offset  
      ELSE
        sampleStep = n_timesteps  
      END IF    
    END IF
Print*, i

  END DO
  !ToDo  : print to file
  
!********************************************************************************
END PROGRAM molecular_dynamics
!********************************************************************************
