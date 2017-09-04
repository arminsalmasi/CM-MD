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

  REAL(dp), ALLOCATABLE :: atomic_masses(:)  &! atomic mass of species 
                         , xyz(:,:)          &! coordinate 
                         , vel(:,:)          &! velocitie
                         , acc(:,:)          &! acceleration
                         , frc(:,:)          &! forces


  INTEGER :: n_particles     &! total number of atoms in the system
           , n_timestps      &! number of timesteps
           , n_species       &! number of species in the system
           , n_samples       &! counter on number of samples
           , samp_offset     &! number of timesteps between two saves
           , i ,j , k        &!
             
  INTEGER, ALLOCATABLE :: atomic_numbers(:)  &! atomic numbers of each species A1, A2, ... 
                        , n_atoms(:)         ! numbers of atoms of each species n1, n2 ,..
    
  TYPE(timestep_saver), DIMENSION(:), ALLOCATABLE :: timestep_samples
  
  TYPE(EAM_data) :: EAMdata, EAMdiff         

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! MAIN
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  WRITE(*,*) 'This is my simple MD code, implemneted with fortran 90'

!********************************************************************************
! Initialization   
!********************************************************************************
  ! Read from input: all vectores atomic number ordered, ascending   
  ! Lenght unit: angestrom
  CALL get_input(t, dt, temperature, volume, n_timestps, &
                 samp_offset, samp_interval, n_samples, &
                 n_particles, n_species, n_atoms, atomic_numbers, atomic_masses)
!write(*,*) 'number of atoms: ', n_atoms(:)
!write(*,*) 'atomic numbers: ', atomic_numbers(:)
!write(*,*) 'atomic_masses: ',atomic_masses(1:size(atomic_masses))
!write(*,*) 't: ', t, 'dt: ', dt, 'T: ', temperature, 'V:', volume, &
!           'sampling offset: ', samp_offset,  'number of timesteps: ', n_timestps
!write(*,*) 'number of particles: ', n_particles, 'number of species:', n_species           

  ! allocate and initialize saver
  ALLOCATE(timestep_samples(0:n_samples))                           
    DO i = 0, n_samples
      ALLOCATE(timestep_samples(i)%pos_save(n_particles,3)) ! N*3 matrix
      ALLOCATE(timestep_samples(i)%vel_save(n_particles,3)) ! N*3 matrix
      ALLOCATE(timestep_samples(i)%acc_save(n_particles,3)) ! N*3 matrix
      ALLOCATE(timestep_samples(i)%frc_save(n_particles,3)) ! N*3 matrix
      timestep_samples(i)%pos_save(:,:) = 0 
      timestep_samples(i)%vel_save(:,:) = 0       
      timestep_samples(i)%acc_save(:,:) = 0
      timestep_samples(i)%frc_save(:,:) = 0
    END DO
    
  !allocate and initialize timestep variables
  ALLOCATE(xyz(n_particles,3)) ! N*3 matrix
  ALLOCATE(vel(n_particles,3)) ! N*3 matrix
  ALLOCATE(acc(n_particles,3)) ! N*3 matrix
  ALLOCATE(frc(n_particles,3)) ! N*3 matrix
  xyz(:,:) = 0 
  vel(:,:) = 0       
  acc(:,:) = 0
  frc(:,:) = 0
!********************************************************************************
  ! read EAM potential data from Ni-Co eam.alloy file
  ! unites: Angestrom, ev
  CALL do_get_EAMPotData(EAMdata)
  ! calculate derivitives of EAM nested functions
  CALL do_calc_drEAMParam(EAMdata, EAMdiffi) 
!********************************************************************************
  ! randomize xyz of atoms in cell 0 
  CALL do_rand_xyz(xyz, volume)
  
  ! randomize velocities of atoms in cell zero
  CALL do_rand_vel(vel, temperature, atomic_masses)          

  ! fix center of mass
  CALL do_fix_centerOfmass(vel, atomic_masses)
  ! Scale velocities with initial temperature
  
  CALL do_scale_vel(vel, temperature, atomic_masses)                               
 
  ! ToDo: save sampel 0       
 
!********************************************************************************
!time step loop
!********************************************************************************
  DO i= 1 , 2!n_timestps      
      ! velocityverlet
      CALL do_velverlet( xyz, vel, acc, frc, atomic_masses, dt , EAMdata,EAMdiff)
      ! ToDo: save sample timestep
  END DO
 
  !ToDo  : print to file
  
!********************************************************************************
END PROGRAM molecular_dynamics
!********************************************************************************
