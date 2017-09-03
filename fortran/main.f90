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
            , drho, dr, cutoff       ! EAM potentail parameters

  REAL(dp), ALLOCATABLE :: atomic_masses(:)  &! atomic mass of species 
                         , xyz(:,:)          &! coordinate 
                         , vel(:,:)          &! velocitie
                         , acc(:,:)          &! acceleration
                         , frc(:,:)          &! forces
                         , F1   (:)          &   
                         , rho1 (:)          &
                         , F2   (:)          &
                         , rho2 (:)          &
                         , phi11(:)          &
                         , phi21(:)          &
                         , phi22(:)          &
                         , EAM(:)            &
                         , dFdr1(:)          &
                         , drhodr1(:)        &
                         , dFdr2(:)          &
                         , drhodr2(:)        &
                         , dphidr11(:)       &
                         , dphidr21(:)       &
                         , dphidr22(:)       &
                         , dEAMdr(:)

  INTEGER :: n_particles     &! total number of atoms in the system
           , n_timestps      &! number of timesteps
           , n_species       &! number of species in the system
           , n_samples       &! counter on number of samples
           , samp_offset     &! number of timesteps between two saves
           , i ,j , k        &!
           , nrho, nr         ! AEM potentials 
             
  INTEGER, ALLOCATABLE :: atomic_numbers(:)  &! atomic numbers of each species A1, A2, ... 
                        , n_atoms(:)         ! numbers of atoms of each species n1, n2 ,..
    
  TYPE(timestep_saver), DIMENSION(:), ALLOCATABLE :: timestep_samples
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! MAIN
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  WRITE(*,*) 'This is my simple MD code, implemneted with fortran 90'
! Read from input: 
! all vectores atomic number ordered, ascending   
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
  
  ! read EAM potential data from Ni-Co eam.alloy file
  CALL do_get_EAMPotData(nrho, drho, nr, dr, cutoff, &
                         F1, rho1, F2, rho2, phi11, phi21, phi22)

  CALL do_calc_drEAMParam(dr, drho, F1, rho1, F2, rho2, phi11, phi21, phi22, &
                dFdr1, drhodr1, dFdr2, drhodr2, dphidr11, dphidr21, dphidr22) 
 
  DO i= 0 , 1!n_timestps    
    IF (i == 0) then          
    ! randomize xyz of atoms in cell 0 
!print *, tmstp(i)%xyz    
      CALL do_rand_xyz(xyz, volume, n_particles)
!write(*,*) xyz(:,1)
!write(*,*) xyz(:,2)
!write(*,*) xyz(:,3)
    ! randomize velocities of atoms in cell zero
      CALL do_rand_vel(vel, temperature, n_particles, atomic_masses)          
!write(*,*) vel(:,1)	
!write(*,*) vel(:,2)	
!write(*,*) vel(:,3)	
!CALL do_calcT(temperature_temp, vel, n_particles, atomic_masses )
!write(*,*) 'calculated T= ', temperature_temp, 'designated T= ', temperature

    ! fix center of mass
    CALL do_fix_centerOfmass(vel, n_particles, atomic_masses)
    ! Scale velocities with initial temperature
  
    CALL do_scale_vel(vel, temperature, n_particles, atomic_masses)                               
    ! ToDo: question? should values travers to i+1?       
    ! save sample time step
  
    ELSE
      CALL do_calcEAM(EAM, xyz, atomic_masses, n_atoms,   &
                   nrho, drho, nr, dr, cutoff,            &
                   F1, rho1, F2, rho2, phi11, phi21, phi22)
      ! velocityverlet
      CALL do_velverlet( n_particles, xyz, vel, acc, frc, &
                   EAM, atomic_masses, dt )  
      ! update forces
      !frc(1,:) =1+1.5* tmstp(i)%frc(1,:)
      ! ToDo: save sample timestep
    END IF
  END DO
 
  !ToDo  : print to file
END PROGRAM molecular_dynamics



!  call do_alloc_dtyp(timestep_samples, n_samples+2)
  ! loops overtime steps

  !CALL do_loop_tstps(tmstp)
