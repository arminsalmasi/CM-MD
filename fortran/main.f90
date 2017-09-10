!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Program MD
! a fortran code for a simple Molocular Dynamics simulation of Ni, Co alloy using
! EAM potentials from NIST EAM repository: 
! https://www.ctcms.nist.gov/potentials/Co-Ni.html
! Potential data is tabulated by Y. Mishin (George Mason Univ.) on 17 Sept. 2013 
! and was posted on 17 Jan. 2014. This version is compatible with LAMMPS and The reference 
! information was updated on 26 Aug. 2015 
! Format: EAM/alloy setfl 
! https://www.ctcms.nist.gov/potentials/Download/Ni-Al-Co-YM13/Mishin-Ni-Co-2013.eam.alloy
! This code is implemented by Armin Salmasi at KTH, Stockholm,Sweden in 2017
! Feel free to use any part of this code.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM md
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  USE io    ! Input output
  USE glbl  ! Global variables and data types
  USE util  ! md functions and subroutins
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  IMPLICIT NONE
  REAL(dp) :: t, dt, Temp, vol &! simulation, time timestep [sec], vol [Ang^3]
            , smpInt           &! time between saving data (purpes = speed up) [sec]
            , E, E0            &! Total energy ::  time dependent , initial [ev]
            , Kin               ! Kinetic energy [ev]

  REAL(dp), ALLOCATABLE :: AtMas(:) &! atomic mass of species [gr]
                         , xyz(:,:) &! coordinate [ang]
                         , vel(:,:) &! velocitie [ang/sec]
                         , acc(:,:) &! acceleration [ang/sec^2]
                         , frc(:,:) &! forces []
                         , EAM(:)    ! EAM potentials [ev]      

  INTEGER :: nPart   &! total number of atoms in the system
           , nTstps  &! number of timesteps
           , nSpcis  &! number of species in the system
           , nSamp   &! number of samples
           , sampOfs &! number of timesteps between two samples except for the last one
           , samptp  &! tstp of sampling
           , sampK   &! sample counter  
           , tstp     ! timestep counter 

  INTEGER, ALLOCATABLE :: AtNums(:)    &! atomic numbers of each species A1, A2, ... 
                        , atNumsVec(:) &!  A1,A1,A1,....; A2,A2,A2... 
                        , nAtms(:)      ! numbers of atoms of each species n1, n2 ,..

  TYPE(timestep_sample), ALLOCATABLE :: samp(:) ! keeps samples from timesteps
 
  TYPE(EAM_data) :: EAMdata ! keeps parameters for calculation of EMA potentials        
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! MAIN
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Read from input: all vectores atomic number ordered, ascending   
! Lenght unit: angestrom
  CALL get_input(t, dt, Temp, vol, nTstps, smpInt, &
                 nPart, nSpcis, nAtms, AtNums, atNumsVec, AtMas)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Allocate and initialize all to Zero
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  sampOfs = Floor(smpInt / dt * 1.0)
  nSamp = 1 + FLOOR(nTstps / sampOfs * 1.0)
  CALL do_Alloc(xyz, vel, frc, acc, EAM, samp, nSamp, nPart)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! read EAM potential data from Ni-Co eam.alloy file
! unites: Angestrom, ev
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  CALL do_get_EAMPotData(EAMdata)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!time step loop
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  samptp = 0
  sampK = 0

  DO tstp= 0, nTstps
    IF (tstp==0) THEN
      CALL do_rand_xyz(xyz, vol)
      CALL do_rand_vel(vel, Temp, AtMas)          
      CALL do_fix_centerOfmass(vel, AtMas)
    ELSE
      CALL do_velverlet(xyz, vel, acc, frc, AtMas, dt)
    END IF
! Thermostat - Scale velocity with designated T
    CALL do_scaleVel(vel, Temp, AtMas)

! Perodic boundary - Recursive
    CALL do_FixXYZ(xyz(1:nPart,1:3), vol)

! Calculate force
    CALL do_calcFrc(frc, EAM, xyz, AtMas, AtNums, atNumsVec, nAtms, EAMdata)

! Calculate total kinetic energy
    CALL do_calcKin(Kin, AtMas, vel)

! calculate total energy
    E = SUM(EAM(:)) + Kin

! save E0 to calculate Error 
    IF (tstp==0) THEN
      E0 = E
    END IF

! save samples    
    IF ( tstp==samptp ) THEN
      samp(sampK)%xyz(:,:) = xyz(:,:)    
      samp(sampK)%vel(:,:) = vel(:,:)
      samp(sampK)%acc(:,:) = acc(:,:)
      samp(sampK)%frc(:,:) = frc(:,:)
      samp(sampK)%EAM(:)   = EAM(:)
      samp(sampK)%kin      = kin
      samp(sampK)%E        = E
      samp(sampK)%error    = (E-E0)/E0
      samp(sampK)%time     = tstp * dt
   
      IF (samptp+sampOfs <= nTstps) THEN
        samptp = samptp + sampOfs  
      ELSE
        samptp = nTstps  
      END IF    
     
      sampK = sampK + 1
      WRITE(*,*), 'timestep ', tstp, ' from ', nTstps, 'dt (sec)= ', dt &
                , 't (sec)= ', nTstps * dt
    END IF
  END DO
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! write sampes to md files
CALL put_SampsToFiles(samp, nPart, nSpcis, atNums, nAtms, Temp, vol, t, dt, nTstps &
                            , smpInt, sampK)

WRITE(*,*) 'Run post.m to see some results'
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END PROGRAM md
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

