!##########################################################  
! Module util (md_utilities) part of md program
! Contains:
!  do_alloc : allocate datatypes and allocatbale parameters of main
!  do_rand_xyz : randomize initial positions of atoms in the box
!  do_rand_vel : randomize initial velocity of atoms using temperature
!  do_fix_CenterOfmas : fixes center of mass of the system
!  do_scaleVel : sacle velocity with temperature : thermostat
!  do_velVerlet : velocity verlet algorithm
!  do_calcFrc : calcualte forces from EAM potentials 
!  do_FixXYZ : periodic boundry condition
!  do_calcKin : calcualte kinetic energy of the system
!  do_calcT : calcualte temperature of the system
!  do_calcEAMpot : calcualte single point EAM potential of Ni-Co system
!  do_genrand : generate a vector of rendom numbers
!  randNormal2 : pic random numbers from a normal distrition
!Armin salmasi 2017-09-09
!##########################################################  
MODULE util
!##########################################################  
  USE glbl
!##########################################################  
  REAL, PRIVATE    :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0
  PRIVATE          :: integral
!##########################################################  
CONTAINS   
!##########################################################  
!##########################################################  
SUBROUTINE do_Alloc(xyz, vel, frc, acc, EAM, samp, nSamp, nPart)
!##########################################################  
! Initial allocation of arrays and data types:
! xyz, velocity, force, acceleration, 
! EAM Potentials array
! and sampel holder
! by nsamples and nParticles
!##########################################################  
  IMPLICIT NONE

  REAL(dp), ALLOCATABLE :: xyz(:,:) &! coordinate 
                         , vel(:,:) &! velocitie
                         , acc(:,:) &! acceleration
                         , frc(:,:) &! forces
                         , EAM(:)    ! EAM potentials      

  INTEGER :: nPart   &! total number of atoms in the system
           , nSamp   &! number of samples
           , ix        ! counters

  TYPE(timestep_sample), ALLOCATABLE :: samp(:)

  ALLOCATE(samp(0:nSamp))                           
  DO ix = 0, nSamp
    ALLOCATE(samp(ix)%xyz(nPart,3)) ! N*3 matrix
    ALLOCATE(samp(ix)%vel(nPart,3)) ! N*3 matrix
    ALLOCATE(samp(ix)%acc(nPart,3)) ! N*3 matrix
    ALLOCATE(samp(ix)%frc(nPart,3)) ! N*3 matrix
    ALLOCATE(samp(ix)%EAM(nPart))   ! N*1 matrix
    samp(ix)%xyz(:,:) = 0D0 
    samp(ix)%vel(:,:) = 0D0        
    samp(ix)%acc(:,:) = 0D0 
    samp(ix)%frc(:,:) = 0D0 
    samp(ix)%EAM(:)   = 0D0 
    samp(ix)%kin      = 0D0 
    samp(ix)%E        = 0D0 
    samp(ix)%error    = 0D0 
  END DO
    
! Allocate and initialize timestep variables
  ALLOCATE(xyz(nPart,3)) ! N*3 matrix
  ALLOCATE(vel(nPart,3)) ! N*3 matrix
  ALLOCATE(acc(nPart,3)) ! N*3 matrix
  ALLOCATE(frc(nPart,3)) ! N*3 matrix
  ALLOCATE(EAM(nPart))   ! N*1 matrix
  xyz(:,:) = 0D0  
  vel(:,:) = 0D0        
  acc(:,:) = 0D0 
  frc(:,:) = 0D0  

END SUBROUTINE do_Alloc

!##########################################################
!##########################################################
SUBROUTINE do_rand_xyz(xyz, V)
!##########################################################   
! randomize xyz positions of atoms inside the given volume 
!##########################################################   
     INTEGER :: ix &! counter
              , np  ! number of particles
     
     REAL(dp)  :: r(3)     & ! random number 
                , xyz(:,:) & ! position n*3
                , V          ! volume

     np = SIZE(xyz,1)
     DO ix = 1, np
       CALL do_genRand(r, 3, ix)
       xyz(ix,:) = r(:) * ( V **(1.0/3.0))
     END DO
   
END SUBROUTINE do_rand_xyz
   
!##########################################################
!##########################################################
SUBROUTINE do_rand_vel(v, T, m)
!##########################################################     
! randomize initial velocity of atoms with temperature
!##########################################################     
     IMPLICIT NONE

     REAL(dp)  :: v(:,:) &! velocities
                , r(3)   &! random numbers
                , vstd   &! standard deviation of velocity
                , m(:)   &! atomic masses in gram  
                , T       ! temperature

     INTEGER :: jx  &! counter
              , np   ! number of particles
     
     np = SIZE(v,1)
     DO jx = 1 , np       
       r(:) = rand_normal2(3) ! open source from roseta code
       vstd = sqrt( T * kb / m(jx) )
       v(jx,:) = 0 + r(:) * vstd   
     END Do
   
END SUBROUTINE do_rand_vel
   
!##########################################################
!##########################################################
SUBROUTINE do_fix_centerOfMass(v, m)
!##########################################################     
! fix center of mass of the system
!##########################################################     
  IMPLICIT NONE
  
  INTEGER :: ix &
           , np  ! number of particles
  
  REAL(dp) :: v(:,:)                &! velocity
            , vcm_tmp(SIZE(v,1), 3) &! swap velocity
            , vcm = 0               &! velocity of center of mass
            , m(:)                   ! atomic masses

  np = SIZE(v,1)
  vcm_tmp(:,:)=0
  DO ix = 1, 3
    vcm_tmp(:,ix) = v(:,ix) * m(:)
  END DO
  vcm = SUM(vcm_tmp) / SUM(m)
  v(:,1:3) = v(:,1:3) - vcm
 
END SUBROUTINE do_fix_centerOfMass
   
!##########################################################
SUBROUTINE do_scaleVel(v, T, m )
!##########################################################      
! scale velocity with temperature of the system: thermostat
! partially taken from web
!##########################################################      
  IMPLICIT NONE
  
  INTEGER :: k, ix &! counters
           , np     ! number of particles
  
  REAL(dp) :: v(:,:)  &! velocity
           , T       &! temperature
           , Tscale  &! temperature scaling parameter
           , Vscale  &! velocity scaling parameter
           , Ttemp   &! swap: temprature
           , m(:)     ! atomic masses
   
  np = SIZE(v,1)
  k = 1 
  CALL do_calcT(Ttemp, v, m)
  Tscale = ( (T - Ttemp) / T )  * 100
  DO WHILE ( (ABS(Tscale)>0.001) .AND. (k<1000) )
     Vscale = sqrt( T / Ttemp)
     DO ix = 1 , 3
       v(:,ix) = v(:,ix) * Vscale
     END DO 
     CALL do_calcT(Ttemp, v , m)
     Tscale = ( (T - Ttemp) / T )  * 100
     k= k + 1  
  END DO

END SUBROUTINE do_scaleVel                   

!##########################################################
SUBROUTINE do_velVerlet(xyz, v, a, f, m, dt)  
!##########################################################
! Velocity Verlet algorithm
!##########################################################
  IMPLICIT NONE
  
  INTEGER :: ix ! counter 
  
  REAL(dp) :: xyz(:,:)    &! xyz positions of all atoms
            , v(:,:)      &! xyz velocities of all atoms 
            , a(:,:)      &! xyz accelerations of all atoms
            , f(:,:)      &! xyz forces on all atoms
            , m(:)        &! atomic masses
            , dt           ! timestep                         
  
  REAL(dp), ALLOCATABLE :: mdiv(:)

  ALLOCATE(mdiv(SIZE(m)))

  mdiv(:) = 1 / m(:)
  DO ix = 1 , 3
    xyz(:,ix) = xyz(:,ix) + v(:,ix) * dt + 0.5 * a(:,ix) * dt * dt
    v(:,ix)   = v(:,ix) + 0.5 * ((f(:,ix) * mdiv(:)) + a(:,ix)) * dt
    a(:,ix)   = f(:,ix) * mdiv(:)
  END DO

END SUBROUTINE do_velVerlet 
  
!########################################################## 
!########################################################## 
SUBROUTINE do_calcFrc(F, U, xyz, m, AtN, atNVec, nAt, EAMdata)
!##########################################################    
! calculate forces on atoms using EAM potntials on Ni-Co
!##########################################################    
  IMPLICIT NONE

  REAL(dp), ALLOCATABLE :: F(:,:)    &! Forces
                         , U(:)      &! Potentials
                         , xyz(:,:)  &! positions
                         , m(:)      &! atomic masses
                         , Ud(:)     &! difrentiated potential    
                         , xyzd(:,:)  ! difrentiated positions
  REAL(dp) :: drdiv
  
  TYPE(EAM_data) :: EAMData 

  INTEGER :: nAt(:)         &!number of atoms of each epecies
           , AtN(:)         &!atomic numbers A1,A2
           , atNVec(:)      & ! A1,A1,... A2,A2....
           , np             &
           , ix , jx

  np = SIZE(xyz,1)
  ALLOCATE(Ud(np))
  ALLOCATE(xyzd(np,3))
  
  U(:) = 0D0
  F(:,:) = 0D0
  
  drdiv = 1D0/EAMdata%dr
  DO ix = 1, np
    CALL do_calcEAMPot( U(ix), ix, xyz, m, AtN, atNVec, nAt, EAMdata)
    DO jx = 1, 3
      xyzd(:,:) = xyz(:,:)
      xyzd(ix,jx) = xyz(ix,jx) + EAMdata%dr
      CALL do_calcEAMPot( Ud(ix), ix, xyzd, m, AtN, atNVec, nAt, EAMdata)
      F(ix,jx) = (Ud(ix) - U(ix)) * drdiv
    END DO
  END DO

END SUBROUTINE do_calcFrc

!########################################################## 
!########################################################## 
RECURSIVE SUBROUTINE do_FixXYZ(xyz, V)
!########################################################## 
! periodic boundary condition :: recursive
! toss atoms which leave the box inside 
!########################################################## 
  IMPLICIT NONE

  REAL(dp) :: xyz(3) &! positions
            , V      &! volume
            , L       ! lenght of the box side

  INTEGER :: ix ! counter
    
  L = V ** (1.0/3.0)
    DO ix=1, 3
      IF (xyz(ix)>(L/2)) THEN
        xyz(ix) = xyz(ix)-L;
        IF (xyz(ix)>(L/2)) THEN
          CALL do_FixXYZ(xyz,V)
        END IF 
      ELSE
        IF (xyz(ix)<(-L/2)) THEN
          xyz(ix) = xyz(ix)+L;
          IF (xyz(ix)<(-L/2)) THEN
            CALL do_FixXYZ(xyz,V)
          END IF        
        END IF
      END IF
    END DO

END SUBROUTINE do_FixXYZ

!########################################################## 
SUBROUTINE do_calcKin(K, m ,v)
!########################################################## 
! calculate kinetic energy of the system
!########################################################## 
  IMPLICIT NONE
 
  REAL(dp) :: k     &! kinetic energy
            , m(:)  &! mass
            , v(:,:) ! velocity

  INTEGER :: jx ! counter

  k = 0D0
  DO jx= 1, 3
    K = k + 0.5 * sum(m(:) * v(:,jx)**2)
  END DO  
  
END SUBROUTINE do_calcKin  

!########################################################## 
!########################################################## 
SUBROUTINE do_calcT(T, v, m)
!##########################################################    
! calculate temperature from velocity
!##########################################################    
  IMPLICIT NONE

  REAL(dp):: v(:,:)           &! velocities
           , v_sqr(SIZE(v,1)) &! square of sime of xyz velocities
           , T                &! temperature 
           , m(:)              ! atomic masses  

  INTEGER :: dof &! degrees of freedom
           , np   ! number of particles

  np = SIZE(v,1)
  v_sqr(:) = v(:,1)**2 + v(:,2)**2 + v(:,3)**2 
  dof = 3 * np - 3
  T = sum(m(:) * v_sqr(:))
  T = T / (kb * dof) 

END SUBROUTINE do_calcT

!##########################################################
!##########################################################
SUBROUTINE do_calcEAMPot( U, ix, xyz, m, AtN, atNVec, nAt, EAMdata)
!##########################################################
! calculate EAM potentials of Ni -Co based on NIST databse
! http://lammps.sandia.gov/doc/pair_eam.html
!https://www.ctcms.nist.gov/potentials/Co.html
!https://www.ctcms.nist.gov/potentials/Download/Ni-Al-Co-YM13/Mishin-Ni-Al-Co-2013.eam.alloy
!##########################################################
  IMPLICIT NONE
  
  REAL(dp), ALLOCATABLE :: xyz(:,:)  &! positions
                         , m(:)       ! atomic masses
  
  TYPE(EAM_data) :: EAMData ! EAM potential parameters for Ni and Co
                  
  REAL(dp) :: U      &! point potential
            , Rho    &! atomic electron density 
            , FxRho  &! embedding function * 
            , Phi    &! pair potential interaction
            , r      
  
  INTEGER :: ix, jx, kx     &! counters
           , np             &! bumber of particles 
           , ridx, rhoidx   &! distance index, electron density index    
           , nAt(:)         &! number of atoms of each epecies
           , AtN(:)         &! atomic numbers A1,A2
           , atNVec(:)       ! A1,A1,... A2,A2....

  np = SIZE(xyz,1)
  U = 0D0
  Rho = 0D0
  Phi = 0D0
  FxRho = 0D0
   
  DO jx = 1 , np
    IF (ix/=jx) THEN
      r = sqrt( sum( (xyz(ix,:) - xyz(jx,:) ) ** 2 ) )
      ridx = floor(r / (EAMData%dr))
      IF ( ( ridx <= 0) .or. (ridx>floor(EAMData%cutoff/EAMData%dr)) ) THEN
         phi = 0D0
         Rho = 0D0  
      ELSE
        ! F1, Rho1 :: 28
        ! F2, Rho2 :: 27
        ! Phi11 = 28-28
        ! phi21 = 27-28
        ! phi22 = 27-27		  
        IF (atNVec(ix)==atN(1)) THEN                    ! ix = 27
          Rho= Rho + EAMData%rho1(ridx)          ! rho1 = 28
          IF (atNVec(jx)==atN(1)) THEN                    ! jx = 27
                                                           ! phi22 = 27-27
              Phi = Phi + EAMData%phi22(ridx)/ridx/EAMData%dr        
            ELSE                                            ! jx = 28
                                                              ! phi21 = 27-28			  
              Phi = Phi + EAMData%phi21(ridx)/ridx/EAMData%dr        
            ENDIF 
          ELSE                                          ! ix = 28
            Rho= Rho + EAMData%rho2(ridx)         ! rho2 = 27 
            IF (atNVec(jx)==atN(1)) THEN                    ! jx = 27
                                                              ! phi21 = 28-27  			  
              Phi = Phi + EAMData%phi21(ridx)/ridx/EAMData%dr        
            ELSE                                            ! jx = 28
                                                              ! phi11= 28-28  
              Phi = Phi + EAMData%phi11(ridx)/ridx/EAMData%dr        
            END IF
          END IF
      END IF
    END IF
  END DO
 

 
  rhoidx = floor(rho/(EAMData%drho)) 
  IF ( ( rhoidx <= 0) .or. (rhoidx>floor(EAMData%cutoff/EAMData%dr)) ) THEN
     FxRho = 0D0
  ELSE
    IF (atNVec(ix)==atN(1)) THEN        ! 27 
      FxRho = EAMData%F2(rhoidx)          ! F2(rhoidx)  
    ELSE                                ! 28
      FxRho = EAMData%F1(rhoidx)          ! F1(rhoidx)
    END IF
  END IF
  
  U =  FxRho + 0.5 * Phi
  
END SUBROUTINE do_calcEAMPot

!##########################################################
!##########################################################
SUBROUTINE do_genRand(r,m,j)
!##########################################################
! Partly taken from web
! generates a vector r(m) with random numbrs 
!##########################################################
  IMPLICIT NONE
  
  INTEGER :: i, j, n  &! counters 
           , m        &! number of random numbers
           , clock     ! system clock 
  
  INTEGER, ALLOCATABLE :: seed(:)
  
  REAL(dp) :: r(m) ! vector of random numbers
    
  CALL RANDOM_SEED(SIZE=n)
  ALLOCATE(seed(n))
  CALL SYSTEM_CLOCK(COUNT=clock)
  clock = clock * j 
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT=seed)
  call RANDOM_NUMBER(r)
  DEALLOCATE(seed)
  
END SUBROUTINE do_genRand

!########################################################## 
!########################################################## 
FUNCTION rand_normal2(n) 
!##########################################################   
! https://rosettacode.org/wiki/Random_numbers#Fortran 
! Based on Box–Muller transform generates random numbers from a normal distribution
! More ifor in : 
! http://www.alanzucconi.com/2015/09/16/how-to-sample-from-a-gaussian-distribution/
! https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
!##########################################################   
  IMPLICIT NONE
  INTEGER  :: i               & ! Loop variabel
            , n                 ! Lenght of array = number of random samples / from input
  REAL(dp) :: array(n)        & 
            , rand_normal2(n) & ! Function Return :: dp random numbers
            , pi              & ! Pi number
            , temp            & ! just a holder
            , mean = 0.0      & ! http://edoras.sdsu.edu/doc/matlab/techdoc/ref/randn.html 
            , sd = 1.0          ! http://edoras.sdsu.edu/doc/matlab/techdoc/ref/randn.html
    
  pi = 4.0*ATAN(1.0)
  CALL do_genRand(array,n,1)
  ! Now convert to normal distribution
  DO i = 1, n-1, 2
    temp = sd * SQRT(-2.0*LOG(array(i))) * COS(2*pi*array(i+1)) + mean
    array(i+1) = sd * SQRT(-2.0*LOG(array(i))) * SIN(2*pi*array(i+1)) + mean
    array(i) = temp
  END DO
  ! Check mean and standard deviation
  mean = SUM(array)/n
  sd = SQRT(SUM((array - mean)**2)/n)     
  !WRITE(*, "(A,F8.6)") "Mean = ", mean
  !WRITE(*, "(A,F8.6)") "Standard Deviation = ", sd      
  rand_normal2(:) = array(:)

END FUNCTION rand_normal2

!##########################################################   
!##########################################################   
END MODULE util
!##########################################################   
!##########################################################   
