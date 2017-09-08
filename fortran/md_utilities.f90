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
  IMPLICIT NONE

  REAL(dp), ALLOCATABLE :: xyz(:,:) &! coordinate 
                         , vel(:,:) &! velocitie
                         , acc(:,:) &! acceleration
                         , frc(:,:) &! forces
                         , EAM(:)    ! EAM potentials      

  INTEGER :: nPart   &! total number of atoms in the system
           , nSamp   &! number of samples
           , ix        ! counters

  TYPE(timestep_sample), DIMENSION(:), ALLOCATABLE :: samp

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
     IMPLICIT NONE
     REAL(dp)  :: v(:,:) &
                , r(3)     &
                , vstd     &
                , m(:)&
                , T
     INTEGER :: jx  &
              , np
     
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
     IMPLICIT NONE
     
     INTEGER :: ix  &
              , np
     
     REAL(dp) :: v(:,:)                         &
               , vcm_tmp(SIZE(v,1), SIZE(v,2))  &
               , vcm = 0                        &
               , m(:) 

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
     IMPLICIT NONE
     
     INTEGER :: k  &
              , ix  &
              , np
     
     REAL(dp) :: v(:,:)   &
               , T          &
               , Tscale     &
               , Vscale   &
               , Ttemp      &
               , m(:)
       
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
     IMPLICIT NONE
     
     INTEGER :: i    & 
              , nRho &
              , nR  
     
     REAL(dp) :: xyz(:,:)    &
               , v(:,:)      &
               , a(:,:)      &
               , f(:,:)      &
               , m(:)        &
               , dt                                   
     
     REAL(dp), ALLOCATABLE :: mdiv(:)

     ALLOCATE(mdiv(SIZE(m)))
     mdiv(:) = 1 / m(:)
     DO i = 1 , 3
       xyz(:,i) = xyz(:,i) + v(:,i) * dt + 0.5 * a(:,i) * dt * dt
       v(:,i) = v(:,i) + 0.5 * ( (f(:,i) * mdiv(:)) + a(:,i) ) * dt
       a(:,i) = f(:,i) * mdiv(:)
     END DO
  
END SUBROUTINE do_velVerlet 
  
!########################################################## 
!########################################################## 
SUBROUTINE do_calcFrc(F, U, xyz, m, AtN, atNVec, nAt, EAMdata)
!##########################################################    
  IMPLICIT NONE
  REAL(dp), ALLOCATABLE :: F(:,:)  &
                         , U(:)      &
                         , xyz(:,:)  &
                         , m(:)      &
                         , Ud(:)     &    
                         , xyzd(:,:)
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
  IMPLICIT NONE
  REAL(dp) :: xyz(3) &! positions
            , V        &! volume
            , L         ! lenght
  INTEGER :: ix 
    
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

!########################################################## 
!########################################################## 
SUBROUTINE do_calcKin(K, m ,v)
!########################################################## 
  IMPLICIT NONE
 
  REAL(dp) :: k     &! kinetic energy
            , m(:)  &! mass
            , v(:,:) ! velocity
  INTEGER :: jx

  k = 0D0
  DO jx= 1, 3
    K = k + 0.5 * sum(m(:) * v(:,jx)**2)
  END DO  
  
END SUBROUTINE do_calcKin  

!Print*, '(2)', kin		
!    kin =0
!    Kin =  sum(0.5 * AtMas(1:nPart) * sum(vel(1:nPart,1:3)**2)	  	) / npart
!Print*, '(3)', kin	

!########################################################## 
!########################################################## 
SUBROUTINE do_calcT(T, v, m)
!##########################################################    
    IMPLICIT NONE
    REAL(dp):: v(:,:)    &
             , v_sqr(SIZE(v,1)) &
             , T        &
             , m(:)
    INTEGER :: dof &
             , np  !

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
  IMPLICIT NONE
  
  REAL(dp), ALLOCATABLE :: xyz(:,:)  &
                         , m(:)      
  
  TYPE(EAM_data) :: EAMData
                  
  REAL(dp) :: U, Rho, FxRho, Phi, r
  
  INTEGER :: ix, jx, kx, np &
           , ridx, rhoidx   &    
           , nAt(:)         &!number of atoms of each epecies
           , AtN(:)         &!atomic numbers A1,A2
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
      IF ( (0 < ridx) .AND. (ridx<=floor(EAMData%cutoff/EAMData%dr)) ) THEN
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
  IF (atNVec(ix)==atN(1)) THEN        ! 27 
    FxRho = Rho * EAMData%F2(rhoidx)    ! F2  
  ELSE                                ! 28
    FxRho = Rho * EAMData%F1(rhoidx)    ! F1
  END IF

  U =  FxRho + 0.5 * Phi
  
END SUBROUTINE do_calcEAMPot

!##########################################################
!##########################################################
SUBROUTINE do_genRand(r,m,j)
!##########################################################
  IMPLICIT NONE
  
  INTEGER :: i     &
           , j     &
           , m     &
           , n     &
           , clock
  
  INTEGER, ALLOCATABLE :: seed(:)
  
  REAL(dp) :: r(m)
    
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
  ! Based on Box–Muller transform
  ! More ifor in : 
  ! http://www.alanzucconi.com/2015/09/16/how-to-sample-from-a-gaussian-distribution/
  ! https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
  IMPLICIT NONE
  INTEGER  :: i               & ! Loop variabel
            , n                 ! Lenght of array = number of random samples / from input
  REAL(dp) :: array(n)        &
            , rand_normal2(n) & ! Function Return :: dp random numbers
            , pi              & ! Pi number
            , temp            & ! just a holder
            , mean = 0.0      & ! from http://edoras.sdsu.edu/doc/matlab/techdoc/ref/randn.html //equlized with randn function in matlab 
            , sd = 1.0          ! from http://edoras.sdsu.edu/doc/matlab/techdoc/ref/randn.html //equlized with randn function in matlab 
    
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


END MODULE util


!!!~!##########################################################   
!!!~!##########################################################   
!!!~SUBROUTINE do_fdmDiff(dfx, f, dx)
!!!~!##########################################################    
!!!~  IMPLICIT NONE
!!!~  REAL(dp) :: f(:) &
!!!~            , dx   &
!!!~            , h
!!!~  
!!!~  REAL(dp), ALLOCATABLE :: dfx(:)
!!!~  INTEGER :: L , ix
!!!~  
!!!~  L = SIZE(f,1)
!!!~  ALLOCATE(dfx(1:L))
!!!~  dfx(:) = 0.0
!!!~  h = 1 / (2*dx)  
!!!~  DO ix = 2, L-1
!!!~    dfx(ix) = (f(ix+1)-f(ix-1)) * h
!!!~  END DO
!!!~  dfx(1) = dfx(2) 
!!!~  dfx(L) = dfx(L-1)
!!!~
!!!~END SUBROUTINE do_fdmDiff
!!!~!##########################################################    
