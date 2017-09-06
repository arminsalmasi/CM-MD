MODULE utilities
!##########################################################  
  USE global
!##########################################################  
  REAL, PRIVATE    :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0
  PRIVATE          :: integral
!##########################################################  
CONTAINS   
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
SUBROUTINE do_scale_vel(v, T, m )
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

END SUBROUTINE do_scale_vel                   

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
SUBROUTINE do_calcEAMfrc(frc, U,  xyz, m, nAt, EAMdata, EAMDiff)
!##########################################################    
  IMPLICIT NONE
  REAL(dp), ALLOCATABLE :: U(:)      &   
                         , Frc(:,:)  &
                         , xyz(:,:)  &
                         , m(:)      &
                         , rho(:)    &
                         , dotrho(:) &
                         , FxRho(:)  &
                         , dotFxdotRho(:)   &
                         , Phi(:)    &
                         , dotPhi(:)
                         
  TYPE(EAM_data) :: EAMData &
                  , EAMDiff
                  
  REAL(dp) :: r ,dr(3), ffactor(3)
  
  INTEGER :: ix, jx, k, np       &
           , ridx, rhoidx
           
  INTEGER, ALLOCATABLE :: nAt(:) !number of atoms of each epecies

  np = SIZE(xyz,1)
  ALLOCATE(rho(np))
  ALLOCATE(dotrho( np))
  ALLOCATE(FxRho(np))
  ALLOCATE(dotFxdotRho(np))
  ALLOCATE(phi(np))
  ALLOCATE(dotphi(np))
  
  U(:) = 0
  Rho(:) = 0.0
  Phi(:) = 0.0
  dotRho(:) = 0.0
  dotPhi(:) = 0.0
  FxRho(:) = 0.0
  dotFxdotRho(:) =0.0
   
  
  DO ix = 1 , np
    DO jx = 1 , np
      IF (ix/=jx) THEN
        r = sqrt( sum( (xyz(ix,:) - xyz(jx,:) ) ** 2 ) )
        ridx = floor(r / (EAMData%dr))
        IF (ridx<=floor(EAMData%cutoff/EAMData%dr) .AND. ridx>0) THEN
           IF (ix>nAt(1)) THEN ! 28 ! note: Rho1 -- correlated with F2
             rho(ix)= rho(ix) + EAMData%rho1(ridx) 
             IF (jx>nAt(1)) THEN !28-28 ! note: first 5000==28-28== phi1
               phi(ix) = phi(ix) + EAMData%phi11(ridx)
             ELSE  !28-27 ! note 2nd 5000 phi2
               phi(ix) = phi(ix) + EAMData%phi21(ridx)
             ENDIF 
           ELSE !27 ! note: rho2 -- correlated with F1
             Rho(ix)= Rho(ix) + EAMData%rho2(ridx) 
             IF (jx>nAt(1)) THEN !27-28 ! note: 2rd 5000 phi2
               phi(ix) = phi(ix) + EAMData%phi21(ridx)
             ELSE !27-27 ! note: 3rd 5000 phi 3
               phi(ix) = phi(ix) + EAMData%phi22(ridx)
             END IF
           END IF
        END IF
      END IF
    END DO
    rhoidx = floor(rho(ix)/(EAMData%drho)) 
    IF (ix>nAT(1)) THEN ! 28 -- rho1
      FxRho(ix) = rho(ix) * EAMData%F2(rhoidx)  ! F27 -- F2
!      dotFxdotRho(ix) = dotrho(ix) * EAMDiff%F2(rhoidx)  ! F27 -- F2
    ELSE ! 27 -- rho2
      FxRho(ix) = Rho(ix) * EAMData%F1(rhoidx)  ! F28 -- F1
!      dotFxdotRho(ix) = dotRho(ix) * EAMDiff%F1(rhoidx)  ! F28 -- F1
    END IF
    U(ix) =  FxRho(ix) + 0.5 * Phi(ix) 
  END DO

  frc(:,:) = 0

  DO ix = 1, (np-1)
    Do jx = (ix+1), np 
      dr(:) = xyz(ix,:) - xyz(jx,:)
      ffactor = (U(ix) - U(jx)) / dr(:)    
      frc(ix,:) = frc(ix,:) +  ffactor(:)
      frc(jx,:) = frc(jx,:) +  ffactor(:)
    END DO
  END DO

END SUBROUTINE do_calcEAMfrc

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
!########################################################## 
RECURSIVE SUBROUTINE do_FixXYZ(xyz, V)
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
SUBROUTINE do_get_EAMPotData(EAMdata)
!##########################################################    
    IMPLICIT NONE
    INTEGER :: nrho       &
             , nr         &
             , ix, jx, kx
    REAL(dp), ALLOCATABLE :: allData(:,:)    &
                           , allData_line(:) 
    TYPE(EAM_data) :: EAMdata 
    
    EAMdata%nrho = 5000 
    EAMdata%drho = 1.2999078e-03 
    EAMdata%nr = 5000 
    EAMdata%dr = 1.2999078e-03 
    EAMdata%cutoff = 6.499539
    ALLOCATE(allData(1:7000, 1:5))
    ALLOCATE(allData_line(1: SIZE(allDATA,1)* SIZE(allDATA,2)))
    ALLOCATE(EAMdata%F1(nrho))
    ALLOCATE(EAMdata%rho1(nr))
    ALLOCATE(EAMdata%F2(nrho))
    ALLOCATE(EAMdata%rho2(nr))
    ALLOCATE(EAMdata%phi11(nr))
    ALLOCATE(EAMdata%phi21(nr))
    ALLOCATE(EAMdata%phi22(nr))
    OPEN(UNIT=22,FILE="pot-Ni-Co-old.dat",FORM="FORMATTED",STATUS="OLD",ACTION="READ")
      DO ix = 1 , 7000
        READ(22,*) allData(ix,:)
      END DO
    CLOSE(UNIT=22)
    kx=1
    DO ix = 1 , 7000
      DO jx = 1 , 5
        allData_line(kx) = allData(ix,jx)
        kx=kx+1
      END DO
    END DO
    EAMdata%F1 = allData_line(1:5000)
    EAMdata%rho1 = allData_line(5001:10000)
    EAMdata%F2 = allData_line(10001:15000)
    EAMdata%rho2 = allData_line(15001:20000)
    EAMdata%phi11 = allData_line(20001:25000)
    EAMdata%phi21 = allData_line(25001:30000)
    EAMdata%phi22 = allData_line(30001:35000)

END SUBROUTINE do_get_EAMPotData

!##########################################################   
!##########################################################   
SUBROUTINE do_fdmDiff(dfx, f, dx)
!##########################################################    
  IMPLICIT NONE
  REAL(dp) :: f(:) &
            , dx   &
            , h
  
  REAL(dp), ALLOCATABLE :: dfx(:)
  INTEGER :: L , ix
  
  L = SIZE(f,1)
  ALLOCATE(dfx(1:L))
  dfx(:) = 0.0
  h = 1 / (2*dx)  
  DO ix = 2, L-1
    dfx(ix) = (f(ix+1)-f(ix-1)) * h
  END DO
  dfx(1) = dfx(2) 
  dfx(L) = dfx(L-1)

END SUBROUTINE do_fdmDiff

!##########################################################
!##########################################################   
SUBROUTINE do_calc_drEAMData(EAMdata, EAMdiff)
!##########################################################
  IMPLICIT NONE
                         
  TYPE(EAM_data) :: EAMdata, EAMdiff 
 
  CALL do_fdmDiff(EAMDiff%F1    ,EAMdata%F1    , EAMdata%drho)  
  CALL do_fdmDiff(EAMDiff%rho1  ,EAMdata%rho1 , EAMdata%dr  ) 
  CALL do_fdmDiff(EAMDiff%F2    ,EAMdata%F2   , EAMdata%drho)     
  CALL do_fdmDiff(EAMDiff%rho2  ,EAMdata%rho2 , EAMdata%dr  ) 
  CALL do_fdmDiff(EAMDiff%phi11 ,EAMdata%phi11, EAMdata%dr  ) 
  CALL do_fdmDiff(EAMDiff%phi21 ,EAMdata%phi21, EAMdata%dr  ) 
  CALL do_fdmDiff(EAMDiff%phi22 ,EAMdata%phi22, EAMdata%dr  ) 

END SUBROUTINE do_calc_drEAMData
!##########################################################
 
END MODULE utilities
