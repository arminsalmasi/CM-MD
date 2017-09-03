MODULE utilities
!##########################################################  
  USE global
  
!##########################################################  
  REAL, PRIVATE    :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0
  PRIVATE          :: integral
  
!##########################################################  
CONTAINS   

!##########################################################
SUBROUTINE do_rand_xyz(xyz, V, np)
!##########################################################   
     INTEGER :: j
     REAL(dp) :: r(3)     &! random number 
               , xyz(:,:) &! position n*3
               , V         ! volume
   
     DO j = 1, np
       call do_genRand(r,3,j)
!WRITE(*,*) 'r = ', r
       xyz(j,:) = r(:) * ( V **(1.0/3.0))
!PRINT*, 'xyz [Ang]= ', xyz(j,:) - (V **(1.0/3.0))
     END DO
   
END SUBROUTINE do_rand_xyz
   
!##########################################################
SUBROUTINE do_rand_vel(vel , tmp, np , masses)
!##########################################################     
     IMPLICIT NONE
     
     REAL(dp)  :: vel(:,:) &
                , r(3), vstd &
                , masses(:)&
                , tmp

     INTEGER :: j &
              , np
     
     DO j = 1 , np       
       r(:) = rand_normal2(3) ! open source from roseta code
       vstd = sqrt( tmp * kb / masses(j) )
       vel(j,:) = 0 + r(:) * vstd   
     END Do
   
END SUBROUTINE do_rand_vel
   
!##########################################################
SUBROUTINE do_fix_centerOfMass(vel, np, masses)
!##########################################################     
     IMPLICIT NONE

     INTEGER :: j  &
              , np
  
     REAL(dp) :: vel(:,:)                           &
               , vcm_tmp(SIZE(vel,1), SIZE(vel,2))  &
               , vcm = 0                            &
               , masses(:) !&, veltmp(SIZE(vel,1), SIZE(vel,2))

     vcm_tmp(:,:)=0
     DO j = 1, 3
       vcm_tmp(:,j) = vel(:,j) * masses(:)
     END DO
     vcm = SUM(vcm_tmp) / SUM(masses)
!PRINT*,vcm								   
     vel(:,1:3) = vel(:,1:3) - vcm
!PRINT*, vel-veltmp
    
END SUBROUTINE do_fix_centerOfMass
   
!##########################################################
SUBROUTINE do_scale_vel(vel, T, np, masses )
!##########################################################      
     IMPLICIT NONE
     
     INTEGER :: k  &
              , j  &
              , np
  
     REAL(dp) :: vel(:,:)   &
               , T          &
               , Tscale     &
               , velscale   &
               , Ttemp      &
               , masses(:)  !	&, veltemp(SIZE(vel,1), SIZE(vel,2))
     
!veltemp(:,1:3) = vel(:,1:3)
       k = 1 
       CALL do_calcT(Ttemp, vel, np , masses)
!PRINT*, Ttemp
       Tscale = ( (T - Ttemp) / T )  * 100
!PRINT *, 'Ttemp start at itteration (k) = ',k ,'is', Ttemp, 'T = ', T, 'Tscale= ', Tscale
       DO WHILE ( (ABS(Tscale)>0.001) .and. (k<1000) )
          velscale = sqrt( T / Ttemp)
!Print *, 'velscale at itteration (k) = ', k, 'is ', velscale
          DO j = 1 , 3
            vel(:,j) = vel(:,j) * velscale
          END DO 
          CALL do_calcT(Ttemp, vel, np , masses)
          Tscale = ( (T - Ttemp) / T )  * 100
!PRINT *, 'Ttemp end at itteration (k) = ',k ,'is', Ttemp, 'T = ', T, 'Tscale= ', Tscale
          k= k + 1  
       end DO
!CALL do_calcT(Ttemp, vel, np , masses)
!PRINT*, Ttemp
!veltemp = veltemp - vel
!PRINT *, veltemp (:,1)
!PRINT *, veltemp (:,2)
!PRINT *, veltemp (:,3)
END SUBROUTINE do_scale_vel                   

!##########################################################
SUBROUTINE do_velVerlet(np, r, v, a, f, U, m, dt);
!##########################################################
     
     IMPLICIT NONE              
     
     INTEGER :: i  &
                , np
  
     REAL(dp) :: r(:,:)  &
               , v(:,:)  &
               , a(:,:)  &
               , f(:,:)  &
               , m(:)    &
               , U(:)    &
               , dt      

     REAL(dp), ALLOCATABLE :: dUdr(:,:)
    
      
 
     ALLOCATE(dUdr(SIZE(r,1), 3))
     
     DO i = 1 , 3
        v(:,i) = v(:,i) + 0.5 * ( f(:,i) / m(:) + a(:,i)) * dt;
     END DO
     
     DO i=1 , 3 
       r(:,i) = r(:,i) + v(:,i) * dt + 0.5 * a(:,i) * dt * dt;
     END DO   
     
     dudr(:,:) = 0
!PRINT*, 'size(U)',size(U,1)
     !DO i= 1, 3
     !  dUdr(:,1) =  U(:, 1) - U(:,1)  !/ &
!             ( r(1:(SIZE(r,1)-1),i) - r(2:(SIZE(r,1)),i) )   
     !END DO
  
!PRINT*, U(:)

!     DO i = 1 , 3 
!        a(:,i) =  frc(:,i) / m(:);
!     END DO
   
END SUBROUTINE do_velVerlet 
  
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
SUBROUTINE do_calcT(temp, vel, np, masses)
!##########################################################    
    IMPLICIT NONE
    
    REAL(dp):: vel(:,:)    &
             , vel_sqr(np) &
             , temp        &
             , masses(:)

    INTEGER :: dof &
             , np  ! & , i
   
    vel_sqr(:) = vel(:,1)**2 + vel(:,2)**2 + vel(:,3)**2 
!do i = 1, size(vel_sqr)
!  print *, vel(i,:)
!  print *, vel_sqr(i)      
!end do
    dof = 3 * np - 3
    temp = sum(masses(:) * vel_sqr(:))
    temp = temp / (kb * dof) 
!write(*,*) temp   
END SUBROUTINE do_calcT
!########################################################## 
SUBROUTINE do_get_EAMPotData(nrho, drho, nr, dr, cutoff, F_rho_1, rho_r_1, F_rho_2, rho_r_2, phi_r11, phi_r21, phi_r22)
!##########################################################    
    IMPLICIT NONE
    
    INTEGER :: nrho       &
             , nr         &
             , ix, jx, kx
    REAL(dp), ALLOCATABLE :: allData(:,:)    &
                           , allData_line(:) &
                           , F_rho_1(:)&
                           , rho_r_1(:)&
                           , F_rho_2(:)&
                           , rho_r_2(:)&
                           , phi_r11(:)&
                           , phi_r21(:)&
                           , phi_r22(:)

    REAL(dp) :: drho  &
              , dr    & 
              , cutoff 
    
    nrho = 5000 
    drho = 1.2999078e-03 
    nr = 5000 
    dr = 1.2999078e-03 
    cutoff = 6.499539
    
    ALLOCATE(allData(1:7000, 1:5))
    ALLOCATE(allData_line(1: SIZE(allDATA,1)* SIZE(allDATA,2)))
    ALLOCATE(F_rho_1(nrho))
    ALLOCATE(rho_r_1(nr))
    ALLOCATE(F_rho_2(nrho))
    ALLOCATE(rho_r_2(nr))
    ALLOCATE(phi_r11(nr))
    ALLOCATE(phi_r21(nr))
    ALLOCATE(phi_r22(nr))
    
    OPEN(UNIT=22,FILE="pot-Ni-Co-old.dat",FORM="FORMATTED",STATUS="OLD",ACTION="READ")
      DO ix = 1 , 7000
        READ(22,*) allData(ix,:)
      END DO
!write(*,*) allData(1,:)
!write(*,*) allData(7000,:)
    CLOSE(UNIT=22)

    kx=1
    do ix = 1 , 7000
      do jx = 1 , 5
        allData_line(kx) = allData(ix,jx)
        kx=kx+1
      end do
    end do
!print*,  allData_line , size(allData_line)
    F_rho_1 = allData_line(1:5000)
    rho_r_1 = allData_line(5001:10000)
    F_rho_2 = allData_line(10001:15000)
    rho_r_2 = allData_line(15001:20000)
    phi_r11 = allData_line(20001:25000)
    phi_r21 = allData_line(25001:30000)
    phi_r22 = allData_line(30001:35000)

END SUBROUTINE do_get_EAMPotData
!########################################################## 
SUBROUTINE do_calcEAM(EAM, xyz, masses,nAt,       &
                   nrho, drho, nr, dr, cutoff,    &
                   F1, rho1, F2, rho2, phi11, phi21, phi22)
!##########################################################    
  IMPLICIT NONE
  
  REAL(dp), ALLOCATABLE :: EAM(:)    &
                         , xyz(:,:)  &
                         , masses(:) &
                         , r2d(:,:)  &
                         , r1d(:)    &
                         , F1(:)     &
                         , rho1 (:)  &
                         , F2(:)     &
                         , rho2 (:)  &
                         , phi11(:)  &
                         , phi21(:)  &
                         , phi22(:)  &
                         , rho_sum(:)&
                         , f_rho(:)  &
                         , phi_sum(:)

  REAL(dp) :: drho   &
            , dr     &
            , cutoff &
            , r
                             
  INTEGER :: ix,jx, k &
           , nrho, nr &
           , ridx, rhoidx     
  
  INTEGER, ALLOCATABLE :: nAt(:) !number of atoms of each epecies
  
  ALLOCATE(EAM(SIZE(xyz,1)))
  ALLOCATE(rho_sum( SIZE(xyz,1)))
  ALLOCATE(f_rho(SIZE(xyz,1)))
  ALLOCATE(phi_sum(SIZE(xyz,1)))

  DO ix = 1 , SIZE(xyz,1)
    rho_sum(ix) = 0.0
    phi_sum(ix) = 0.0
    DO jx = 1 , SIZE(xyz,1)
      IF (ix/=jx) THEN
        r = sqrt( sum( (xyz(ix,:) - xyz(jx,:) ) ** 2 ) )
        ridx = floor(r / (dr))
!PRINT*, 'ridx ix jx =(',ix,',', jx ,') =',  ridx
        IF (ridx<=floor(cutoff/dr) .AND. ridx>0) THEN
           IF (ix>nAt(1)) THEN ! 28 ! note: Rho1 -- correlated with F2
             rho_sum(ix)= rho_sum(ix) + rho1(ridx) 
             IF (jx>nAt(1)) THEN !28-28 ! note: first 5000==28-28== phi1
               phi_sum(ix) = phi_sum(ix) + phi11(ridx)
             ELSE  !28-27 ! note 2nd 5000 phi2
               phi_sum(ix) = phi_sum(ix) + phi21(ridx)
             ENDIF 
           ELSE !27 ! note: rho2 -- correlated with F1
             rho_sum(ix)= rho_sum(ix) + rho2(ridx) 
             IF (jx>nAt(1)) THEN !27-28 ! note: 2rd 5000 phi2
               phi_sum(ix) = phi_sum(ix) + phi21(ridx)
             ELSE !27-27 ! note: 3rd 5000 phi 3
               phi_sum(ix) = phi_sum(ix) + phi22(ridx)
             END IF
           END IF
        END IF
      END IF
    END DO
!PRINT*, 'rho_sum ix [ev]= ',ix , 'is', rho_sum(ix)
    rhoidx = floor(rho_sum(ix)/(drho)) 
!PRINT*, 'rhoidx ix = ',ix ,'=',  rhoidx
    IF (ix>nAT(1)) THEN ! 28 -- rho1
      F_rho(ix) = rho_sum(ix) * F2(rhoidx)  ! F27 -- F2
    ELSE ! 27 -- rho2
      F_rho(ix) = rho_sum(ix) * F1(rhoidx)  ! F28 -- F1
    END IF
    EAM(ix) = (f_rho(ix) + 0.5 * phi_sum(ix) ) 
!PRINT*, 'EAM [ev] = ', EAM(ix)
  END DO 
  
END SUBROUTINE do_calcEAM
!##########################################################   
SUBROUTINE do_fdmDiff(dfx, f, dx)
!##########################################################    
  IMPLICIT NONE
  
  REAL(dp) :: f(:) &
            , dx
  
  REAL(dp), ALLOCATABLE :: dfx(:)
  INTEGER :: L 
  
  L = SIZE(f,1)

  ALLOCATE(dfx(1:L))
  
  dfx(:) = 0.0
  dfx(1:L-1) = (f(1:L-1) - f(2:L) ) / dx 
  dfx(L)=dfx(L-1)
!write(*,*), dfx(:)
END SUBROUTINE do_fdmDiff
!##########################################################
SUBROUTINE do_calc_drEAMParam(dr, drho, F1, rho1, F2, rho2,   & 
                              phi11, phi21, phi22,            & 
                              dFdr1, drhodr1, dFdr2, drhodr2, &
                              dphidr11, dphidr21, dphidr22)

  REAL(dp), ALLOCATABLE :: F1   (:)          & 
                         , rho1 (:)          &
                         , F2   (:)          &
                         , rho2 (:)          &
                         , phi11(:)          &
                         , phi21(:)          &
                         , phi22(:)          &
                         , dFdr1(:)          &
                         , drhodr1(:)        &
                         , dFdr2(:)          &
                         , drhodr2(:)        &
                         , dphidr11(:)       &
                         , dphidr21(:)       &
                         , dphidr22(:)       
                         
  REAL(dp) :: dr, drho
  
  CALL do_fdmDiff(dFdr1    ,F1   , drho)  
  CALL do_fdmDiff(drhodr1  ,rho1 , dr  ) 
  CALL do_fdmDiff(dFdr2    ,F2   , drho)     
  CALL do_fdmDiff(drhodr2  ,rho2 , dr  ) 
  CALL do_fdmDiff(dphidr11 ,phi11, dr  ) 
  CALL do_fdmDiff(dphidr21 ,phi21, dr  ) 
  CALL do_fdmDiff(dphidr22 ,phi22, dr  ) 
!write(*,*), dfdr1(:) 
!write(*,*), '*****' 
!write(*,*), drhodr1(:) 
!write(*,*), '*****' 
!write(*,*), dfdr2(:) 
!write(*,*), '*****' 
!write(*,*), drhodr2(:) 
!write(*,*), '*****' 
!write(*,*), dphidr11(:)
!write(*,*), '*****' 
!write(*,*), dphidr21(:)
!write(*,*), '*****' 
!write(*,*), dphidr22(:)

END SUBROUTINE do_calc_drEAMParam
!##########################################################

 
END MODULE utilities
