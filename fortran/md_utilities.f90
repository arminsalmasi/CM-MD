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
			   
					!write(*,*) xyz(:,:)
     DO j = 1, np
       call do_genRand(r,3,j)
					!write(*,*) 'r = ', r
 	   xyz(j,:) = r(:) * ( V **(1.0/3.0))
					!write(*,*) 'xyz = ', xyz(j,:)
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
									!veltmp(:,j) = vel(:,j) 
     END DO
     vcm = SUM(vcm_tmp) / SUM(masses)
                   				   !print *,vcm								   
	   vel(:,1:3) = vel(:,1:3) - vcm
								   !print*, vel-veltmp
    
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
								!print*, Ttemp
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
SUBROUTINE do_velVerlet(np, u , xyz, vel, acc, frc, masses, dt);
!##########################################################
     
     IMPLICIT NONE              
     
     INTEGER :: i  &
	          , np
			  
     REAL(dp) :: xyz(:,:)  &
	           , vel(:,:)  &
	           , acc(:,:)  &
	           , frc(:,:)  &
	           , masses(:) &
	           , dt 	 
	 
	 
	 DO i=1 , 3 
       xyz(:,i) = xyz(:,i) + vel(:,i) * dt + 0.5 * acc(:,i) * dt * dt;
     END DO   
   
     DO i = 1 , 3
        vel(:,i) = vel(:,i) + 0.5 * ( frc(:,i) / masses(:) + acc(:,i)) * dt;
     END DO
   
     DO i = 1 , 3 
        acc(:,i) =  frc(:,i) / masses(:);
     END DO
   
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
    
  ! CALL RANDOM_NUMBER(array) ! Uniform distribution
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
    
    REAL(dp)::  vel(:,:), vel_sqr(np), temp, masses(:)
    INTEGER :: dof , i, np
   
    vel_sqr(:) = vel(:,1)**2 + vel(:,2)**2 + vel(:,3)**2 
							!do i = 1, size(vel_sqr)
							!	print *, vel(i,:)
							!	print *, vel_sqr(i)      
						    !end do
	dof = 3 * np - 3
	temp = sum(masses(:) * vel_sqr(:))
    temp = temp / (kb * dof) 
							!write(*,*) temp   
  END SUBROUTINE do_calcT
!########################################################## 

END MODULE utilities







 
!  FUNCTION random_normal(n)
!  
!  ! Adapted from the following Fortran 77 code
!  !      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!  !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!  !      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
!  
!  !  The function random_normal() returns a normally distributed pseudo-random
!  !  number with zero mean and unit variance.
!  
!  !  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  !  and J.F. Monahan augmented with quadratic bounding curves.
!  
!    IMPLICIT NONE
!  
!  ! local variables
!    INTEGER  :: n , i
!    REAL(dp) :: random_normal(n)
!    REAL(dp)     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
!                r1 = 0.27597, r2 = 0.27846, u, v, x, y, q , mean , sd
!  
!  !     Generate P = (u,v) uniform in rectangle enclosing acceptance region
!  
!    DO i =  1 , n
!      DO
!        CALL RANDOM_NUMBER(u)
!        CALL RANDOM_NUMBER(v)
!        v = 1.7156 * (v - half)
!    
!    !     Evaluate the quadratic form
!        x = u - s
!        y = ABS(v) - t
!        q = x**2 + y*(a*y - b*x)
!    
!    !     Accept P if inside inner ellipse
!        IF (q < r1) EXIT
!    !     Reject P if outside outer ellipse
!        IF (q > r2) CYCLE
!    !     Reject P if outside acceptance region
!        IF (v**2 < -4.0*LOG(u)*u**2) EXIT
!    
!      END DO
!
!      random_normal(i) = v/u
!
!    END DO
!  
!  
!  
!  ! Check mean and standard deviation                 
!    mean = SUM(random_normal)/n              
!    sd = SQRT(SUM((random_normal - mean)**2)/n)
!  
!    WRITE(*, "(A,F8.6)") "Mean = ", mean            
!    WRITE(*, "(A,F8.6)") "Standard Deviation = ", sd
!    
!    
!  !     Return ratio of P's coordinates as the normal deviate
!  
!    RETURN
!  
!  END FUNCTION random_normal
!
!  !##########################################################





























































