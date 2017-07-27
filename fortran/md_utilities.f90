MODULE utilities
  
  use datastructure
  !use random
    
  REAL, PRIVATE    :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0
  PRIVATE          :: integral
  
  CONTAINS  
   
  !##########################################################
   
   SUBROUTINE do_genRand(r,m,j)
    
    !ToDo : no newwd to pass m
    ! 
    ! 
 
     IMPLICIT NONE
     
     INTEGER :: i, j, m, n, clock
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
  
  FUNCTION random_normal(n)
  
  ! Adapted from the following Fortran 77 code
  !      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
  !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
  !      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
  
  !  The function random_normal() returns a normally distributed pseudo-random
  !  number with zero mean and unit variance.
  
  !  The algorithm uses the ratio of uniforms method of A.J. Kinderman
  !  and J.F. Monahan augmented with quadratic bounding curves.
  
    IMPLICIT NONE
  
  ! local variables
    INTEGER  :: n , i
    REAL(dp) :: random_normal(n)
    REAL(dp)     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
                r1 = 0.27597, r2 = 0.27846, u, v, x, y, q , mean , sd
  
  !     Generate P = (u,v) uniform in rectangle enclosing acceptance region
  
    DO i =  1 , n
      DO
        CALL RANDOM_NUMBER(u)
        CALL RANDOM_NUMBER(v)
        v = 1.7156 * (v - half)
    
    !     Evaluate the quadratic form
        x = u - s
        y = ABS(v) - t
        q = x**2 + y*(a*y - b*x)
    
    !     Accept P if inside inner ellipse
        IF (q < r1) EXIT
    !     Reject P if outside outer ellipse
        IF (q > r2) CYCLE
    !     Reject P if outside acceptance region
        IF (v**2 < -4.0*LOG(u)*u**2) EXIT
    
      END DO

      random_normal(i) = v/u

    END DO
  
  
  
  ! Check mean and standard deviation                 
    mean = SUM(random_normal)/n              
    sd = SQRT(SUM((random_normal - mean)**2)/n)
  
    WRITE(*, "(A,F8.6)") "Mean = ", mean            
    WRITE(*, "(A,F8.6)") "Standard Deviation = ", sd
    
    
  !     Return ratio of P's coordinates as the normal deviate
  
    RETURN
  
  END FUNCTION random_normal

  !########################################################## 
  
  FUNCTION rand_normal2(n) 
  
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
 
    RETURN
  END FUNCTION rand_normal2
  
  !########################################################## 
  SUBROUTINE do_calcT(t,vel)
    
    IMPLICIT NONE
    
    REAL(dp)::  vel(:,:), vel_sqr,t
    INTEGER :: dgr_frdm , i
    
    dgr_frdm = 3 * N_atms - 3
    t = 0
    DO i = 1 , N_atms
      vel_sqr = vel(i,1)**2 + vel(i,2)**2 + vel(i,3)**2 
      !print *, vel_sqr      
      t = t +  (atm_masses(i)* vel_sqr) / (kb * dgr_frdm) 
    END DO      
    
    RETURN
   
  END SUBROUTINE do_calcT




END MODULE utilities





































































