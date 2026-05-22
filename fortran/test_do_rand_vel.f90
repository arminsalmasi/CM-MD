PROGRAM test_do_rand_vel
  USE glbl
  USE util
  IMPLICIT NONE

  INTEGER, PARAMETER :: np = 10000
  REAL(dp) :: v(np, 3)
  REAL(dp) :: m(np)
  REAL(dp) :: T
  REAL(dp) :: expected_vstd
  REAL(dp) :: calculated_mean, calculated_std
  REAL(dp) :: v_flat(np * 3)
  INTEGER :: i, j

  ! Initialize variables
  T = 300.0_dp
  m(:) = 1.0e-26_dp

  expected_vstd = SQRT( T * kb / m(1) )

  ! Call subroutine to test
  CALL do_rand_vel(v, T, m)

  ! Flatten the array to calculate overall mean and std
  DO i = 1, np
    DO j = 1, 3
      v_flat((i-1)*3 + j) = v(i, j)
    END DO
  END DO

  calculated_mean = SUM(v_flat) / (np * 3)
  calculated_std = SQRT(SUM((v_flat - calculated_mean)**2) / (np * 3))

  WRITE(*, "(A,F15.6)") "Target std = ", expected_vstd
  WRITE(*, "(A,F15.6)") "Calculated mean = ", calculated_mean
  WRITE(*, "(A,F15.6)") "Calculated std = ", calculated_std

  ! Basic assertion
  ! The standard deviation of the mean for N samples is sigma / sqrt(N).
  ! We have 30000 samples, so standard error is roughly 643 / 173 = ~3.7.
  ! Setting the tolerance for the mean to 5% of std is ~32, which is very safe.
  IF (ABS(calculated_mean) > expected_vstd * 0.05_dp) THEN
    WRITE(*,*) "FAIL: Mean is too far from 0"
    STOP 1
  END IF

  IF (ABS(calculated_std - expected_vstd) > expected_vstd * 0.05_dp) THEN
    WRITE(*,*) "FAIL: Standard deviation does not match target"
    STOP 1
  END IF

  WRITE(*,*) "SUCCESS: do_rand_vel generates correct distribution"

END PROGRAM test_do_rand_vel
