PROGRAM test_rand_normal2
  USE glbl
  USE util
  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 1000000
  REAL(dp) :: result_array(n)
  REAL(dp) :: mean, sd
  REAL(dp) :: expected_mean = 0.0_dp
  REAL(dp) :: expected_sd = 1.0_dp
  REAL(dp) :: tolerance = 0.05_dp

  ! Call the function to test
  result_array = rand_normal2(n)

  ! Calculate empirical mean and standard deviation
  mean = SUM(result_array) / n
  sd = SQRT(SUM((result_array - mean)**2) / n)

  WRITE(*, "(A,F8.6)") "Empirical Mean = ", mean
  WRITE(*, "(A,F8.6)") "Empirical Standard Deviation = ", sd

  ! Check if they are within tolerance
  IF (ABS(mean - expected_mean) > tolerance) THEN
    WRITE(*, *) "Test FAILED: Mean is out of tolerance."
    STOP 1
  END IF

  IF (ABS(sd - expected_sd) > tolerance) THEN
    WRITE(*, *) "Test FAILED: Standard Deviation is out of tolerance."
    STOP 1
  END IF

  WRITE(*, *) "Test PASSED: Mean and Standard Deviation are within tolerance."

END PROGRAM test_rand_normal2
