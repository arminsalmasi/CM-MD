PROGRAM test_do_genRand
  USE glbl
  USE util
  IMPLICIT NONE

  INTEGER, PARAMETER :: m = 100000
  REAL(dp) :: r1(m), r2(m), r3(m), r4(1)
  REAL(dp) :: mean, variance, expected_mean, expected_var
  LOGICAL :: passed = .TRUE.
  INTEGER :: i
  REAL :: x

  PRINT *, "-------------------------------------------"
  PRINT *, "Running tests for do_genRand..."
  PRINT *, "-------------------------------------------"

  ! Test 1: Sequence Repeatability (Deterministic seed using j=0)
  ! With j=0, clock * j is always 0, eliminating system clock variation
  CALL do_genRand(r1, m, 0)
  CALL do_genRand(r2, m, 0)

  IF (ANY(r1 /= r2)) THEN
    PRINT *, "FAIL: Repeatability test failed: sequences are different for j=0"
    passed = .FALSE.
  ELSE
    PRINT *, "PASS: Repeatability test passed: sequences are identical for j=0"
  END IF

  ! Test 2: Sequence Non-Repeatability (Different seeds using j>0 and delay)
  CALL do_genRand(r1, m, 1)
  ! Small delay to ensure clock ticks
  x = 0.0
  DO i = 1, 10000000
    x = x + 1.0
  END DO
  CALL do_genRand(r3, m, 2)

  IF (ALL(r1 == r3)) THEN
    PRINT *, "FAIL: Seed variation test failed: sequences are identical for different j and time"
    passed = .FALSE.
  ELSE
    PRINT *, "PASS: Seed variation test passed: sequences are different for different j and time"
  END IF

  ! Test 3: Value Range (Values should be strictly within [0.0, 1.0))
  IF (ANY(r1 < 0.0_dp) .OR. ANY(r1 >= 1.0_dp)) THEN
    PRINT *, "FAIL: Value range test failed: values outside [0.0, 1.0)"
    passed = .FALSE.
  ELSE
    PRINT *, "PASS: Value range test passed: all values within [0.0, 1.0)"
  END IF

  ! Test 4: Distribution Properties (Uniform [0, 1))
  mean = SUM(r1) / m
  variance = SUM((r1 - mean)**2) / m
  expected_mean = 0.5_dp
  expected_var = 1.0_dp / 12.0_dp

  ! Use a tolerance for statistical tests
  IF (ABS(mean - expected_mean) > 0.01_dp) THEN
    PRINT *, "FAIL: Distribution test failed: mean is not ~0.5 (", mean, ")"
    passed = .FALSE.
  ELSE
    PRINT *, "PASS: Distribution test passed: mean is ~0.5"
  END IF

  IF (ABS(variance - expected_var) > 0.01_dp) THEN
    PRINT *, "FAIL: Distribution test failed: variance is not ~1/12 (", variance, ")"
    passed = .FALSE.
  ELSE
    PRINT *, "PASS: Distribution test passed: variance is ~1/12"
  END IF

  ! Test 5: Edge case m=1
  CALL do_genRand(r4, 1, 0)
  IF (r4(1) < 0.0_dp .OR. r4(1) >= 1.0_dp) THEN
    PRINT *, "FAIL: Edge case m=1 test failed"
    passed = .FALSE.
  ELSE
    PRINT *, "PASS: Edge case m=1 test passed"
  END IF

  PRINT *, "-------------------------------------------"
  IF (passed) THEN
    PRINT *, "ALL do_genRand TESTS PASSED"
  ELSE
    PRINT *, "SOME do_genRand TESTS FAILED"
    STOP 1
  END IF

END PROGRAM test_do_genRand
