PROGRAM test_do_calcKin
  USE glbl
  USE util
  IMPLICIT NONE

  REAL(dp) :: k
  REAL(dp), ALLOCATABLE :: m(:), v(:,:)
  LOGICAL :: all_passed

  all_passed = .TRUE.

  PRINT *, "========================================="
  PRINT *, "Running tests for do_calcKin"
  PRINT *, "========================================="

  ! Test 1: Zero velocity
  ALLOCATE(m(2))
  ALLOCATE(v(2,3))
  m(1) = 2.0_dp
  m(2) = 3.0_dp
  v = 0.0_dp
  CALL do_calcKin(k, m, v)
  IF (ABS(k) < 1e-10) THEN
    PRINT *, "  [PASS] Test 1: Zero velocity"
  ELSE
    PRINT *, "  [FAIL] Test 1: Zero velocity. Expected 0.0, got", k
    all_passed = .FALSE.
  END IF
  DEALLOCATE(m, v)

  ! Test 2: Basic functionality (positive velocities)
  ALLOCATE(m(2))
  ALLOCATE(v(2,3))
  m(1) = 2.0_dp
  v(1,1) = 1.0_dp
  v(1,2) = 2.0_dp
  v(1,3) = 3.0_dp

  m(2) = 3.0_dp
  v(2,1) = 2.0_dp
  v(2,2) = 0.0_dp
  v(2,3) = 1.0_dp
  CALL do_calcKin(k, m, v)
  ! Expected = 0.5*2*(1^2 + 2^2 + 3^2) + 0.5*3*(2^2 + 0^2 + 1^2)
  !          = 1 * 14 + 1.5 * 5 = 14 + 7.5 = 21.5
  IF (ABS(k - 21.5_dp) < 1e-10) THEN
    PRINT *, "  [PASS] Test 2: Basic functionality"
  ELSE
    PRINT *, "  [FAIL] Test 2: Basic functionality. Expected 21.5, got", k
    all_passed = .FALSE.
  END IF
  DEALLOCATE(m, v)

  ! Test 3: Negative velocities (squared, so should be same as positive)
  ALLOCATE(m(2))
  ALLOCATE(v(2,3))
  m(1) = 2.0_dp
  v(1,1) = -1.0_dp
  v(1,2) = -2.0_dp
  v(1,3) = -3.0_dp

  m(2) = 3.0_dp
  v(2,1) = -2.0_dp
  v(2,2) = 0.0_dp
  v(2,3) = -1.0_dp
  CALL do_calcKin(k, m, v)
  IF (ABS(k - 21.5_dp) < 1e-10) THEN
    PRINT *, "  [PASS] Test 3: Negative velocities"
  ELSE
    PRINT *, "  [FAIL] Test 3: Negative velocities. Expected 21.5, got", k
    all_passed = .FALSE.
  END IF
  DEALLOCATE(m, v)

  ! Test 4: Larger system
  ALLOCATE(m(10))
  ALLOCATE(v(10,3))
  m = 1.0_dp
  v = 2.0_dp
  CALL do_calcKin(k, m, v)
  ! Each particle KE = 0.5 * 1.0 * (2^2 + 2^2 + 2^2) = 0.5 * 12 = 6.0
  ! Total KE for 10 particles = 60.0
  IF (ABS(k - 60.0_dp) < 1e-10) THEN
    PRINT *, "  [PASS] Test 4: Larger system"
  ELSE
    PRINT *, "  [FAIL] Test 4: Larger system. Expected 60.0, got", k
    all_passed = .FALSE.
  END IF
  DEALLOCATE(m, v)

  PRINT *, "========================================="
  IF (all_passed) THEN
    PRINT *, "ALL TESTS PASSED"
  ELSE
    PRINT *, "SOME TESTS FAILED"
    CALL EXIT(1)
  END IF

END PROGRAM test_do_calcKin
