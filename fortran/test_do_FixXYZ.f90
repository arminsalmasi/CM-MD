PROGRAM test_do_FixXYZ
  USE glbl
  USE util
  IMPLICIT NONE

  REAL(dp) :: xyz(3), V, L
  INTEGER :: fail_count

  fail_count = 0
  V = 1000.0_dp
  L = V ** (1.0_dp / 3.0_dp) ! L = 10

  PRINT *, "Testing do_FixXYZ with V=", V, ", L=", L

  ! Case 1: Particle inside the box (no folding expected)
  xyz = [2.0_dp, -3.0_dp, 4.0_dp]
  CALL do_FixXYZ(xyz, V)
  IF (ABS(xyz(1) - 2.0_dp) > 1E-5 .OR. ABS(xyz(2) - (-3.0_dp)) > 1E-5 .OR. ABS(xyz(3) - 4.0_dp) > 1E-5) THEN
    PRINT *, "FAIL: Case 1 - Particle inside box moved. Got:", xyz
    fail_count = fail_count + 1
  ELSE
    PRINT *, "PASS: Case 1"
  END IF

  ! Case 2: Particle slightly outside the box, > L/2 (needs single folding)
  ! Expected: 6.0 -> 6.0 - 10.0 = -4.0
  xyz = [6.0_dp, 0.0_dp, 0.0_dp]
  CALL do_FixXYZ(xyz, V)
  IF (ABS(xyz(1) - (-4.0_dp)) > 1E-5) THEN
    PRINT *, "FAIL: Case 2 - Single folding > L/2 failed. Got:", xyz
    fail_count = fail_count + 1
  ELSE
    PRINT *, "PASS: Case 2"
  END IF

  ! Case 3: Particle slightly outside the box, < -L/2 (needs single folding)
  ! Expected: -7.0 -> -7.0 + 10.0 = 3.0
  xyz = [0.0_dp, -7.0_dp, 0.0_dp]
  CALL do_FixXYZ(xyz, V)
  IF (ABS(xyz(2) - 3.0_dp) > 1E-5) THEN
    PRINT *, "FAIL: Case 3 - Single folding < -L/2 failed. Got:", xyz
    fail_count = fail_count + 1
  ELSE
    PRINT *, "PASS: Case 3"
  END IF

  ! Case 4: Particle far outside the box, > L/2 (recursive folding)
  ! Expected: 17.0 -> 17.0 - 10.0 = 7.0 -> 7.0 - 10.0 = -3.0
  xyz = [0.0_dp, 0.0_dp, 17.0_dp]
  CALL do_FixXYZ(xyz, V)
  IF (ABS(xyz(3) - (-3.0_dp)) > 1E-5) THEN
    PRINT *, "FAIL: Case 4 - Recursive folding > L/2 failed. Got:", xyz
    fail_count = fail_count + 1
  ELSE
    PRINT *, "PASS: Case 4"
  END IF

  ! Case 5: Particle far outside the box, < -L/2 (recursive folding)
  ! Expected: -24.0 -> -24.0 + 10.0 = -14.0 -> -14.0 + 10.0 = -4.0
  xyz = [-24.0_dp, 0.0_dp, 0.0_dp]
  CALL do_FixXYZ(xyz, V)
  IF (ABS(xyz(1) - (-4.0_dp)) > 1E-5) THEN
    PRINT *, "FAIL: Case 5 - Recursive folding < -L/2 failed. Got:", xyz
    fail_count = fail_count + 1
  ELSE
    PRINT *, "PASS: Case 5"
  END IF

  IF (fail_count == 0) THEN
    PRINT *, "ALL TESTS PASSED!"
  ELSE
    PRINT *, fail_count, " TESTS FAILED."
    STOP 1
  END IF

END PROGRAM test_do_FixXYZ
