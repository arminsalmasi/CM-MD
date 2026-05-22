PROGRAM test_calcEAMPot
  USE glbl
  USE util
  IMPLICIT NONE

  REAL(dp) :: U
  INTEGER :: ix
  REAL(dp), ALLOCATABLE :: xyz(:,:), m(:)
  INTEGER, ALLOCATABLE :: AtN(:), atNVec(:), nAt(:)
  TYPE(EAM_data) :: EAMData
  INTEGER :: test_fails

  test_fails = 0

  ! Initialize EAMData
  EAMData%nR = 10
  EAMData%nRho = 10
  EAMData%dr = 0.5_dp
  EAMData%drho = 0.5_dp
  EAMData%cutoff = 4.0_dp

  ALLOCATE(EAMData%F1(10))
  ALLOCATE(EAMData%rho1(10))
  ALLOCATE(EAMData%F2(10))
  ALLOCATE(EAMData%rho2(10))
  ALLOCATE(EAMData%phi11(10))
  ALLOCATE(EAMData%phi21(10))
  ALLOCATE(EAMData%phi22(10))

  EAMData%F1 = 1.0_dp
  EAMData%rho1 = 1.0_dp
  EAMData%F2 = 2.0_dp
  EAMData%rho2 = 2.0_dp
  EAMData%phi11 = 1.0_dp
  EAMData%phi21 = 2.0_dp
  EAMData%phi22 = 3.0_dp

  ALLOCATE(xyz(2,3))
  ALLOCATE(m(2))
  ALLOCATE(AtN(2))
  ALLOCATE(nAt(2))
  ALLOCATE(atNVec(2))

  AtN(1) = 27
  AtN(2) = 28

  ! --- TEST 1: Two Co atoms interacting ---
  nAt(1) = 2
  nAt(2) = 0
  atNVec(1) = 27
  atNVec(2) = 27

  xyz(1,:) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
  xyz(2,:) = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
  m = 1.0_dp

  CALL do_calcEAMPot( U, 1, xyz, m, AtN, atNVec, nAt, EAMData)

  IF (ABS(U - 3.5_dp) > 1e-6) THEN
    PRINT *, "FAIL: Test 1 (Co-Co interaction). Expected 3.5, got ", U
    test_fails = test_fails + 1
  ELSE
    PRINT *, "PASS: Test 1 (Co-Co interaction)"
  END IF

  ! --- TEST 2: Two Ni atoms interacting ---
  nAt(1) = 0
  nAt(2) = 2
  atNVec(1) = 28
  atNVec(2) = 28

  CALL do_calcEAMPot( U, 1, xyz, m, AtN, atNVec, nAt, EAMData)

  IF (ABS(U - 1.5_dp) > 1e-6) THEN
    PRINT *, "FAIL: Test 2 (Ni-Ni interaction). Expected 1.5, got ", U
    test_fails = test_fails + 1
  ELSE
    PRINT *, "PASS: Test 2 (Ni-Ni interaction)"
  END IF

  ! --- TEST 3: Co and Ni interacting ---
  nAt(1) = 1
  nAt(2) = 1
  atNVec(1) = 27
  atNVec(2) = 28

  CALL do_calcEAMPot( U, 1, xyz, m, AtN, atNVec, nAt, EAMData)

  IF (ABS(U - 3.0_dp) > 1e-6) THEN
    PRINT *, "FAIL: Test 3 (Co-Ni interaction). Expected 3.0, got ", U
    test_fails = test_fails + 1
  ELSE
    PRINT *, "PASS: Test 3 (Co-Ni interaction)"
  END IF

  ! --- TEST 4: Three atoms, one outside cutoff ---
  DEALLOCATE(xyz)
  DEALLOCATE(m)
  DEALLOCATE(atNVec)
  ALLOCATE(xyz(3,3))
  ALLOCATE(m(3))
  ALLOCATE(atNVec(3))

  atNVec(1) = 27
  atNVec(2) = 28
  atNVec(3) = 27
  xyz(1,:) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
  xyz(2,:) = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
  xyz(3,:) = (/ 5.0_dp, 0.0_dp, 0.0_dp /)
  m = 1.0_dp

  CALL do_calcEAMPot( U, 1, xyz, m, AtN, atNVec, nAt, EAMData)

  IF (ABS(U - 3.0_dp) > 1e-6) THEN
    PRINT *, "FAIL: Test 4 (Outside cutoff). Expected 3.0, got ", U
    test_fails = test_fails + 1
  ELSE
    PRINT *, "PASS: Test 4 (Outside cutoff)"
  END IF

  ! --- TEST 5: Same atom test (ix=jx should be ignored) ---
  DEALLOCATE(xyz)
  DEALLOCATE(m)
  DEALLOCATE(atNVec)
  ALLOCATE(xyz(1,3))
  ALLOCATE(m(1))
  ALLOCATE(atNVec(1))

  atNVec(1) = 27
  xyz(1,:) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
  m = 1.0_dp

  CALL do_calcEAMPot( U, 1, xyz, m, AtN, atNVec, nAt, EAMData)

  IF (ABS(U - 0.0_dp) > 1e-6) THEN
    PRINT *, "FAIL: Test 5 (Self interaction). Expected 0.0, got ", U
    test_fails = test_fails + 1
  ELSE
    PRINT *, "PASS: Test 5 (Self interaction)"
  END IF

  IF (test_fails > 0) THEN
    PRINT *, "Tests failed: ", test_fails
    STOP 1
  ELSE
    PRINT *, "All tests passed!"
  END IF

END PROGRAM test_calcEAMPot
