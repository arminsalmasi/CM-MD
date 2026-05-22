PROGRAM test_md_utilities
  USE glbl
  USE util
  IMPLICIT NONE

  REAL(dp), ALLOCATABLE :: v(:,:)
  REAL(dp), ALLOCATABLE :: m(:)
  REAL(dp) :: T_calc, K_calc, T_expected, rel_err
  INTEGER  :: dof, np

  np = 2
  ALLOCATE(v(np,3))
  ALLOCATE(m(np))

  v = 0.0_dp
  m = 0.0_dp

  m(1) = 2.0_dp
  m(2) = 3.0_dp
  v(1,:) = [1.0_dp, 2.0_dp, 3.0_dp]
  v(2,:) = [4.0_dp, 5.0_dp, 6.0_dp]

  ! Calculate T directly using the subroutine
  CALL do_calcT(T_calc, v, m)

  ! Calculate K inline to verify against
  ! From do_calcKin equivalent:
  K_calc = 0.5_dp * sum(m * (v(:,1)**2 + v(:,2)**2 + v(:,3)**2))

  dof = 3 * np - 3
  T_expected = (2.0_dp * K_calc) / (kb * REAL(dof, dp))

  rel_err = ABS(T_calc - T_expected) / MAX(ABS(T_expected), 1e-10_dp)

  IF (rel_err < 1e-10_dp) THEN
     PRINT *, "PASS: do_calcT matches kinetic energy definition"
  ELSE
     PRINT *, "FAIL: do_calcT does not match kinetic energy definition"
     PRINT *, "Calculated T: ", T_calc
     PRINT *, "Expected T:   ", T_expected
     STOP 1
  END IF

  ! Test case for v = 0
  v = 0.0_dp
  CALL do_calcT(T_calc, v, m)
  IF (ABS(T_calc) < 1e-15_dp) THEN
     PRINT *, "PASS: do_calcT returns 0 for zero velocities"
  ELSE
     PRINT *, "FAIL: do_calcT should return 0 for zero velocities"
     STOP 1
  END IF

  ! Test case for a larger array with random values
  ! Since there's no random generator available easily,
  ! we just use deterministic varying data
  DEALLOCATE(v)
  DEALLOCATE(m)
  np = 100
  ALLOCATE(v(np,3))
  ALLOCATE(m(np))

  m = 10.0_dp
  v(:,1) = 5.0_dp
  v(:,2) = -2.0_dp
  v(:,3) = 1.0_dp

  CALL do_calcT(T_calc, v, m)
  K_calc = 0.5_dp * sum(m * (v(:,1)**2 + v(:,2)**2 + v(:,3)**2))
  dof = 3 * np - 3
  T_expected = (2.0_dp * K_calc) / (kb * REAL(dof, dp))

  rel_err = ABS(T_calc - T_expected) / MAX(ABS(T_expected), 1e-10_dp)

  IF (rel_err < 1e-10_dp) THEN
     PRINT *, "PASS: do_calcT matches K definition for larger system"
  ELSE
     PRINT *, "FAIL: do_calcT does not match K definition for larger system"
     PRINT *, "Calculated T: ", T_calc
     PRINT *, "Expected T:   ", T_expected
     STOP 1
  END IF

  PRINT *, "All tests passed!"

END PROGRAM test_md_utilities
