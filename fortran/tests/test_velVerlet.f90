PROGRAM test_velVerlet
  USE glbl
  USE util
  IMPLICIT NONE

  REAL(dp), ALLOCATABLE :: xyz(:,:), v(:,:), a(:,:), f(:,:), m(:)
  REAL(dp), ALLOCATABLE :: xyz_expected(:,:), v_expected(:,:), a_expected(:,:)
  REAL(dp) :: dt
  INTEGER :: n_atoms
  REAL(dp) :: tol = 1.0e-12_dp
  INTEGER :: ix, jx
  LOGICAL :: passed

  n_atoms = 2
  dt = 0.1_dp

  ALLOCATE(xyz(n_atoms, 3), v(n_atoms, 3), a(n_atoms, 3), f(n_atoms, 3), m(n_atoms))
  ALLOCATE(xyz_expected(n_atoms, 3), v_expected(n_atoms, 3), a_expected(n_atoms, 3))

  ! Initialize state
  ! Atom 1
  xyz(1,:) = [0.0_dp, 0.0_dp, 0.0_dp]
  v(1,:) = [1.0_dp, 0.0_dp, 0.0_dp]
  a(1,:) = [2.0_dp, 0.0_dp, 0.0_dp]
  f(1,:) = [4.0_dp, 0.0_dp, 0.0_dp]
  m(1) = 2.0_dp

  ! Atom 2
  xyz(2,:) = [1.0_dp, 1.0_dp, 1.0_dp]
  v(2,:) = [0.0_dp, 1.0_dp, 0.0_dp]
  a(2,:) = [0.0_dp, 3.0_dp, 0.0_dp]
  f(2,:) = [0.0_dp, 9.0_dp, 0.0_dp]
  m(2) = 3.0_dp

  ! Expected values after velocity Verlet step
  ! xyz_new = xyz + v*dt + 0.5*a*dt^2
  ! a_new = f/m
  ! v_new = v + 0.5*(a + a_new)*dt

  ! Atom 1 Expected
  xyz_expected(1,:) = [0.0_dp + 1.0_dp*0.1_dp + 0.5_dp*2.0_dp*0.01_dp, 0.0_dp, 0.0_dp] ! [0.11, 0, 0]
  a_expected(1,:) = [4.0_dp/2.0_dp, 0.0_dp, 0.0_dp] ! [2.0, 0, 0]
  v_expected(1,:) = [1.0_dp + 0.5_dp*(2.0_dp + 2.0_dp)*0.1_dp, 0.0_dp, 0.0_dp] ! [1.2, 0, 0]

  ! Atom 2 Expected
  xyz_expected(2,:) = [1.0_dp, 1.0_dp + 1.0_dp*0.1_dp + 0.5_dp*3.0_dp*0.01_dp, 1.0_dp] ! [1.0, 1.115, 1.0]
  a_expected(2,:) = [0.0_dp, 9.0_dp/3.0_dp, 0.0_dp] ! [0.0, 3.0, 0.0]
  v_expected(2,:) = [0.0_dp, 1.0_dp + 0.5_dp*(3.0_dp + 3.0_dp)*0.1_dp, 0.0_dp] ! [0.0, 1.3, 0.0]

  ! Call the subroutine
  CALL do_velVerlet(xyz, v, a, f, m, dt)

  ! Check the results
  passed = .TRUE.

  DO ix = 1, n_atoms
    DO jx = 1, 3
      IF (ABS(xyz(ix, jx) - xyz_expected(ix, jx)) > tol) THEN
        PRINT *, "Mismatch in xyz for atom ", ix, " component ", jx
        PRINT *, "Expected: ", xyz_expected(ix, jx), " Got: ", xyz(ix, jx)
        passed = .FALSE.
      END IF

      IF (ABS(v(ix, jx) - v_expected(ix, jx)) > tol) THEN
        PRINT *, "Mismatch in v for atom ", ix, " component ", jx
        PRINT *, "Expected: ", v_expected(ix, jx), " Got: ", v(ix, jx)
        passed = .FALSE.
      END IF

      IF (ABS(a(ix, jx) - a_expected(ix, jx)) > tol) THEN
        PRINT *, "Mismatch in a for atom ", ix, " component ", jx
        PRINT *, "Expected: ", a_expected(ix, jx), " Got: ", a(ix, jx)
        passed = .FALSE.
      END IF
    END DO
  END DO

  IF (passed) THEN
    PRINT *, "Test passed!"
  ELSE
    PRINT *, "Test failed."
    STOP 1
  END IF

END PROGRAM test_velVerlet
