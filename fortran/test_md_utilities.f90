PROGRAM test_md_utilities
  USE glbl
  USE util
  IMPLICIT NONE

  REAL(dp), ALLOCATABLE :: v(:,:), m(:)
  REAL(dp) :: target_T, calculated_T, tolerance
  INTEGER :: n_particles, i

  tolerance = 1.0D-3

  ! --- Test Case 1: Happy Path ---
  WRITE(*,*) "Running Test Case 1: Happy Path"
  n_particles = 10
  ALLOCATE(v(n_particles, 3))
  ALLOCATE(m(n_particles))

  ! Initialize masses to 1.0
  m(:) = 1.0D0

  ! Initialize velocities to some values
  DO i = 1, n_particles
    v(i, 1) = 0.5D0 * i
    v(i, 2) = -0.2D0 * i
    v(i, 3) = 0.1D0 * i
  END DO

  target_T = 300.0D0 ! Room temperature

  CALL do_scaleVel(v, target_T, m)
  CALL do_calcT(calculated_T, v, m)

  WRITE(*,*) "Target T: ", target_T
  WRITE(*,*) "Calculated T: ", calculated_T

  IF (ABS(calculated_T - target_T) > tolerance) THEN
    WRITE(*,*) "TEST 1 FAILED!"
    STOP 1
  ELSE
    WRITE(*,*) "TEST 1 PASSED!"
  END IF

  DEALLOCATE(v, m)

  ! --- Test Case 2: Edge Case - Very Low Temperature ---
  WRITE(*,*) "Running Test Case 2: Very Low Temperature"
  n_particles = 5
  ALLOCATE(v(n_particles, 3))
  ALLOCATE(m(n_particles))

  ! Initialize masses to 1.0
  m(:) = 2.0D0

  ! Initialize velocities to some values
  DO i = 1, n_particles
    v(i, 1) = 1.0D0
    v(i, 2) = 2.0D0
    v(i, 3) = -1.0D0
  END DO

  target_T = 0.01D0 ! Very low temperature

  CALL do_scaleVel(v, target_T, m)
  CALL do_calcT(calculated_T, v, m)

  WRITE(*,*) "Target T: ", target_T
  WRITE(*,*) "Calculated T: ", calculated_T

  IF (ABS(calculated_T - target_T) > tolerance) THEN
    WRITE(*,*) "TEST 2 FAILED!"
    STOP 1
  ELSE
    WRITE(*,*) "TEST 2 PASSED!"
  END IF

  DEALLOCATE(v, m)

  ! --- Test Case 3: Edge Case - Very High Temperature ---
  WRITE(*,*) "Running Test Case 3: Very High Temperature"
  n_particles = 100
  ALLOCATE(v(n_particles, 3))
  ALLOCATE(m(n_particles))

  ! Initialize masses to 1.0
  m(:) = 0.5D0

  ! Initialize velocities to small values
  DO i = 1, n_particles
    v(i, 1) = 0.01D0
    v(i, 2) = 0.01D0
    v(i, 3) = 0.01D0
  END DO

  target_T = 5000.0D0 ! Very high temperature

  CALL do_scaleVel(v, target_T, m)
  CALL do_calcT(calculated_T, v, m)

  WRITE(*,*) "Target T: ", target_T
  WRITE(*,*) "Calculated T: ", calculated_T

  IF (ABS(calculated_T - target_T) > tolerance) THEN
    WRITE(*,*) "TEST 3 FAILED!"
    STOP 1
  ELSE
    WRITE(*,*) "TEST 3 PASSED!"
  END IF

  DEALLOCATE(v, m)

  WRITE(*,*) "ALL TESTS PASSED SUCCESSFULLY!"

END PROGRAM test_md_utilities
