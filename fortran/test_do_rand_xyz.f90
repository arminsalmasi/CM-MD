PROGRAM test_do_rand_xyz_main
  USE glbl
  USE util
  IMPLICIT NONE

  REAL(dp), ALLOCATABLE :: xyz1(:,:), xyz2(:,:)
  REAL(dp) :: V
  INTEGER :: nPart
  INTEGER :: ix, jx
  REAL(dp) :: diff

  nPart = 100
  V = 1000.0_dp

  ALLOCATE(xyz1(nPart, 3))
  ALLOCATE(xyz2(nPart, 3))

  ! Initialize
  xyz1 = 0.0_dp
  xyz2 = 0.0_dp

  ! Test 1: Boundaries
  CALL do_rand_xyz(xyz1, V)

  DO ix = 1, nPart
     DO jx = 1, 3
        IF (xyz1(ix, jx) < 0.0_dp .OR. xyz1(ix, jx) > V**(1.0_dp/3.0_dp)) THEN
           PRINT *, "Boundary test failed!", xyz1(ix, jx)
           STOP 1
        END IF
     END DO
  END DO
  PRINT *, "Boundary test passed!"

  ! Test 2: Seed
  CALL do_rand_xyz(xyz1, V, 42)
  CALL do_rand_xyz(xyz2, V, 42)

  DO ix = 1, nPart
     DO jx = 1, 3
        IF (ABS(xyz1(ix, jx) - xyz2(ix, jx)) > 1e-10) THEN
           PRINT *, "Seed test failed!", xyz1(ix, jx), xyz2(ix, jx)
           STOP 1
        END IF
     END DO
  END DO
  PRINT *, "Seed test passed!"

END PROGRAM test_do_rand_xyz_main
