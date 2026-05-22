PROGRAM test_cm
  USE glbl
  USE util
  IMPLICIT NONE

  INTEGER, PARAMETER :: np = 2
  REAL(dp) :: v(np, 3), m(np)
  REAL(dp) :: vcm(3)

  ! Initialize
  v(1,:) = [1.0_dp, 0.0_dp, 0.0_dp]
  v(2,:) = [-1.0_dp, 0.0_dp, 0.0_dp]
  m(1) = 2.0_dp
  m(2) = 1.0_dp

  CALL do_fix_centerOfMass(v, m)

  ! Verify total momentum is now zero
  vcm(1) = SUM(v(:,1) * m(:))
  vcm(2) = SUM(v(:,2) * m(:))
  vcm(3) = SUM(v(:,3) * m(:))

  IF (ABS(vcm(1)) > 1e-10 .OR. ABS(vcm(2)) > 1e-10 .OR. ABS(vcm(3)) > 1e-10) THEN
    WRITE(*,*) "TEST FAILED: Net momentum is not zero"
    WRITE(*,*) "vcm = ", vcm
    STOP 1
  ELSE
    WRITE(*,*) "TEST PASSED"
  END IF

END PROGRAM test_cm
