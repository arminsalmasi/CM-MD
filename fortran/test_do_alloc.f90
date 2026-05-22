PROGRAM test_do_alloc
  USE glbl
  USE util
  IMPLICIT NONE

  INTEGER :: nPart, nSamp
  REAL(dp), ALLOCATABLE :: xyz(:,:), vel(:,:), acc(:,:), frc(:,:), EAM(:)
  TYPE(timestep_sample), ALLOCATABLE :: samp(:)
  INTEGER :: ix, jx
  LOGICAL :: passed

  ! Test parameters
  nPart = 10
  nSamp = 5
  passed = .TRUE.

  ! Call do_Alloc subroutine
  CALL do_Alloc(xyz, vel, frc, acc, EAM, samp, nSamp, nPart)

  ! Verify allocation and sizes
  IF (.NOT. ALLOCATED(xyz)) THEN
    PRINT *, "FAIL: xyz is not allocated"
    passed = .FALSE.
  ELSE IF (SIZE(xyz, 1) /= nPart .OR. SIZE(xyz, 2) /= 3) THEN
    PRINT *, "FAIL: xyz size is incorrect. Expected:", nPart, "x 3, got", SIZE(xyz, 1), "x", SIZE(xyz, 2)
    passed = .FALSE.
  ELSE IF (ANY(xyz /= 0D0)) THEN
    PRINT *, "FAIL: xyz is not initialized to zero"
    passed = .FALSE.
  END IF

  IF (.NOT. ALLOCATED(vel)) THEN
    PRINT *, "FAIL: vel is not allocated"
    passed = .FALSE.
  ELSE IF (SIZE(vel, 1) /= nPart .OR. SIZE(vel, 2) /= 3) THEN
    PRINT *, "FAIL: vel size is incorrect."
    passed = .FALSE.
  ELSE IF (ANY(vel /= 0D0)) THEN
    PRINT *, "FAIL: vel is not initialized to zero"
    passed = .FALSE.
  END IF

  IF (.NOT. ALLOCATED(acc)) THEN
    PRINT *, "FAIL: acc is not allocated"
    passed = .FALSE.
  ELSE IF (SIZE(acc, 1) /= nPart .OR. SIZE(acc, 2) /= 3) THEN
    PRINT *, "FAIL: acc size is incorrect."
    passed = .FALSE.
  ELSE IF (ANY(acc /= 0D0)) THEN
    PRINT *, "FAIL: acc is not initialized to zero"
    passed = .FALSE.
  END IF

  IF (.NOT. ALLOCATED(frc)) THEN
    PRINT *, "FAIL: frc is not allocated"
    passed = .FALSE.
  ELSE IF (SIZE(frc, 1) /= nPart .OR. SIZE(frc, 2) /= 3) THEN
    PRINT *, "FAIL: frc size is incorrect."
    passed = .FALSE.
  ELSE IF (ANY(frc /= 0D0)) THEN
    PRINT *, "FAIL: frc is not initialized to zero"
    passed = .FALSE.
  END IF

  IF (.NOT. ALLOCATED(EAM)) THEN
    PRINT *, "FAIL: EAM is not allocated"
    passed = .FALSE.
  ELSE IF (SIZE(EAM) /= nPart) THEN
    PRINT *, "FAIL: EAM size is incorrect."
    passed = .FALSE.
  ELSE IF (ANY(EAM /= 0D0)) THEN
    PRINT *, "FAIL: EAM is not initialized to zero"
    passed = .FALSE.
  END IF

  IF (.NOT. ALLOCATED(samp)) THEN
    PRINT *, "FAIL: samp is not allocated"
    passed = .FALSE.
  ELSE IF (LBOUND(samp, 1) /= 0 .OR. UBOUND(samp, 1) /= nSamp) THEN
    PRINT *, "FAIL: samp bounds are incorrect. Expected: 0 to", nSamp, ", got", LBOUND(samp, 1), "to", UBOUND(samp, 1)
    passed = .FALSE.
  ELSE
    ! Check inside samp elements
    DO ix = 0, nSamp
      IF (.NOT. ALLOCATED(samp(ix)%xyz) .OR. SIZE(samp(ix)%xyz, 1) /= nPart .OR. SIZE(samp(ix)%xyz, 2) /= 3) THEN
        PRINT *, "FAIL: samp(", ix, ")%xyz allocation/size is incorrect"
        passed = .FALSE.
      END IF
      IF (ANY(samp(ix)%xyz /= 0D0)) THEN
        PRINT *, "FAIL: samp(", ix, ")%xyz is not zero"
        passed = .FALSE.
      END IF

      IF (.NOT. ALLOCATED(samp(ix)%EAM) .OR. SIZE(samp(ix)%EAM) /= nPart) THEN
        PRINT *, "FAIL: samp(", ix, ")%EAM allocation/size is incorrect"
        passed = .FALSE.
      END IF
      IF (ANY(samp(ix)%EAM /= 0D0)) THEN
        PRINT *, "FAIL: samp(", ix, ")%EAM is not zero"
        passed = .FALSE.
      END IF

      IF (samp(ix)%kin /= 0D0) THEN
         PRINT *, "FAIL: samp(", ix, ")%kin is not zero"
         passed = .FALSE.
      END IF
    END DO
  END IF

  IF (passed) THEN
    PRINT *, "do_Alloc tests PASSED"
  ELSE
    PRINT *, "do_Alloc tests FAILED"
    CALL EXIT(1)
  END IF

END PROGRAM test_do_alloc
