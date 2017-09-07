module md

use global
use utilities


contains

    
end module md





    
    !subroutine do_alloc_dtyp(tmstp)
    !
    !! allocate tmstp datatype vector
    !! allocate all vector fileds of all cells of tmstp vector
    !
    !  type(tmstp_holder), allocatable :: tmstp(:)
    !  allocate(tmstp(0:floor(t/dt)))                           
    ! 
    !! loops over all timesteps to allocate variables in 
    !! the tmstp elements(datatpe)
    !  DO i = 0,N_tstps 
    !    ALLOCATE(tmstp(i)%xyz(N_atms,3)) ! N*3 matrix
    !    ALLOCATE(tmstp(i)%vel(N_atms,3)) ! N*3 matrix
    !    ALLOCATE(tmstp(i)%acc(N_atms,3)) ! N*3 matrix
    !    ALLOCATE(tmstp(i)%frc(N_atms,3)) ! N*3 matrix
    !
    !    tmstp(i)%xyz(:,:) = 0 
    !    tmstp(i)%vel(:,:) = 0       
    !    tmstp(i)%acc(:,:) = 0
    !    tmstp(i)%frc(:,:) = 0
    !
    !  END DO
    !  
    !end subroutine do_alloc_dtyp
   
   !##########################################################
   !subroutine do_loop_tstps(tmstp)
   !
   !! loops over all timesteps(N_tstp=t/dt)
   !!   in tstp=0:
   !!   in 0<tstp<N_tstp
   !
   !
   !  TYPE(tmstp_holder), DIMENSION(:), ALLOCATABLE :: tmstp
   !  INTEGER :: i ,j , k
   !  !REAL(dp) :: r(3) , vstd, tempT, vcm , vcm_tmp(n_atms,3), vel_scale, tmp_scale
   !
   !  DO i= 0 , N_tstps
   !    IF (i == 0) then
   !                      
   !    ! randomize xyz of atoms in cell 0 
   !                            !print *, tmstp(i)%xyz    
   !      CALL do_rand_xyz(tmstp(i)%xyz)
   !          
   !    ! randomize velocities of atoms in cell zero
   !      CALL do_rand_vel(tmstp(i)%vel)          
   !                            !CALL do_calcT(tempT, tmstp(i)%vel)
   !                            !PRINT *, 'tempT= ', tempT,'T= ', tmp
   !      
   !    ! fix center of mass
   !      CALL do_fix_com(tmstp(i)%vel)
   !     
   !    ! Scale velocities with initial temperature
   !        CALL do_scale_vel(tmstp(i)%vel)
   !                             !print *,tmstp(i)%vel(1,:)
   !    ! ToDo: question? should values travers to i+1?       
   !
   !    ELSE
   !      ! TODO : velocityverlet?
   !      
   !      CALL do_velverlet( tmstp(i)%xyz, tmstp(i)%vel, tmstp(i)%acc, tmstp(i)%frc, atm_masses, dt )  
   !      
   !      !tmstp(i)%frc(1,:) =1+1.5* tmstp(i)%frc(1,:)
   !       
   !      !print *, tmstp(i)%vel(1,:)
   !
   !
   !      ! write tstp data to a file 
   !    END IF
   !                            print *, tmstp(i)%vel(1,:)
   !  END DO
   !  
   !END SUBROUTINE do_loop_tstps
   
   !##########################################################
