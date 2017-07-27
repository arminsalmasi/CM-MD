module md

use datastructure
use utilities


contains
!##########################################################
    
    subroutine do_allocate_dtyp(tmstp)
    
    ! allocate tmstp datatype vector
    ! allocate all vector fileds of all cells of tmstp vector
    
      type(tmstp_holder), allocatable :: tmstp(:)
      allocate(tmstp(0:floor(t/dt)))                           
     
    ! loops over all timesteps to allocate variables in 
    ! the tmstp elements(datatpe)
      do i = 0,N_tstps 
        allocate(tmstp(i)%xyz(N_atms,3)) ! N*3 matrix
        allocate(tmstp(i)%vel(N_atms,3)) ! N*3 matrix
        allocate(tmstp(i)%acc(N_atms,3)) ! N*3 matrix
      end do
      
    end subroutine do_allocate_dtyp
   
   !##########################################################
   
   subroutine do_loop_tstps(tmstp)
  
   ! loops over all timesteps(N_tstp=t/dt)
   !   in tstp=0:
   !   in 0<tstp<N_tstp
   
   
     TYPE(tmstp_holder), DIMENSION(:), ALLOCATABLE :: tmstp
     INTEGER :: i ,j , k
     REAL(dp) :: r(3) , vstd, tempT, vcm , vcm_tmp(n_atms,3), vel_scale, tmp_scale
 
     DO i= 0 , N_tstps
       IF (i == 0) then
                         
         ! randomize xyz of atoms in cell 0 
                               !print *, tmstp(i)%xyz    
         CALL do_rand_xyz(tmstp(i)%xyz)
             
         ! randomize velocities of atoms in cell zero
         CALL do_rand_vel(tmstp(i)%vel)          
                               !CALL do_calcT(tempT, tmstp(i)%vel)
                               !PRINT *, 'tempT= ', tempT,'T= ', tmp
         
         ! fix center of mass
                               !print *,tmstp(i)%vel(1,:)
         CALL do_fix_com(tmstp(i)%vel)
                               !print *,tmstp(i)%vel(1,:)
        
         ! Scale velocities with initial temperature
           CALL do_scale_vel(tmstp(i)%vel)
       ELSE
         ! TODO : do md on each timestep
         ! do_md()
!         print *, tmstp(i)%xyz(j,:)            
         ! save_tstp() 
       END IF
     END DO
     
   END SUBROUTINE do_loop_tstps
   
   !##########################################################
   
   SUBROUTINE do_rand_xyz(xyz)
   
     INTEGER :: j
     REAL(dp) :: r(3), xyz(:,:)
     
     DO j = 1, N_atms
    
       call do_genRand(r,3,j)
       xyz(j,:) = r(:) * ( box_vol**(1.0/3.0))
   
     END DO
   
   END SUBROUTINE do_rand_xyz
   
   !##########################################################
   
   SUBROUTINE do_rand_vel(vel)
     
     IMPLICIT NONE
     
     REAL(dp)  :: vel(:,:) , r(3), vstd 
     INTEGER :: j
     
     DO j = 1 , N_atms
       
       r(:) = rand_normal2(3) ! open source from roseta code
       vstd = sqrt( tmp * kb / atm_masses(j) )
       vel(j,:) = 0 + r(:) * vstd  
     
     END Do
   
   END SUBROUTINE do_rand_vel
   
   !##########################################################
   SUBROUTINE do_fix_com(vel)
     
     IMPLICIT NONE

     INTEGER :: j
     REAL(dp) :: vel(:,:), vcm_tmp(SIZE(vel,1),SIZE(vel,2)), vcm 


     DO j = 1, n_atms
       vcm_tmp(j,:) = vel(j,:) * atm_masses(j)
     END DO
     vcm = SUM(vcm_tmp) / SUM(atm_masses)
                           !print *,vcm
     DO j = 1, n_atms
       vel(j,:) = vel(j,:) - vcm
     END DO
   END SUBROUTINE do_fix_com
   
   !##########################################################
   SUBROUTINE do_scale_vel(vel)
   
     IMPLICIT NONE
     
     INTEGER :: k, j
     real(dp) :: vel(:,:), trmpT, tmp_scale, vel_scale, tempT
     
       k = 1 
       CALL do_calcT(tempT, vel)
       tmp_scale = ( (tmp - tempT) / tmp )  * 100
                           !PRINT *, 'tempT = k:',k ,'=', tempT, "T = ", tmp, "tmp_scale= ", tmp_scale
       DO WHILE ( (ABS(tmp_scale)>0.001) .and. (k<1000) )
          vel_scale = sqrt(tmp / tempT)
                           !Print *, 'vel_scale', k, '=', vel_scale
          DO j = 1 , n_atms
            vel(j,:) = vel(j,:) * vel_scale
          END DO 
          CALL do_calcT(tempT, vel)
          tmp_scale = ( (tmp - tempT) / tmp )  * 100
                            !PRINT *, 'tempT = k:',k ,'=', tempT, "T = ", tmp, "tmp_scale= ", tmp_scale
          k= k + 1  
       end DO

   END SUBROUTINE do_scale_vel                   

   !##########################################################
    
end module md

