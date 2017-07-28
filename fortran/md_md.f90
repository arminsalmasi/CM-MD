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
      DO i = 0,N_tstps 
        ALLOCATE(tmstp(i)%xyz(N_atms,3)) ! N*3 matrix
        ALLOCATE(tmstp(i)%vel(N_atms,3)) ! N*3 matrix
        ALLOCATE(tmstp(i)%acc(N_atms,3)) ! N*3 matrix
        ALLOCATE(tmstp(i)%frc(N_atms,3)) ! N*3 matrix

        tmstp(i)%xyz(:,:) = 0 
        tmstp(i)%vel(:,:) = 0       
        tmstp(i)%acc(:,:) = 0
        tmstp(i)%frc(:,:) = 0
 
      END DO
      
    end subroutine do_allocate_dtyp
   
   !##########################################################
   
   subroutine do_loop_tstps(tmstp)
  
   ! loops over all timesteps(N_tstp=t/dt)
   !   in tstp=0:
   !   in 0<tstp<N_tstp
   
   
     TYPE(tmstp_holder), DIMENSION(:), ALLOCATABLE :: tmstp
     INTEGER :: i ,j , k
     !REAL(dp) :: r(3) , vstd, tempT, vcm , vcm_tmp(n_atms,3), vel_scale, tmp_scale
 
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
         CALL do_fix_com(tmstp(i)%vel)
        
       ! Scale velocities with initial temperature
           CALL do_scale_vel(tmstp(i)%vel)
                                !print *,tmstp(i)%vel(1,:)
       ! ToDo: question? should values travers to i+1?       
   
       ELSE
         ! TODO : velocityverlet?
         
         CALL do_velverlet( tmstp(i)%xyz, tmstp(i)%vel, tmstp(i)%acc, tmstp(i)%frc, atm_masses, dt )  
         
         !tmstp(i)%frc(1,:) =1+1.5* tmstp(i)%frc(1,:)
          
         !print *, tmstp(i)%vel(1,:)


         ! write tstp data to a file 
       END IF
                               print *, tmstp(i)%vel(1,:)
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
   SUBROUTINE do_velVerlet(xyz, vel, acc, frc, masses, dt);
       
     IMPLICIT NONE              
     
     INTEGER :: i
     REAL(dp) :: xyz(:,:), vel(:,:), acc(:,:), frc(:,:), masses(:), dt 
     
     DO i=1 , 3 
       xyz(:,i) = xyz(:,i) + (vel(:,i) * dt) + (0.5 * acc(:,i) * dt**2);
     END DO   
  
     DO i = 1 , N_atms ! = size(atm_masses)
        vel(i,:) = vel(i,:) + 0.5 * ( frc(i,:) / masses(i) + acc(i,:)) * dt;
     END DO

     DO i = 1 , N_atms 
        acc(i,:) =  frc(i,:) / masses(i);
     END DO

   END SUBROUTINE do_velVerlet 
    
end module md

