module md

use datastructure
use utilities


contains
!##########################################################
    
    subroutine do_allocate_dtyp(tmstp)
    
    ! allocate tmstp datatype vector
    ! allocate all vector fileds of all cells of tmstp vector
    
      type(tmstp_holder), dimension(:), allocatable :: tmstp
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
         DO j = 1, N_atms
          
           call do_genRand(r,3,j)
           tmstp(i)%xyz(j,:) = r(:) * ( box_vol**(1.0/3.0))
         
         END DO
                      
         ! randomize velocities of atoms in cell zero
         DO j = 1 , N_atms
           r(:) = rand_normal2(3) ! open source from roseta code
           vstd = sqrt( tmp * kb / atm_masses(j) )
           tmstp(i)%vel(j,:) = 0 + r(:) * vstd  
         END Do  
                               !CALL do_calcT(tempT, tmstp(i)%vel)
                               !PRINT *, 'tempT= ', tempT,'T= ', tmp
         
         ! fix center of mass
                               !print *,tmstp(i)%vel(1,:)
         DO j = 1, n_atms
           vcm_tmp(j,:) = tmstp(i)%vel(j,:) * atm_masses(j)
         END DO
         vcm = sum(vcm_tmp) / sum(atm_masses)
                               !print *,vcm
         DO j = 1, n_atms
           tmstp(i)%vel(j,:) = tmstp(i)%vel(j,:) - vcm
         END DO
                              !print *,tmstp(i)%vel(1,:)
        
         ! Scale velocities with initial temperature
           k = 1 
           CALL do_calcT(tempT, tmstp(i)%vel)
           tmp_scale = ( (tmp - tempT) / tmp )  * 100
                               !PRINT *, 'tempT = k:',k ,'=', tempT, "T = ", tmp, "tmp_scale= ", tmp_scale
           DO WHILE ( (ABS(tmp_scale)>0.001) .and. (k<1000) )
              
              vel_scale = sqrt(tmp / tempT)
                               !Print *, 'vel_scale', k, '=', vel_scale
              DO j = 1 , n_atms
                tmstp(i)%vel(j,:) = tmstp(i)%vel(j,:) * vel_scale
              END DO 
           
              CALL do_calcT(tempT, tmstp(i)%vel)
              tmp_scale = ( (tmp - tempT) / tmp )  * 100
              
                               !PRINT *, 'tempT = k:',k ,'=', tempT, "T = ", tmp, "tmp_scale= ", tmp_scale
              k= k + 1  
             
           end DO
                           
       ELSE
         ! TODO : do md on each timestep
         ! do_md()
!         print *, tmstp(i)%xyz(j,:)            
         ! save_tstp() 
       END IF
     END DO
     
   END SUBROUTINE do_loop_tstps
   
   !##########################################################

end module md

