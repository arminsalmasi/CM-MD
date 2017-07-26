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
     INTEGER :: i ,j 
     REAL(dp) :: r(3) , vstd, tempT
 
     DO i= 0 , N_tstps
       IF (i == 0) then
         ! TODO : intialize tstp 0 
         ! call do_rand_xyz()
         DO j = 1, N_atms
           
           ! randomize xyz of atoms in cell 0 
           call do_genRand(r,3,j)
           r(:)=0
           tmstp(i)%xyz(j,:) = r(:) * ( box_vol**(1.0/3.0))
           
           ! randomize velocities of atoms in cell zero
             !(:) = random_normal(5) !Stanford algorithm
             !print *, 'first algorithm', r
             r(:) = 0             
             r(:) = rand_normal2(3) ! open source from roseta code
                          
             vstd = sqrt( tmp * kb / atm_masses(j) )
             tmstp(i)%vel(j,:) = 0 + r(:) * vstd  
             !print *, tmstp(i)%vel(j,:) , 'vstd', vstd
            
           
      
         END DO
             print *,'velocity', tmstp(i)%vel 
             CALL do_calcT(tempT, tmstp(i)%vel)
             PRINT *, tempT, tmp
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

