module utilities

  use datastructure

  contains  
   !##########################################################
   subroutine do_allocate_dtyp(t,dt)
   ! input t: simulation time 
   ! input dt: prescribed lenght of timesteps
   ! allocates vector of datatype by t/dt=# of timesteps
   
      integer :: N_tstps & ! number of timesteps
               , i         ! counter variable
      
      real(dp):: t  & ! simulation time
               , dt   ! lenght of timesteps   
      
      ! allocate the tmstp vector by N_tstp = t/dt
      N_tstps = floor(t/dt)
            
      allocate(tmstp(N_tstps))
      
      ! loop over all timesteps to allocate variables in 
      ! the tmstp elements(datatpe)
      do i = 1, size(tmstp) 
        allocate(tmstp(i)%xyz(N_atms,3)) ! N*3 matrix
        allocate(tmstp(i)%vel(N_atms,3)) ! N*3 matrix
        allocate(tmstp(i)%acc(N_atms,3)) ! N*3 matrix
      end do
      
   end subroutine do_allocate_dtyp
   
   !##########################################################
   
   subroutine do_loop_tstps(t, dt)
   ! input t: simulation time 
   ! input dt: prescribed lenght of timesteps
   ! loops over all timesteps(N_tstp=t/dt)
   !   in tstp=0:
   !   in 0<tstp<N_tstp
	
	  integer N_tstp

    real(dp):: t  & ! simulation time
             , dt   ! lenght of timesteps

	  N_tstp = floor(t/dt)
    
    !loop over all timesteps
	  do i= 0 , N_tstp
	    if (i == 0) then
	      write(*,*) i
		  ! TODO : intialize tstp 0 
	    else
	      write(*,*) i
		  ! TODO : do_md
	    end if
	  end do
	  
	end subroutine do_loop_tstps
  
  !##########################################################

end module utilities
