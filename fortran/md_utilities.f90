module utilities

  use datastructure

  contains  
!##########################################################
   subroutine do_allocate_dtyp(t,dt)
      
      integer :: N_tstps &
               , i
      
      real(dp):: t       &
               , dt   
      
      N_tstps = floor(t/dt)
          
      allocate(tmstp(N_tstps))
      
      do i = 1, size(tmstp)
        allocate(tmstp(i)%xyz(N_atms,3))
        allocate(tmstp(i)%vel(N_atms,3))
        allocate(tmstp(i)%acc(N_atms,3))
      end do
      
    end subroutine do_allocate_dtyp
!##########################################################
    subroutine do_loop_tstps(t, dt)
	
	  integer N_tstp

      real(dp):: t       &
               , dt   

	  N_tstp = floor(t/dt)
	
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
