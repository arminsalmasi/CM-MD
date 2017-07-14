module utilities

  use datastructure

  contains  

   subroutine allocate_dtyp(t,dt)
      
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
      
    end subroutine allocate_dtyp

end module utilities
