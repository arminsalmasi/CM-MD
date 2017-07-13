module monitor

  use datastructure
  
  contains
  
   subroutine get_input(t, dt, temp, pot_typ, nums_atms, &
                         N, atm_nums, atm_masses)
     
     integer :: pot_typ &
              , N_atms  &
              , N_spcs  &
              , N_tstps &
              ,i
     
     real(dp):: t    &
          , dt   &
          , Tmp        
 
     integer, allocatable :: atm_nums(:) &
                           , nums_atms(:)  !atomicnumberslist 

     real(dp), allocatable :: atm_masses(:)



     N_spcs = 2                
     
     allocate(atm_nums(N_spcs))
     allocate(nums_atms(N_spcs))
     allocate(atm_masses(N_spcs))

     
     t = 1e-10
     
     dt = 1e-15

     Tmp = 0.1
     
     pot_typ = 1

     nums_atms = (/5, 5/)
     
     N_atms = sum(nums_atms)
     
     atm_nums = (/22, 27/)

     atm_masses = (/7.9485017e-23, 9.7860864e-23/)
     
     N_tstps = floor(t/dt)
     write(*,*) N_tstps 
     
     allocate(tmstp(N_tstps))
     
     
     do i = 1, size(tmstp)
       allocate(tmstp(i)%xyz(N_atms,3))
       allocate(tmstp(i)%vel(N_atms,3))
       allocate(tmstp(i)%acc(N_atms,3))
    end do

   end subroutine get_input


end module monitor
