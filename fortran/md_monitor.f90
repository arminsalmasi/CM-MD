module monitor

  use datastructure
  
  contains
  
   subroutine get_input(t, dt, tmp, pot_typ, atm_nums, nums_atms, atm_masses, box_vol)
     
     integer :: pot_typ
  
     real(dp):: t    &
              , dt   &
              , Tmp  &
              , box_vol       
 
     integer, allocatable :: atm_nums(:) &  ! list of atomic numbers of all species
                           , nums_atms(:)   ! list of numbers of atoms of each species

     real(dp), allocatable :: atm_masses(:) ! list of atmic mass of each especies  

     N_spcs = 2 ! number of species \must be read from the monitor                           
     
     allocate(atm_nums(N_spcs)) 
     allocate(nums_atms(N_spcs))
     allocate(atm_masses(N_spcs))

     
     t = 1e-10
     
     dt = 1e-15

     Tmp = 0.1
     
     pot_typ = 1

     nums_atms = (/5, 5/)
     
     atm_nums = (/22, 27/)

     atm_masses = (/7.9485017e-23, 9.7860864e-23/)
    
   end subroutine get_input


end module monitor
