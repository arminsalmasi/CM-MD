module monitor

  use datastructure
  
  contains
  !##########################################################
  subroutine get_input(t, dt, tmp, pot_typ, atm_nums, nums_atms, atm_masses, box_vol)
  
  ! output: t: simulation time [s]
  !         dt: lenght of each timescale [s]
  !         tmp: temperature [k]
  !         pot_typ: type of potential to use:
  !                  0-Lennard-Jones or &&
  !                  1-Embedded Atom Method (EAM) potentials 
  !         atm_nums: list of atomic numbers of all species
  !         nums_atms: list of numbers of atoms of each species
  !         atm_masses: list of atmic mass of each especies 
  !         box_vol: volume of the cubic simulation domain [m^3]
  ! intent: reads simulation parameters from the keyboard
  
  ! ToDo: implement read from Input
  
  
     
     integer :: pot_typ
  
     real(dp):: t    &  
              , dt   &   
              , Tmp  &   
              , box_vol       
 
     integer, allocatable :: atm_nums(:) &  
                           , nums_atms(:)    

     real(dp), allocatable :: atm_masses(:)   

     ! ToDo: must be read from the monitor                           
     N_spcs = 2 ! number of species 
     
     
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
   !##########################################################

end module monitor
