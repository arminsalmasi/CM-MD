module monitor

  use datastructure
  
  contains
!##########################################################
    subroutine get_input()
     
    ! intent: reads simulation parameters from the keyboard
      ! N_spcs: number of species in the system 
      ! t: simulation time [s]
      ! dt: lenght of each timescale [s]
      ! tmp: temperature [k]
      ! pot_typ: type of potential to use:
        ! 0-Lennard-Jones or &&
        ! 1-Embedded Atom Method (EAM) potentials 
      ! atm_nums: list of atomic numbers of all species
      ! nums_atms: list of numbers of atoms of each species
      ! atm_masses: list of atmic mass of each especies 
      ! box_vol: volume of the cubic simulation domain [m^3]
      ! N_tstps: number of timesteps =t/dt
      ! N_atms : total number of atoms  inthe system
  
    ! ToDo: implement read from Input
        
      N_spcs = 2 ! number of species 
         
      allocate(atm_nums(N_spcs)) 
      allocate(nums_atms(N_spcs))
      
      t = 1e-10 
      
      dt = 1e-15
      
      tmp = 0.1
     
      pot_typ = 1

      atm_nums = (/22, 27/)
      
      nums_atms = (/5, 5/)
      
          
      N_atms = sum(nums_atms)
      
      allocate(atm_masses(N_atms))
      
      DO i = 1 , nums_atms(1)
         atm_masses(i)  = 7.9485017e-23
      END DO
      DO i = nums_atms(2)+1 , nums_atms(1)*2
        atm_masses(i) = 9.7860864e-23
      END DO
    
      box_vol = 1e-15   
   
      N_tstps = floor(t/dt)
      
       
    end subroutine get_input
!##########################################################

end module monitor
