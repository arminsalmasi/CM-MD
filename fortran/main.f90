! Program MD
! MD is a fortran code for a simple MD simulation.
! It is possible to use eather Lennard-Jones or &&
! Embedded Atom Method (EAM) potentials.
!
! case study is a monoatomic metal in this case Ti
! 
! The code is developed by Armin Salmasi at KTH, Stockholm, Sweden.
! First code commit (commit#2-2): 2017-07-11
! 
! Feel free to use any part of tyhe code.


Program MD
  
  use datastructure
  use monitor 
  
  write(*,*) 'This is a simple fortran MD code'
  call get_input(t, dt, temp, pot_typ, nums_atms, &
                         N, atm_nums, atm_masses)

end program MD
