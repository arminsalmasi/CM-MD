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
  use utilities 
  
  write(*,*) 'This is a simple fortran MD code'

! Read from input: 
  ! t, dt, tmp, pot_typ, atm_nums, nums_atms, atm_masses, box_vol
  call get_input()

! allocate tmstp datatype vector
! allocate all vector fileds of all cells of tmstp vector
  call do_allocate_dtyp()

! loops overtime steps
  call do_loop_tstps()

    
end program MD
