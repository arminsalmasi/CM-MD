module datastructure
  
  Public 
   
  ! constants for real size defginitions
  integer, parameter :: sp = kind(1.0),&
                        dp = selected_real_kind(2*precision(1.0_sp)), &
                        qp = selected_real_kind(2*precision(1.0_dp))
    
  real(dp) :: t      &! time step lenght [real-dp]
            , dt     &! simulation time [real-dp]
            , Tmp    &! temperature [real-dp] 
            , box_vol !volume of the simulation box
   
  integer :: pot_typ ! potentiel type [int] : 
                       ! 1-Lennard-Jones
                       ! 2-Embeded Atomic Method
  
  integer :: N_atms  &! total number of atoms in the system [int]
           , N_tstps &! number of timesteps = t/dt
           , N_spcs   ! number of species in the system

  integer,allocatable :: nums_atms(:) 
  ! list of numbers of atoms of each species [int(numberOfSpecies)]  
    ! ordered by atomic number / ascending
  
  integer, allocatable :: atm_nums(:) 
  ! list of atomic number of species [int(numberOfSpecies)] 
    ! ordered / ascending
  
  real(dp), allocatable :: atm_masses(:)  
  ! atomic mass of species [real-dp(numberOfSpecies)]  
    ! ordered by atomic number / ascending   
  
  
  ! type tmstp_holder holds all data of a single timestep 
  type tmstp_holder 
    real(dp), allocatable :: xyz(:,:) ! coordinate [real-dp(N*3)]
    real(dp), allocatable :: vel(:,:) ! velocitie [real-dp(N*3)]
    real(dp), allocatable :: acc(:,:) ! acceleration [real-dp(N*3)]
  end type tmstp_holder
  
  
  ! vector of tmstp_hlolders
  !type(tmstp_holder), dimension(0:), allocatable :: tmstp 
  
end module datastructure
