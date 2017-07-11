module datastructure
  
  integer, parameter ::                         &
  sp = kind(1.0),                               &
  dp = selected_real_kind(2*precision(1.0_sp)), &
  qp = selected_real_kind(2*precision(1.0_dp))
    

  real(dp) :: t,   & ! time step lenght [real-dp]
              dt,  & ! simulation time [real-dp]
              Tmp    ! temperature [real-dp] 
  
  integer :: pot_type ! potentiel type [int] : 
                      !   1-Lennard-Jones
                      !   2-Embeded Atomic Method

  integer :: N ! total number of atoms in the system [int]

  integer, dimension(:),allocatable :: N_vec   ! list of numbers of atoms of each species [int(numberOfSpecies)]  
                                               ! ordered by atomic number / ascending
  
  integer, dimension(:),allocatable :: atm_num ! list of atomic number of species [int(numberOfSpecies)] 
                                               ! ordered / ascending
  
  real(dp), dimension(:), allocatable :: mass  ! atomic mass of species [real-dp(numberOfSpecies)]  
                                               ! ordered by atomic number / ascending   

  type tmstp_holder ! type holds all data of a timestep 
    real(dp), dimension(:,:,:), allocatable :: xyz ! coordinate [real-dp(N*3)]
    real(dp), dimension(:,:,:), allocatable :: vel ! velocitie [real-dp(N*3)]
    real(dp), dimension(:,:,:), allocatable :: acc ! acceleration [real-dp(N*3)]
  end type tmstp_holder

  type(tmstp_holder), dimension(:), allocatable :: tmstp ! list of all time steps [0 to t/dt]


end module datastructure
