module datastructure
  
  Public ! if the module is used: all of following vartiable are accecebale 
   
  ! constants for real size defginitions
  integer, parameter ::                         &
  sp = kind(1.0),                               &
  dp = selected_real_kind(2*precision(1.0_sp)), &
  qp = selected_real_kind(2*precision(1.0_dp))
    
  real(dp) :: t   &    ! time step lenght [real-dp]
            , dt  &    ! simulation time [real-dp]
            , Tmp &    ! temperature [real-dp] 
            ,  box_vol !volume of the simulation box
   
  integer :: pot_typ ! potentiel type [int] : 
                      !   1-Lennard-Jones
                      !   2-Embeded Atomic Method
  
  integer :: N ! total number of atoms in the system [int]
  
  integer,allocatable :: nums_atms(:)   ! list of numbers of atoms of each species [int(numberOfSpecies)]  
                                        ! ordered by atomic number / ascending
  
  integer, allocatable :: atm_nums(:)   ! list of atomic number of species [int(numberOfSpecies)] 
                                        ! ordered / ascending
  
  real(dp), allocatable :: atm_masses(:)  ! atomic mass of species [real-dp(numberOfSpecies)]  
                                          ! ordered by atomic number / ascending   
  
  ! type holds all data of a timestep 
  type tmstp_holder 
    real(dp), allocatable :: xyz(:,:) ! coordinate [real-dp(N*3)]
    real(dp), allocatable :: vel(:,:) ! velocitie [real-dp(N*3)]
    real(dp), allocatable :: acc(:,:) ! acceleration [real-dp(N*3)]
  end type tmstp_holder
  
  type(tmstp_holder), dimension(:), allocatable :: tmstp ! list of all time steps [0 to t/dt]
  
  
end module datastructure
