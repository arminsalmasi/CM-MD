module monitor

USE global
!###############################################################################################
  contains
!###############################################################################################
    subroutine get_input(t, dt, temp, V, nTstps, sampInt, nP, nSp, nAt, atNums, atNumsVec, atMas)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Declerations
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT NONE     
  
      REAL(dp):: t        &! simulation time
               , dt       &! timestep 
               , temp     &! tempereture 
               , V        &! Volume of the simulation box
               , sampInt   ! time between saving data (purpes = speed up)
  
      REAL(dp), ALLOCATABLE :: atmas(:)  ! atomic mass of species 
       
      INTEGER :: nP        &! total number of atoms in the system
               , nSp       &! number of species in the system
               , nTstps     ! total number of timesteps
         
      INTEGER, ALLOCATABLE :: atNums(:)   &! contains atomic numbers of species A1, A2, ... 
                            , atNumsVec(:) &! A1,A1,A1,....; A2,A2,A2... 
							, nAt(:)       ! contains numbers of atoms of each species n1, n2 ,..
							


      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Subroutine MAIN
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      t = 1e-10 
      dt = 1e-15
      temp = 273
      V = (6.6e-6 / 6.0221415e23) * 10 * 1e30 !m3 / avogadro# * #atoms * #toAngstrom^3 
!PRINT*, 'volume of 10 atoms[Ang^3] = ', V, 'parameter[Ang] = ', V**(1.0/3.0)
      nTstps = floor(t/dt)
      sampInt = 1e-14 
  
      nSp = 2 
      
	  ALLOCATE(nAt(nSp))
      nAt = (/5, 5/)  
	  np = SUM(nAt)
      
	  ALLOCATE(atNums(nSp)) 
	  ALLOCATE(atNumsVec(np)) 
      atNums = (/27, 28/) ! 27=Co, 28=Ni
      atNumsVec(1:nAt(1)) = atNums(1)
	  atNumsVec(nAt(1)+1: nP) = atNums(2)
      
	  ALLOCATE(atMas(nP))
      atMas(1:nAt(1))  = 9.786093775e-23 ! Cobalt
      atMas(nAt(1)+1: nP) = 9.746275017e-23 ! gr - Ni
  
    end subroutine get_input
!###############################################################################################

end module monitor


