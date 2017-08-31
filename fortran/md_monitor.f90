module monitor

USE global
!###############################################################################################
  contains
!###############################################################################################
    subroutine get_input(t, dt, temp, V, nTstps, sampOff, sampInt, nSamps, nP, nSp, nAt, atNums, atMas)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Declerations
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT NONE     
  
      REAL(dp), INTENT(OUT):: t        &! simulation time
                            , dt       &! timestep 
                            , temp     &! tempereture 
                            , V        &! Volume of the simulation box
  	                        , sampInt   ! time between saving data (purpes = speed up)
			
  
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: atmas(:)  ! atomic mass of species 
      	
       
      INTEGER, INTENT(OUT) :: nP	    &! total number of atoms in the system
                            , nSp   	&! number of species in the system
						    , sampOff   &! timestep offset between saving data 
							, nTstps    &! total number of timesteps
							, nSamps     ! number of samples
							

         
      INTEGER, ALLOCATABLE, INTENT(OUT) :: atNums(:)  &! contains atomic numbers of species A1, A2, ... 
                                         , nAt(:)	   ! contains numbers of atoms of each species n1, n2 ,..
	

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Subroutine MAIN
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      
	  t = 1e-10 
      dt = 1e-15
      temp = 273
      V = 1e-15 
	  nTstps = floor(t/dt)
      sampInt = 1e-14 
	  sampOff = floor(sampInt/dt)
	  nSamps = floor(sampInt/(dt * sampOff)) ! +2 for first and last timestep??
	  
      nSp = 2 
      allocate(nAt(nSp))
      nAt = (/5, 5/)  
	  allocate(atNums(nSp)) 
	  atNums = (/27, 28/) ! 27=Co, 28=Ni
	  nP = sum(nAt)
	  allocate(atMas(nP))
      7.9485017e-23
	  atMas(1:nAt(1))  = 9.786093775e-23 ! Cobalt
      atMas(nAt(1)+1: nP) = 9.746275017e-23 ! gr - Ni
	  
    end subroutine get_input
!###############################################################################################

end module monitor


