MODULE io
!###############################################################################################
USE glbl
!###############################################################################################
CONTAINS
!###############################################################################################
!###############################################################################################
SUBROUTINE get_input(t, dt, Temp, V, nTstps, sampInt, nP, nSp, nAt, atNums, atNumsVec, atMas)
!###############################################################################################
  IMPLICIT NONE     
  
  REAL(dp):: t       &! simulation time
           , dt      &! timestep 
           , Temp    &! Tempereture 
           , sampInt &! time between saving data (purpes = speed up)
           , V        ! Volume of the simulation box

  REAL(dp), ALLOCATABLE :: atMas(:)  ! atomic mass of species 
   
  INTEGER :: nP      &! total number of atoms in the system
           , nSp     &! number of species in the system
           , nTstps  &! total number of timesteps
           , nCo     &
           , nNi     
     
  INTEGER, ALLOCATABLE :: atNums(:)    &! contains atomic numbers of species A1, A2, ... 
                        , atNumsVec(:) &! A1,A1,A1,....; A2,A2,A2... 
                        , nAt(:)        ! contains numbers of atoms of each species n1, n2 ,..

  WRITE(*,*) " Program MD  "
  WRITE(*,*) " a fortran code for a simple Molocular Dynamics simulation of Ni, Co alloy using          "
  WRITE(*,*) " EAM potentials from NIST EAM repository:                                                 " 
  WRITE(*,*) " https://www.ctcms.nist.gov/potentials/Co-Ni.html                                         "
  WRITE(*,*) " Potential data is tabulated by Y. Mishin (George Mason Univ.) on 17 Sept. 2013           "
  WRITE(*,*) " and was posted on 17 Jan. 2014. This version is compatible with LAMMPS and The reference "
  WRITE(*,*) " information was updated on 26 Aug. 2015                                                  "
  WRITE(*,*) " Format: EAM/alloy setfl                                                                  "
  WRITE(*,*) " https://www.ctcms.nist.gov/potentials/Download/Ni-Al-Co-YM13/Mishin-Ni-Co-2013.eam.alloy "
  WRITE(*,*) " This code is implemented by Armin Salmasi at KTH, Stockholm,Sweden in 2017               "
  WRITE(*,*) " Feel free to use any part of this code.                                                  "  
  WRITE(*,*) "******************************************************************************************"
  t  = 1e-10         
  dt = 1e-15         
  Temp =273.15       
  sampInt =1e-14     
  V = 109.595564966  
  nCo = 5
  nNi = 5

  WRITE(*,*) ""
  WRITE(*,*) 'enter number of Co atoms/5/:'
  READ(*,*) , nCo ! make sure it is integer
  WRITE(*,*) ""
  WRITE(*,*) 'enter number of Ni atoms/5/:'
  READ(*,*) , nNi ! make sure it is integer
  WRITE(*,*) ""
  WRITE(*,*) 'enter volume of the cell in Angestrom/109.595564966/:'
  READ(*,*) , V  != (6.6e-6 / 6.0221415e23) * 10 * 1e30 !m3 / avogadro# * #atoms * #toAngstrom^3 
  V = V * 1D0 ! make sure it is dp float
  WRITE(*,*) ""
  WRITE(*,*) 'enter simulation time in seconds/1e-10/:'
  READ(*,*) , t ! make sure it is dp float
  t= t *1D0
  WRITE(*,*) ""
  WRITE(*,*) 'enter timestep lenght in seconds/1e-15/:'
  READ(*,*) , dt 
  t= t *1d0 ! make sure it is dp float
  WRITE(*,*) ""
  WRITE(*,*) 'enter interval between saves in seconds/1e-13/:'
  READ(*,*) , sampInt
  sampInt = sampInt * 1D0
  WRITE(*,*) ""

  nTstps = floor(t/dt)
  WRITE(*,*) 'number of timesteps is: ', nTstps
  
  nSp = 2 
  
  ALLOCATE(nAt(nSp))
  nAt = (/nCo, nNi/)  
  np = SUM(nAt)
      
  ALLOCATE(atNums(nSp)) 
  ALLOCATE(atNumsVec(np)) 
  atNums = (/27, 28/)                   ! 27=Co, 28=Ni
  atNumsVec(1:nAt(1))     = atNums(1)
  atNumsVec(nAt(1)+1: nP) = atNums(2)
      
  ALLOCATE(atMas(nP))
  atMas(1:nAt(1))     = 9.786093775e-23 ! Co
  atMas(nAt(1)+1: nP) = 9.746275017e-23 ! gr - Ni
  
END SUBROUTINE get_input
!###############################################################################################
!###############################################################################################
SUBROUTINE do_get_EAMPotData(EAMdata)
!###############################################################################################
   IMPLICIT NONE
   INTEGER :: nrho       &
            , nr         &
            , ix, jx, kx
   REAL(dp), ALLOCATABLE :: allData(:,:)    &
                          , allData_line(:) 
   TYPE(EAM_data) :: EAMdata 
   
   EAMdata%nrho = 5000 
   EAMdata%drho = 1.2999078e-03 
   EAMdata%nr = 5000 
   EAMdata%dr = 1.2999078e-03 
   EAMdata%cutoff = 6.499539

   ALLOCATE(allData(1:7000, 1:5))
   ALLOCATE(allData_line(1: SIZE(allDATA,1)* SIZE(allDATA,2)))

   ALLOCATE(EAMdata%F1(nrho))
   ALLOCATE(EAMdata%rho1(nr))
   ALLOCATE(EAMdata%F2(nrho))
   ALLOCATE(EAMdata%rho2(nr))
   ALLOCATE(EAMdata%phi11(nr))
   ALLOCATE(EAMdata%phi21(nr))
   ALLOCATE(EAMdata%phi22(nr))
   OPEN(UNIT=22,FILE="pot-Ni-Co-old.dat",FORM="FORMATTED",STATUS="OLD",ACTION="READ")
     DO ix = 1 , 7000
       READ(22,*) allData(ix,:)
     END DO
   CLOSE(UNIT=22)
   kx=1
   DO ix = 1 , 7000
     DO jx = 1 , 5
       allData_line(kx) = allData(ix,jx)
       kx=kx+1
     END DO
   END DO
   EAMdata%F1    = allData_line(1:5000)
   EAMdata%rho1  = allData_line(5001:10000)
   EAMdata%F2    = allData_line(10001:15000)
   EAMdata%rho2  = allData_line(15001:20000)
   EAMdata%phi11 = allData_line(20001:25000)
   EAMdata%phi21 = allData_line(25001:30000)
   EAMdata%phi22 = allData_line(30001:35000)

END SUBROUTINE do_get_EAMPotData

!###############################################################################################
END MODULE io


