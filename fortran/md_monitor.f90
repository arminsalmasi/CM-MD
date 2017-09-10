!###############################################################################################
!###############################################################################################
! Module io - program md
! Input output method of md NiCo EAM fortran code
! Armin Salmasi, 2017-09-09
! Contains:
!   Subroutine get_input which reads input parameters of md simulation
!   Subroutine get_EAMPotData which reads parameters of EAM potential data from a cleaned file
!###############################################################################################
!###############################################################################################
MODULE io
!###############################################################################################
USE glbl ! Module global contains holestic public variables an data types
!###############################################################################################
CONTAINS
!###############################################################################################
!###############################################################################################
SUBROUTINE get_input(t, dt, Temp, V, nTstps, sampInt, nP, nSp, nAt, atNums, atNumsVec, atMas)
!###############################################################################################
! Read input data from monitor: 
! simulation time, timestep, interval between saving samples: seconds
! temperature: kelvin
! number of atoms, number of atomes of each typs (nCo and nNi) 
! Volume : in Angstrom^3
!###############################################################################################
  IMPLICIT NONE     
  
  REAL(dp):: t       &! simulation time [sec]
           , dt      &! timestep [sec]
           , Temp    &! Tempereture [kelvin]
           , sampInt &! time between saving data (purpes = speed up) [sec]
           , V        ! Volume of the simulation box [angstrom^3]

  REAL(dp), ALLOCATABLE :: atMas(:)  ! atomic mass of species [gram]
   
  INTEGER :: nP      &! total number of atoms in the system
           , nSp     &! number of species in the system
           , nTstps  &! total number of timesteps
           , nCo     &! number of Co atoms
           , nNi      ! number of Ni atoms
     
  INTEGER, ALLOCATABLE :: atNums(:)    &! contains atomic numbers of species A1, A2, ... 
                        , atNumsVec(:) &! A1,A1,A1,....; A2,A2,A2... 
                        , nAt(:)        ! contains numbers of atoms of each species n1, n2 ,..

  WRITE(*,*) "******************************************************************************************"
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
  WRITE(*,*) 'enter simulation temperature in kelvin/273.15/:'
  READ(*,*) , Temp
  Temp = Temp * 1D0
  WRITE(*,*) ""

  WRITE(*,*) ""
  WRITE(*,*) 'enter number of Co atoms/5/:'
  READ(*,*) , nCo ! make sure it is integer
  
  WRITE(*,*) ""
  WRITE(*,*) 'enter number of Ni atoms/5/:'
  READ(*,*) , nNi ! make sure it is integer
  
  WRITE(*,*) ""
  WRITE(*,*) 'Use ThermoCalc to calulate molar volume of the composition in m^3'
  WRITE(*,*) 'molar volume * (total number of atoms * 10^30 / avogadro number)'
  WRITE(*,*) 'enter volume of the cell in Ang^3 /109.595564966/:'
  READ(*,*) , V  != (6.6e-6 / 6.0221415e23) * 10 * 1e30 !m3 / avogadro# * #atoms * #toAngstrom^3 
  V = V * 1D0 ! make sure it is dp float
  
  WRITE(*,*) ""
  WRITE(*,*) 'enter simulation time in seconds/1e-10/:'
  READ(*,*) , t ! make sure it is dp float
  t= t *1D0
  
  WRITE(*,*) ""
  WRITE(*,*) 'enter timestep lenght in seconds/1e-15/:'
  READ(*,*) , dt 
  dt= dt *1d0 ! make sure it is dp float
  
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
! Read EAM potential data of Co-Ni system from pot-Ni-Co-old.dat file
! Note that headers of the original .Alloy file (From NIST data base) is removed from this file
! File contains:
! 5000 elements F of Ni,
! 5000 elemnts Rho of Ni,
! 5000 elements F of Co, 
! 5000 elements Rho of Co,
! 5000 elements phi if Ni-Ni
! 5000 elements phi of Ni-co
! 5000 elements phi of Co-Co
!###############################################################################################
   IMPLICIT NONE
   INTEGER :: nrho       &
            , nr         &
            , ix, jx, kx
   
   REAL(dp), ALLOCATABLE :: allData(:,:)    &
                          , allData_line(:) 
   
   TYPE(EAM_data) :: EAMdata 
   
   EAMdata%nrho   = 5000 
   EAMdata%drho   = 1.2999078e-03 
   EAMdata%nr     = 5000 
   EAMdata%dr     = 1.2999078e-03 
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

! read all data from the file to a 35000*5 matrix 
   OPEN(UNIT=22,FILE="pot-Ni-Co-old.dat",FORM="FORMATTED",STATUS="OLD",ACTION="READ")
     DO ix = 1 , 7000
       READ(22,*) allData(ix,:)
     END DO
   CLOSE(UNIT=22)

! make the 35000*5 matrix llinear
   kx=1
   DO ix = 1 , 7000
     DO jx = 1 , 5
       allData_line(kx) = allData(ix,jx)
       kx=kx+1
     END DO
   END DO

! put data in variables to pass to the main   
   EAMdata%F1    = allData_line(1:5000)
   EAMdata%rho1  = allData_line(5001:10000)
   EAMdata%F2    = allData_line(10001:15000)
   EAMdata%rho2  = allData_line(15001:20000)
   EAMdata%phi11 = allData_line(20001:25000)
   EAMdata%phi21 = allData_line(25001:30000)
   EAMdata%phi22 = allData_line(30001:35000)

END SUBROUTINE do_get_EAMPotData

!###############################################################################################
!###############################################################################################
!###############################################################################################
!###############################################################################################
SUBROUTINE put_SampsToFiles(samp, nPart, nSpcis, atNums, nAtms, Temp, vol, t, dt, nTstps &
                            , smpInt, sampK)
!###############################################################################################
! Write simultion results from timestep samples to files:
! dat.md: # of atoms, # of species, # Atoms of each species, # of samples   
! xyz.md: line#, sample#, atom#, sample time, x, y, z    of (atom#)
! vel.md: line#, sample#, atom#, sample time, Vx, Vy, Vz of (atom#)   
! acc.md: line#, sample#, atom#, sample time, ax, ay, az of (atom#)   
! frc.md: line#, sample#, atom#, sample time, fx, fy, fz of (atom#)   
! pot.md: line#, sample#, atom#, sample time, EAM        of (atom#)   
! erg.md: line#, sample#, sample time, total energy, kinetic energy, potential energy, error  
! grl.md: general information
!###############################################################################################
  IMPLICIT NONE

  INTEGER :: jx, ix, l, k , nPart, nSpcis, nTstps, sampK, atNums(:), nAtms(:)

  REAL(dp) :: Temp, vol, t, dt, smpInt
  
  TYPE(timestep_sample) :: samp(:)
  open(unit=27, file="dat.md", action="write", status="replace")
  open(unit=26, file="grl.md", action="write", status="replace")
  open(unit=25, file="erg.md", action="write", status="replace")
  open(unit=24, file="pot.md", action="write", status="replace")
  open(unit=23, file="frc.md", action="write", status="replace")
  open(unit=22, file="acc.md", action="write", status="replace")
  open(unit=21, file="vel.md", action="write", status="replace")
  open(unit=20, file="xyz.md", action="write", status="replace")
  
  write(26,*), '******************************************************************************'
  write(26,*), 'All values are sorted by atomic numbers of species (ascending)'
  write(26,*), '******************************************************************************'
  write(26,*), 'Number of atoms = ', nPart
  write(26,*), 'Number of species = ', nSpcis
  write(26,*), 'Atomic numbers of species = ', atNums
  write(26,*), 'Numbers of atoms of each species = ', nAtms
  write(26,*), '******************************************************************************'
  write(26,*), 'Simulation Temp [K] = ', Temp
  write(26,*), 'vol of the box[Ang^3] = ' , vol
  write(26,*), 'Simulation time [sec]= ', t
  write(26,*), 'Time step [sec] = ', dt
  write(26,*), 'Number of time steps = ', nTstps
  write(26,*), '******************************************************************************'
  write(26,*), 'Sampling interval [sec] = ', smpInt
  write(26,*), 'Number of Samples (starting from 0) / per atom = ', sampK
  write(26,*), '******************************************************************************'
  write(26,*), 'File Format :: xyz.md, vel.md. acc.md, frc.md'
  write(26,*), "line#	 ", "sample# ","time[sec]	 ",&
               "x [Ang]	 " ,"y [Ang]	 ","z [Ang] "
  write(26,*), '******************************************************************************'
  write(26,*), 'File Format :: pot.md - potential energy of each atom'
  write(26,*), "line#	 ", "sample#	 ","time[sec]	 ", "Pot [ev]	 "
  write(26,*), '******************************************************************************'
  write(26,*), 'File Format :: Energies.md - Energy of the system'
  write(26,*), "line#	 ", "sample#	 ","time[sec]	 ", "Total Energy [ev]	 ", &
                "kinetic Energy [ev]	 ","Potential Energy [ev]	 ", "Relative erroror"
  write(26,*), '******************************************************************************'
  
  
  k = 1
  DO jx = 1, nPart
    DO ix = 1, sampK
      WRITE(20,*), k, ix, jx, samp(ix)%time, samp(ix)%xyz(jx,1:3)
      WRITE(21,*), k, ix, jx, samp(ix)%time, samp(ix)%vel(jx,1:3)
      WRITE(22,*), k, ix, jx, samp(ix)%time, samp(ix)%acc(jx,1:3)
      WRITE(23,*), k, ix, jx, samp(ix)%time, samp(ix)%frc(jx,1:3)
      WRITE(24,*), k, ix, jx, samp(ix)%time, samp(ix)%eam(jx)
      k= k+1
    END DO
  END DO

  k = 1
  DO ix = 1, sampK
    WRITE(25,*), k, ix, samp(ix)%time,  samp(ix)%E, samp(ix)%Kin, samp(ix)%E - samp(ix)%Kin  , samp(ix)%error
    k = k+1                   
  END DO
  
  write(27,*), nPart, nSpcis, nAtms, sampk
  
  CLOSE(UNIT=20)
  CLOSE(UNIT=21)
  CLOSE(UNIT=22)
  CLOSE(UNIT=23)
  CLOSE(UNIT=24)
  CLOSE(UNIT=25)
  CLOSE(UNIT=26)
  CLOSE(UNIT=27)
  
END SUBROUTINE put_SampsToFiles  

!###############################################################################################
!###############################################################################################
END MODULE io
!###############################################################################################

