! Program MD
! a fortran code for a simple Molocular Dynamics simulation of Ni, Co alloy using
! EAM potentials from NIST EAM repository: 
! https://www.ctcms.nist.gov/potentials/Co-Ni.html
! Potential data is tabulated by Y. Mishin (George Mason Univ.) on 17 Sept. 2013 
! and was posted on 17 Jan. 2014. This version is compatible with LAMMPS and The reference 
! information was updated on 26 Aug. 2015 
! Format: EAM/alloy setfl 
! https://www.ctcms.nist.gov/potentials/Download/Ni-Al-Co-YM13/Mishin-Ni-Co-2013.eam.alloy
! ! This code is implemented by Armin Salmasi at KTH, Stockholm,Sweden in 2017
! ! Feel free to use any part of this code.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM md
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Anounces
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  USE io
  USE glbl
  USE util
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declerations - main local variables
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  IMPLICIT NONE
  REAL(dp) :: t, dt, Temp, vol &! simulation, time timestep, vol
            , smpInt           &! time between saving data (purpes = speed up)
            , E, E0            &! Total energy ::  time dependent , initial
            , Kin               ! Kinetic energy

  REAL(dp), ALLOCATABLE :: AtMas(:) &! atomic mass of species 
                         , xyz(:,:) &! coordinate 
                         , vel(:,:) &! velocitie
                         , acc(:,:) &! acceleration
                         , frc(:,:) &! forces
                         , EAM(:)    ! EAM potentials      

  INTEGER :: nPart   &! total number of atoms in the system
           , nTstps  &! number of timesteps
           , nSpcis  &! number of species in the system
           , nSamp   &! number of samples
           , sampOfs &! number of timesteps between two samples except for the last one
           , samptp  &! tstp of sampling
           , sampK   &! sample counter  
           , tstp ,ix, jx , k, l ! counters

  INTEGER, ALLOCATABLE :: AtNums(:)    &! atomic numbers of each species A1, A2, ... 
                        , atNumsVec(:) &!  A1,A1,A1,....; A2,A2,A2... 
                        , nAtms(:)      ! numbers of atoms of each species n1, n2 ,..

  TYPE(timestep_sample), DIMENSION(:), ALLOCATABLE :: samp

  TYPE(EAM_data) :: EAMdata        
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! MAIN
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Read from input: all vectores atomic number ordered, ascending   
! Lenght unit: angestrom
  CALL get_input(t, dt, Temp, vol, nTstps, smpInt, &
                 nPart, nSpcis, nAtms, AtNums, atNumsVec, AtMas)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Allocate and initialize all to Zero
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  sampOfs = Floor(smpInt / dt * 1.0)
  nSamp = 1 + FLOOR(nTstps / sampOfs * 1.0)
  CALL do_Alloc(xyz, vel, frc, acc, EAM, samp, nSamp, nPart)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! read EAM potential data from Ni-Co eam.alloy file
! unites: Angestrom, ev
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  CALL do_get_EAMPotData(EAMdata)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!time step loop
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  samptp = 0
  sampK = 0
  DO tstp= 0, nTstps
    IF (tstp==0) THEN
      CALL do_rand_xyz(xyz, vol)
      CALL do_rand_vel(vel, Temp, AtMas)          
      CALL do_fix_centerOfmass(vel, AtMas)
    ELSE
      CALL do_velverlet(xyz, vel, acc, frc, AtMas, dt)
    END IF

! Thermostat - Scale velocity with designated T
    CALL do_scaleVel(vel, Temp, AtMas)

! Perodic boundary - Recursive
    CALL do_FixXYZ(xyz(1:nPart,1:3), vol)

! Calculate force
    CALL do_calcFrc(frc, EAM, xyz, AtMas, AtNums, atNumsVec, nAtms, EAMdata)

! Calculate total kinetic energy
    CALL do_calcKin(Kin, AtMas, vel)

! sample initial values
    E = sum(EAM(:)) + Kin
    IF (tstp==0) THEN
      E0 = E0
    END IF
    
    IF ( tstp==samptp ) THEN
      samp(sampK)%xyz(:,:) = xyz(:,:)    
      samp(sampK)%vel(:,:) = vel(:,:)
      samp(sampK)%acc(:,:) = acc(:,:)
      samp(sampK)%frc(:,:) = frc(:,:)
      samp(sampK)%EAM(:) = EAM(:)
      samp(sampK)%kin = kin
      samp(sampK)%E = E
      samp(sampK)%error = (E-E0)/E0
      samp(sampK)%time = tstp * dt
      IF (samptp+sampOfs <= nTstps) THEN
        samptp = samptp + sampOfs  
      ELSE
        samptp = nTstps  
      END IF    
        sampK = sampK + 1
        WRITE(*,*), 'timestep ', tstp, ' from ', nTstps, 'dt (sec)= ', dt, 't (sec)= ', nTstps * dt
    END IF
  END DO


open(unit=27, file="data.md", action="write", status="replace")
open(unit=26, file="readme.md", action="write", status="replace")
open(unit=25, file="energies.md", action="write", status="replace")
open(unit=24, file="pot.md", action="write", status="replace")
open(unit=23, file="frc.md", action="write", status="replace")
open(unit=22, file="acc.md", action="write", status="replace")
open(unit=21, file="vel.md", action="write", status="replace")
open(unit=20, file="xyz.md", action="write", status="replace")

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
l = 1
do jx = 1, npart
  do ix = 0, sampk-1
    write(20,*), k, ix, samp(ix)%time,  samp(ix)%xyz(jx,1:3)
    write(21,*), k, ix, samp(ix)%time,  samp(ix)%vel(jx,1:3)
    write(22,*), k, ix, samp(ix)%time,  samp(ix)%acc(jx,1:3)
    write(23,*), k, ix, samp(ix)%time,  samp(ix)%frc(jx,1:3)
    write(24,*), k, ix, samp(ix)%time,  samp(ix)%eam(jx)
    k= k+1
  end do
  write(25,*), l, ix, samp(ix)%time,  samp(ix)%E, &
               samp(ix)%Kin,                                  &
               samp(ix)%E - samp(ix)%Kin,    &
               samp(ix)%error                   
  l = l + 1 
end do

write(27,*), nPart, nSpcis, nAtms, sampk

CLOSE(UNIT=20)
CLOSE(UNIT=21)
CLOSE(UNIT=22)
CLOSE(UNIT=23)
CLOSE(UNIT=24)
CLOSE(UNIT=25)
CLOSE(UNIT=26)
CLOSE(UNIT=27)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END PROGRAM md
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

