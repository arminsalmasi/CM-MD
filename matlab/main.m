% Created by : Armin Salmasi.
% Date Start: 2017-06-09.
% For : Project in course: computational methods in materials science.
% Pupropse : to calculate MD simple case for a ternary alloy.
% Test case : Co Fe Ni Cr [0.25 0.25 0.25 0.25]. 

%% Main body MD Project. 
  
  clear variables
  clc
  close all

%% set global parameter
    global NA KB X s periodicData
      % Avogadro number [mol^-1].
      NA = 6.02214179e23; 
      % Boltzmann constant [m^2.kg.s^-2.K^-1]. 
      KB = 1.38064852e-23; 
      % s = Randstream by system clock
      s = RandStream.setGlobalStream(RandStream('mt19937ar',...
        'seed',sum(100*clock)));
      % specification of elements from periodic table:
        % .textData = AtomicNumber ElementName
        % .Data = AtomicR IonicR CovalentR VanderWaalsR CrystalR MolarMass
      periodicData = importdata('atmdata.txt');
      
%% read initial data - IO 
%Input
    
%Output:
  % N  = number of atoms in cell [#].
  % T  = temperature [K].
  % t  = simualtion time [sec].
  % dt = lenght of timesteps [#].
  % aTyp = array of types of atoms / orderd by atomic number [# ...].
  % X = concentration of each atom type / ordered by atomic number. 
   %[#real ...].
  % Vol = equilibrium volume of N atoms at T from TC in [m^3].
  % ND = number of dimensions.
    
  [ N aXs T t dt aNums V ND] = getInitData;

%% create data structure - DO    
%Input:
  % N  = number of atoms in cell [#].
  % T  = temperature [K].
  % t  = simualtion time [sec].
  % dt = lenght of timesteps [real: seconds].
  % aTyp = array of types of atoms / orderd by atomic number [# ...].
  % X = concentration of each atom type / ordered by atomic number. 
   %[#real ...].
  % Vol = equilibrium volume of N atoms at T from TC in [m^3].
  % ND = number of dimensions.    
%Output:
  % timeSpaceVec : first cell of array of struct:
  % timeSpaceVec.tstep  = timestep #.        [tstp         : integer].
  % timeSpaceVec.aIdx   = Index of atoms.    [idx          : integer].
  % timeSpaceVec.aTyps  = Atomic number.     [atN ...      : integer].
  % timeSpaceVec.aPos   = Pos of Atoms.      [x y z ...    : real].
  % timeSpaceVec.aVel   = Velocities.        [vx vy vz ... : real].
  % timeSpaceVec.aAcc   = Accelerations.     [ax ay az ... : real].
  % timeSpaceVec.aforce = Interatomic Force. [fx fy fz ... : real].
  % timeSpaceVec.aMass  = Atomic Masses.     [atM ...      : real]. 
  % timeSpaceVec.aRadi  =Atomic radiuses.    [atR ...      : real]. 
  % timeSpaceVec.boxV   = Volume of box.     [boxVol       : real]
  % timeSpaceVec.boxL   = Length of box.     [boxLenght    : real].
  % timeSpaceVec.NP     = # of atoms in box. [Natm         : integer].
  
  [ DataStruct ] = doMakeDataStruct( N, ND, aNums, aXs, T);
    
%% Loop over time steps - Do 
  % 1- doInit: initialaize first timestep.
   %a-random distribution of different types of atom - atomic weight and
    %position
   %b-assign random distribution of velocities atomes bsed on Gausian
    %distribution which is generated based on the temperature and 
    %kinetic energy 
    % fixing the center of mass
    % rescaling the velocities to fit with the prescribed temperature
  % 2- doMD loops over all time steps and return the results : timeSpaceVec
   %doMD contains: a-upadate of forces b-velocities c-coordinates
    %d-thermostat
%Input:
  % timeSpaceVec : first cell of array of struct
  % V = equilibrium volume of N atoms at T from TC in [m^3]
  % ND = number of dimensions
  % N  = number of atoms in cell [#]
  % T  = temperature [K]
  % t  = simualtion time [sec]
%Output:
   % timeSpaceVec : array of struct containign the MD results

  [timeSpaceVec] = doLoop( DataStruct, V, N, ND, t, dt, T);
    
%% export data
