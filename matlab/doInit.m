%% functin doInit initialize data structure with random numbers
% created by Armin Salamsi
% start date 2017-06-09

%create data structure - data structure is not initialaized
%Input:
    % N  = number of atoms in cell [#]
    % T  = temperature [K]
    % t  = simualtion time [sec]
    % dt = lenght of timesteps [real: seconds]
    % aTyp = array of types of atoms / orderd by atomic number [# ...]
    % X = concentration of each atom type / ordered by atomic number 
        %[#real ...]
    % Vol = equilibrium volume of N atoms at T from TC in [m^3]
    % ND = number of dimensions 
    % NA = Avogadro number [mol^-1].
    % KB = Boltzmann constant [m^2.kg.s^-2.K^-1].
%Output:
    % timeSpaceVec : first cell of array of struct:
    % timeSpaceVec.tstep  = timestep #         [tstp         : integer].
    % timeSpaceVec.aIdx   = Index of atoms     [idx          : integer].
    % timeSpaceVec.aTyps  = Atomic number      [atN ...      : integer].
    % timeSpaceVec.aPos   = Pos of Atoms       [x y z ...    : real].
    % timeSpaceVec.aVel   = Velocities.        [vx vy vz ... : real].
    % timeSpaceVec.aAcc   = Accelerations.     [ax ay az ... : real].
    % timeSpaceVec.aforce = Interatomic Force. [fx fy fz ... : real].
    % timeSpaceVec.aMass  = Atomic Masses.     [atM ...      : real]. 
    % timeSpaceVec.aRadi  =Atomic radiuses.    [atR ...      : real]. 
    % timeSpaceVec.boxV   = Volume of box.     [boxVol       : real]
    % timeSpaceVec.boxL   = Length of box.     [boxLenght    : real]
    % timeSpaceVec.NP     = # of atoms in box. [Natm         : integer]

function [a] = doInit(a, V, N, ND, tstep, T);
    global NA
    %% set timestep number to 1
      a(1).tstep = tstep;
      
    %% set volume of the space domain
      a(1).boxV = V* N/ NA;
    
    %% set lenght of the space domain
      a(1).boxL = (V* N/ NA)^(1/ ND);
      
    % randomize distribution of N atoms by types of atoms [atomic number]
      [a(1)] = doRandTypsAtoms(a);
%plot(timeSpaceVec.aTyps,'*')

      %% randomize positions of atoms in the box
      [a(1)] = doRandPosAtoms(a);
%scatter3(timeSpaceVec.aPos(:,1),timeSpaceVec.aPos(:,2)...
%,timeSpaceVec.aPos(:,3))

    
    %% randomize velocities : prescribe a random velocity vector based on
      % the Temperature and  standard deviation of velocities 
      % Vrand = rand(1) * sum(sqrt( T * KB  / Mass(i))
      [a] = doRandVel(a);
      %plotVel(timeSpaceVec);  %plot velocity distribution and sum 
      doTCalc(a) % calculate temperature from velocities
      
    
    %% fix center of mass :Calcualte velocity of the center of mass 
      % VCM = sum(mi*Vi)/sum(mi) then Fix any center-of-mass drift 
      % V= V - VCoM

      [a] = doFixCoM(a);
      %plotVel(timeSpaceVec);   %plot velocity distribution and sum  
      doTCalc(a) % calculate temperature from velocities     
    
    %% scale velocities by the prescribed temperature factor = T/T0
      [a] = doScaleVel(a);
      %plotVel(timeSpaceVec); %plot velocity distribution and sum   
      doTCalc(a) % calculate temperature from velocities
    
    %% calculate initial forces and accelerations!
    
      
      
end


