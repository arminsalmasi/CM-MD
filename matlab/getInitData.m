% get IO initial data from Input
% created by Armin Salamsi
% start date 2017-06-09
% read initial data - IO 
    % N = number of atoms in cell
    % T = temperature
    % t = simualtion time
    % dt = lenght of timesteps
    % idxA = array of types of atoms / orderd by atomic number
    % X = concentration of each atom type / ordered by atomic number
    % eqVol = equilibrium volume of N atoms at T from TC
    % NA (global) = Avogadro number [mol^-1].
function [ N X T t dt aTyp V ND] = getInitData
    global NA 
    
    % N = number of atoms in the cell
    N = 100; 
    % X = concentration in mole fraction
    X = [1]; %Argon  %[0.25 0.25 0.25 0.25];  
    % T = temperature in Kelvin
    T = 298.15;%298.15; 
    % t = Simulation time in seconds
    t = 1e-12; 
    %dt = leght of timestep timesteps
    dt = 1e-15; 
    % aTyp = Type of atoms orderd with atomic number Cr Fe Co Ni 
    aTyp = [18]; %Argon %[24 26 27 28]; 
    % V = equilibrium Volum of N atom - read from TC at:
    % T = 298.15, N=1, X(FE)=0.25, X(CO)=0.25, X(NI)=0.25, P=1E5
    V = 6.9601534e-6 * N / NA; 
    % Number of dimenstions
    ND = 3;
  
%x = input(prompt)
%str = input(prompt,'s')
end