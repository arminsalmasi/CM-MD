%% functin dorandTypsAtoms - creats a vector which contains randomly
  % distributed N atoms based on types of atoms (atomic number):read from X
  % vector
% created by Armin Salamsi
% start date 2017-06-12

%% Input:
  % timeSpaceVec = first cell of array of datastructure
  % s [global] = rand stream crerated based on system clock
  % aTyp [global] = types of atoms in the systm
  % X [global] = concentration matrix
  % periodicData[global] data from priodic table 
    % .textData = AtomicNumber ElementName
    % .Data = AtomicR IonicR CovalentR VanderWaalsR CrystalR MolarMass
%% Output:
  % timeSpaceVec = first cell of array of datastructure
%%
function [timeSpaceVec] = doRandTypsAtoms(timeSpaceVec)
  global s periodicData
  
  % elN = ordered distribution of atoms based on concentration
    % vector
  % oAV = ordered (by type) vector = distribution of N atomes 
  % elN = Initialize ordered distribution of N atoms based on Types
    elN = round (timeSpaceVec.aXs  * timeSpaceVec.N);      
  % create Ordered Atom Vector
    oAV=[]; 
    for i = 1 : size(timeSpaceVec.aNums,2)
      oAV = [oAV ones(1,elN(1,i))*timeSpaceVec.aNums(1,i) ];
    end
   % Shuffle oAV to Random Atom Vector based on Types
    rAV =  oAV(randperm(s,size(oAV,2))); 
   % Report function value
    timeSpaceVec.aTyps = rAV';
    for  i = 1 : timeSpaceVec.N
      % set atomic mass of atoms    
      timeSpaceVec.aMass (i ,1)  = periodicData.data(rAV(1,i),6);
      % set atomic radi of atoms    
      timeSpaceVec.aRadi (i ,1)  = periodicData.data(rAV(1,i),2);
    end
    
end
