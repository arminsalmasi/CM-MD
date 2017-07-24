%% function doLoop : loops over time steps. 
% created by Armin Salamsi
% start date 2017-06-09

% 1- doInit: initialaize first timestep.
  %a-random distribution of different types of atom - atomic weight and
    %position
  %b-assign random distribution of velocities atomes bsed on Gausian
    %which is generated based on the temperature and kinetic energy
% 2- doMD loops over all time steps and return the results in timeSpaceVec
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

function [tSV]=doLoop(DataStruct ,V ,N ,ND ,t ,dt ,T) 
  %% initialize data structure with random numbers
    % set volume
    % set box lenght
    % randomize atoms
    % randomize positions
    % randomize velocities
    % scale velocities by energy
  %% Input: 
    % V =volume of the cell [m^3]
    % N = number of atoms in the cell
    % ND = number of dimensions
    % tstp = timestep number
    % T = teperature [K]
  %% output:
    % timeSpaceVec : array of struct containign the MD results

      
      
  %%
    for stp = 1 : floor(t/dt)-1
      if stp == 1
          [ tSV(1) ] = doInit(DataStruct, V, N, ND, 0, T); %tSC = timespaceVector
          tSV(1).tstep = stp;
          %tSV(stp).aforce = doPotCalcLJ(tSV(stp), 1.645e-21 , 3.405e-10);
          %???
      else
        tSV(stp) = tSV(stp-1);
        
        tSV(stp).tstep = stp;
        
        [tSV(stp).aPos, tSV(stp).aVel, tSV(stp).aAcc] =...
        doVelVerlet( tSV(stp-1).aPos, tSV(stp-1).aVel, tSV(stp-1).aAcc,...
        tSV(stp-1).aforce, tSV(stp-1).aMass, dt );
        
        %tSV(stp).aforce = doPotCalcLJ(tSV(stp), 1.645e-21 , 3.405e-10);
        tSV(stp).aforce = LJ_Force(tSV(stp).aPos, tSV(stp).boxL,...
            1.645e-21 , 3.405e-10 );
        
      end
      %figure(stp)
      %scatter3(tSV(stp).aPos(:,1),tSV(stp).aPos(:,2),tSV(stp).aPos(:,3))
  end
  
%end

%%%%timeSpaceVec(tstep)  = doMD(timeSpaceVec(tstep-1), 0, dt);
%% doMD loop over time steps
      %(except for 1 and calculate the interactions)   
      %doMD(timeSpaceVec, tstep , dt);
      
      
      
%       timeSpaceVec(tstep).tstep = tstep;
%       if tstep == 1 
%         timeSpaceVec(1)  = doMD(initDataStruct, 1, dt);
%       else
%       % update timestep index from n to n+1
%         
%       % calculate potentials - update force vector
%         timeSpaceVec(tstep) = doPotCalcLJ(timeSpaceVec(tstep), 1.645e-21 , 3.405e-10 );
%         N = timeSpaceVec(tstep).N;
%         ND = timeSpaceVec(tstep).ND;
%         for i = 1 : timeSpaceVec(tstep).N
%             timeSpaceVec(tstep).aAcc(i,:) = timeSpaceVec(tstep).aforce(i,:) / timeSpaceVec(tstep).aMass(i);
%         end
%         timeSpaceVec(tstep).aVel = timeSpaceVec(tstep-1).aVel + ( 0.5 * (timeSpaceVec(tstep-1).aAcc + timeSpaceVec(tstep).aAcc) ) * dt ;

      