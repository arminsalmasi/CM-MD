%% functin plotVel - plot velocities 

% created by Armin Salamsi
% start date 2017-06-12

% T = sum( Mass(i) *V(i)) /(KB *Nf) ) 

%% Input
  % a = timeSpaceVec = first element of datastructure array

%% Output
  % calculated temperature after correction itterations

%%
function []=plotVel(timeSpaceVec)
  %% plot sorted velocities
  figure
    hold on
      plot(sort(timeSpaceVec.aVel(:,1)));
      plot(sort(timeSpaceVec.aVel(:,2)));
      plot(sort(timeSpaceVec.aVel(:,3)));
      for  i = 1 : timeSpaceVec.N
        b(i) = sqrt( dot(timeSpaceVec.aVel(i,:), timeSpaceVec.aVel(i,:)));
      end
      plot(sort (b))
    hold off
  %% plot unsorted velocities
  figure
    hold on
      plot(timeSpaceVec.aVel(:,1));
      plot(timeSpaceVec.aVel(:,2));
      plot(timeSpaceVec.aVel(:,3));
      plot(b)
    hold off
end