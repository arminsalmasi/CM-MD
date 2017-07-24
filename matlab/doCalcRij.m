%% calculate distance between poit i and point j
function [d x y z] = doCalcR(a, i ,j)
 d = sqrt( (a.aPos(i,1)-a.aPos(j,1))^2 + ...
     (a.aPos(i,2)-a.aPos(j,2))^2 + ...
     (a.aPos(i,3)-a.aPos(j,3))^2 );
 x = a.aPos(i,1)-a.aPos(j,1);
 y = a.aPos(i,2)-a.aPos(j,2);
 z = a.aPos(i,3)-a.aPos(j,3);
end
