% Calculate leonard Jons Potentials
function  [Fcomp] = doPotCalcLJ(a, eps, sig)
  
for i = 1 : a.N-1
  for j= i+1 : a.N
    %% update forces
    if i ~= j
      [rij x y z]= doCalcRij(a, i ,j);
      F = ( 24 * eps / sig ) * (2 * (sig/rij)^13 - (sig/rij)^7 ); 
      Fcomp(i,:) =  F * (a.aPos(i,:) - a.aPos(j,:));
    end
  end
end
  
end