clc, close all, clear variables

load('xyz.md')
load('dat.md')
time(:,1) = xyz(:,4);
pos(:,1:3) = xyz(:,5:7);

nSamps = dat(1,end);
np = dat(1,1);



% path of all particles 
% or selected particle by spcifying i-1
figure
for i = 0 : np-1
   Ps = nSamps * i + 1;
   Pe   = nSamps * (i + 1);
   pos_atm(i+1,:,1:3) = pos(Ps:Pe,1:3); %pos of atom i in all tstps samps
   h = plot3( pos(Ps:Pe,1), pos(Ps:Pe,2), pos(Ps:Pe,3));
   hold on
end
title('path of all particles')
grid on
% pos of all atoms in a given sampel numbers
% [samp_number = 0 == (random positions)]
% max = nSamps
for samp_number =  [1, nSamps]
    samp_number
    figure
    for i = 0 : np-1
        Ps = nSamps * i + samp_number;
        if (i+1)<=dat(1,3)
            h = plot3(pos(Ps,1), pos(Ps,2),pos(Ps,3), 'O', 'MarkerSize',10);
            hold on
                title('pos of all atoms in given timestep sampeles' )
            grid on
        else
            h = plot3(pos(Ps,1), pos(Ps,2),pos(Ps,3),'*', 'MarkerSize',10 ) ;
                title('pos of all atoms in given timestep samples' )
            grid on
        end
    end
end
