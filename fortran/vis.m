clc, close all, clear variables

load('xyz.md')
load('data.md')
time(:,1) = xyz(:,3);
pos(:,1:3) = xyz(:,4:6);

num_samples = data(1,end);
np = data(1,1);
figure
for i = 1 : np
    if i<data(1,3)
        h = plot3(pos((i-1)*num_samples+1,1), pos((i-1)*num_samples+1,1:2),pos((i-1)*num_samples+1,3),'O' ) 
        hold on
    else
        h = plot3(pos((i-1)*num_samples+1,1), pos((i-1)*num_samples+1,1:2),pos((i-1)*num_samples+1,3),'*' )
    end
end

figure
for i = 1 : np
    pos_atm_1(i,:,1:3) = pos((i-1)*num_samples+1:i*num_samples,1:3); 
    h = plot3(pos_atm_1(i,:,1), pos_atm_1(i,:,2),pos_atm_1(i,:,3) ) 
    hold on
end


figure
for i = 1 : np
    if i<data(1,3)
        h = plot3(pos(i*num_samples,1), pos(i*num_samples,1:2),pos(i*num_samples,3),'O' ) 
        hold on
    else
        h = plot3(pos(i*num_samples,1), pos(i*num_samples,1:2),pos(i*num_samples,3),'*' ) 
    end
end


%pos_atm_t(:,1:3) = pos(1:offset:size(pos,1) ,1:3); %% position of all atoms at timestep 1(aka0) time 1
    %pos_atm_t(:,1:3) = pos(5000:offset:size(pos,1) ,1:3); %% position of all atoms at timestep 1(aka0) time 2
    %figure
    %h = plot3(pos_atm_t1(1:15,1), pos_atm_t1(1:15,2),pos_atm_t1(1:15,3) ,'*') 
    %hold on
    %h = plot3(pos_atm_t1(16:20,1), pos_atm_t1(16:20,2),pos_atm_t1(16:20,3), 'o' ) 

    %figure
    %h = plot3(pos_atm_t2(1:15,1), pos_atm_t2(1:15,2),pos_atm_t2(1:15,3) ,'*') 
    %hold on
    %h = plot3(pos_atm_t2(16:20,1), pos_atm_t2(16:20,2),pos_atm_t2(16:20,3), 'o' ) 
    
    a =1