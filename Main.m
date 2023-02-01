clc; clear; close;
% Properties
m1 = 1; % kg
m2 = 0.1; % kg
l1 = .6; % m
l2 = .3; % m
g = 9.807; % m/s2
% Initial Conditions
theta1 = pi*0.95; % rad
theta2 = pi; % rad
pth1 = 0; % kgm2/s
pth2 = 0; % kgm2/s
% Setup
time_0 = 0; % s
time_step = 1e-4; % s
time_end = 5; % s
offset_x = (l1+l2)*1.05;
offset_y = (l1+l2)*1.05;
% Preallocate
timeVec = time_0:time_step:time_end;
u = nan(length(timeVec),4);
u_dot = nan(length(timeVec),4);
loc = nan(length(timeVec),4);
% Assign
u(1,:) = [theta1 theta2 pth1 pth2];
u_dot(1,:) = [0 0 0 0];
loc(1,:) = u(1,:);
w1d_f =@(t1,t2,w1,w2) (-g*(2*m1+m2)*sin(t1)-m2*g*sin(t1-2*t2)-2*sin(t1-t2)*m2*(w2^2*l2+w1^2*l1*cos(t1-t2)))...
    /(l1*(2*m1+m2-m2*cos(2*t1-2*t2)));
w2d_f =@(t1,t2,w1,w2)  (2*sin(t1-t2)*(w1^2*l1*(m1+m2)+g*(m1+m2)*cos(t1)+w2^2*l2*m2*cos(t1-t2))) ...
    /(l2*(2*m1+m2-m2*cos(2*t1-2*t2)));

% DO LOOP
for iTime = 2 : length(timeVec)
    t1 = u(iTime-1,1);
    t2 = u(iTime-1,2);
    w1 = u(iTime-1,3);
    w2 = u(iTime-1,4);
    
    u_dot(iTime,:) = [w1 w2 w1d_f(t1,t2,w1,w2) w2d_f(t1,t2,w1,w2)];
    
   	u(iTime,:) = u(iTime-1,:) + u_dot(iTime-1,:) * time_step;
    
    t1 = u(iTime,1);
    t2 = u(iTime,2);
    w1 = u(iTime,3);
    w2 = u(iTime,4);
    
    loc(iTime,1) = offset_x + l1*sin(t1);
    loc(iTime,2) = offset_y - l1*cos(t1);
    loc(iTime,3) = offset_x + l1*sin(t1) + l2*sin(t2);
    loc(iTime,4) = offset_y - l1*cos(t1) - l2*cos(t2);
end

LocationData = timeseries(loc,timeVec);
tx = 2:iTime;
hold on;
plot(offset_x,offset_y,'ok')
plot(loc(tx,1),loc(tx,2),loc(tx,3),loc(tx,4));
xlim([0 2*offset_x]);
ylim([0 2*offset_y]);
