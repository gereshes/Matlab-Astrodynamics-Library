close all
clear all
clc
%% 
%This script is designed to provide an exmpale how to use different
%functions in the Matlab Astrodynamics Library CR3BP
%
%   Ari Rubinsztejn
%   www.gereshes.com
%   2019.03.06

%% Initialization

[mu,lStar,tStar,sec,prim] = stats2Body("Moon"); % Get some stats for the Earth-Moon System Note the use of a string

ic=[1-mu,.0455,0,-.5322,.2,0];% Define some initial conditions for our simulation (in the rotating reference frame)
ts=[0,20]; %Timespan of simulation in non-dimensional time
opts=odeset('RelTol',1e-12,'AbsTol',1e-12);
j0=jacobiValue3D(ic,mu);%Calcualtes the intial jacobi value


[t,y]=ode45(@(t,x) rotFrame3d(t,x,mu),ts,ic,opts); %Integrate the trajectory

diffJ = jacobiValDiff(y,mu);


%% Plotting the Trajectory and the Change of the Jacobi Value
h0=figure();
h0.Position=[220.200000000000,630,1083.80000000000,420];
subplot(1,2,1)
plot(y(:,1),y(:,2),'DisplayName','Spacecraft Trajectory')
plotOther(mu,j0,1) %Plot other
title('Trajectory in the Earth-Moon System')
xlabel('X-Axis (ND)')
ylabel('Y-Axis (ND)')
legend()
subplot(1,2,2)
plot(t,diffJ)
title('Change in Jacobi Value over time')
xlabel('Time (ND)')
ylabel('\Delta Jacobi Value')
grid on

%% Plotting the Change of Jacobi Value and Distance from Primary Body

h1=figure();
plot(t,abs(diffJ)./max(abs(diffJ)),'DisplayName','Absolute \Delta Jacobi Value','LineWidth',2)
hold on
plot(t,sqrt(((y(:,1)-mu).^2)+(y(:,2).^2))./max(sqrt(((y(:,1)-mu).^2)+(y(:,2).^2))),'--','DisplayName','Normalized Distance From Primary','LineWidth',2)
title('Normalized Absolute \Delta Jacobi Value Vs Normalized Distance From Primary')
xlabel('Time (ND')
legend('Location','southeast')
set(gca, 'YScale', 'log')

%% Plotting the Trajectory Alone
h=figure();
plot(y(:,1),y(:,2),'DisplayName','Spacecraft Trajectory')
plotOther(mu,j0,1) %Plot other
title('Trajectory in the Earth-Moon System')
xlabel('X-Axis (ND)')
ylabel('Y-Axis (ND)')
legend()