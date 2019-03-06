function [XNew] = rotFrame3dSTM(t,X,mu)
% CR3BP rotating frame  dynamics and State Transition matrix function to be used with Matlab ODE suite
% 
%Inputs:
%   t - double - current time
%   X - 42 by 1 vector - Current states [x,y,z,xDot,yDot,zDot,stmvector]'
%   mu- double - Non-Dimensional Mass ratio
%
%Outputs:
%   XNew - 42 by 1 vector - Derivative of curent states
%       [xDot,yDot,zDot,xDotDot,yDotDot,zDotDot,stmvectorDot]'
%
%   Ari Rubinsztejn
%   www.gereshes.com
%   2019.03.06

%% Unpack variables
x=X(1);
y=X(2);
z=X(3);
xDot=X(4);
yDot=X(5);
zDot=X(6);
phiMat=reshape(X(7:end),6,6);

%% Calcualte partial derivatives
d=sqrt(((x+mu)^2)+(y^2)+(z^2));
d3=d^3;
d5=d^5;
r= sqrt(((x-1+mu)^2)+(y^2)+(z^2));
r3=r^3;
r5=r^5;
sigXX= 1-((1-mu)/(d3)) - (mu/(r3))+...
    (3*(1-mu)*((x+mu)^2)/(d5))+...
    (3*mu*((x-1+mu)^2)/(r5));
sigYY=1-((1-mu)/(d3)) - (mu/(r3))+...
    (3*(1-mu)*(y^2)/(d5))+...
    (3*mu*(y^2)/(r5));
sigZZ=-((1-mu)/(d3)) - (mu/(r3))+...
    (3*(1-mu)*(z^2)/(d5))+...
    (3*mu*(z^2)/(r5));
sigXY=(3*(1-mu)*(x+mu)*y/(d5))+...
    (3*mu*(x-1+mu)*y/(r5));
sigXZ=(3*(1-mu)*(x+mu)*z/(d5))+...
    (3*mu*(x-1+mu)*z/(r5));
sigYZ=(3*(1-mu)*z*y/(d5))+...
    (3*mu*z*y/(r5));

%% Form STM
A = [0,0,0,1,0,0;...
     0,0,0,0,1,0;...
     0,0,0,0,0,1;...
    sigXX,sigXY,sigXZ,0,2,0;...
    sigXY,sigYY,sigYZ,-2,0,0;...
    sigXZ,sigYZ,sigZZ,0,0,0];

%% Calculate Updates
xDotDot = (2*yDot)+x - ((1-mu)*(x+mu)/d3)-(mu*(x-1+mu)/(r3));
yDotDot = (-2*xDot)+y - ((1-mu)*y/d3)-(mu*y/r3);
zDotDot = (-z*(1-mu)/d3) - (mu*z/r3);
XNew=[xDot,yDot,zDot,xDotDot,yDotDot,zDotDot]';
phiDotMat=A*phiMat;
phiDotVec=reshape(phiDotMat,36,1);

XNew=[XNew;phiDotVec];
end

