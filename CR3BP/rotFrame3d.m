function [XNew] = rotFrame3d(t,X,mu)
%CR3BP rotating frame dynamics function, to be used with Matlab ODE suite
% 
%Inputs:
%   t - double - current time
%   X - 6 by 1 vector - Current states [x,y,z,xDot,yDot,zDot]'
%   mu- double - Non-Dimensional Mass ratio
%
%Outputs:
%   XNew - 6 by 1 vector - Derivative of curent states
%       [xDot,yDot,zDot,xDotDot,yDotDot,zDotDot]'
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

d=sqrt(((x+mu)^2)+(y^2)+(z^2));
d3=d^3;
r= sqrt(((x-1+mu)^2)+(y^2)+(z^2));
r3=r^3;

%% Calculate Updates
xDotDot = (2*yDot)+x - ((1-mu)*(x+mu)/d3)-(mu*(x-1+mu)/(r3));
yDotDot = (-2*xDot)+y - ((1-mu)*y/d3)-(mu*y/r3);
zDotDot = (-z*(1-mu)/d3) - (mu*z/r3);
XNew=[xDot,yDot,zDot,xDotDot,yDotDot,zDotDot]';
end

