function [jacobiConst] = jacobiValue3D(X,mu)
%This function calculates the jacobi value of a state in the Circular
%Restricted 3-Body Problem (CR3BP)
%
%Inputs
%   X  - 1x6 Vector - Your current states [x,y,z,xDot,yDot,zDot]
%   mu - double     - Your systems nondimenisonal mass ratio
%
%Outputs
%   jacobiConst - double - The jacobi constant
%
%Notes:For background on the CR3BP I suggest going here 
%      https://gereshes.com/2018/11/12/dynamics-of-the-3-body-problem/
%           or
%      https://gereshes.com/category/math/astrodynamics/cr3bp/
%
% Ari Rubinsztejn
% www.gereshes.com
% 2018.11.22

x=X(1);
y=X(2);
z=X(3);
xDot=X(4);
yDot=X(5);
zDot=X(6);
r=sqrt(((x-1+mu).^2)+(y.^2)+(z.^2));
d=sqrt(((x+mu).^2)+(y.^2)+(z.^2));
jacobiConst = (x.^2)+(y.^2) +(2*(1-mu)./d)+(2*mu./r)-((xDot.^2)+(yDot.^2)+(zDot.^2));
end

