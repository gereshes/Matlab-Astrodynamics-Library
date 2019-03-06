function [diffJ] = jacobiValDiff(X,mu)
%This function calculates the Jacobi value difference in the Circular
%Restricted 3-Body Problem (CR3BP)
%
%Inputs
%   X  - n x 6 Vector - Your current states [x,y,z,xDot,yDot,zDot]
%   mu - double     - Your systems nondimenisonal mass ratio
%
%Outputs
%   diffJ - n x 1 Vector - The differnece of the jacobi constant over the
%                       trajectory
%
%Notes:For background on the CR3BP I suggest going here 
%      https://gereshes.com/2018/11/12/dynamics-of-the-3-body-problem/
%           or
%      https://gereshes.com/category/math/astrodynamics/cr3bp/
%
% Ari Rubinsztejn
% www.gereshes.com
% 2019.03.06


[entries,~]=size(X);
j0=jacobiValue3D(X(1,:),mu);
difJ=NaN(entries,1);
for c=1:entries
    diffJ(c)=jacobiValue3D(X(c,:),mu)-j0;
end

end

