function [xPos,yPos,rHill,rHillPrim] = lagrangePoints(mu,itts)
%This function calucaltes the positions of the lagrange points in the rotating frame
%for the Circular Restricted 3-Body Problem (CR3BP) using the hill radius
%approximation
%
%Inputs
%   mu  - double     - Your systems nondimenisonal mass ratio
%   itts- integer    - Number of itterations you 
%
%Outputs
%   xPos      - 5x1 Vec - x Positions of Lagrange points
%   yPos      - 5x1 Vec - y Positions of Lagrange points
%   rHill     - double  - Hill radius of the secondary body
%   rHillPrim - double  - Hill radius of the primary body
%
%Notes:For background on the CR3BP I suggest going here 
%
%      https://gereshes.com/category/math/astrodynamics/cr3bp/
%
% Ari Rubinsztejn
% www.gereshes.com
% 2019.03.06

rHill = (mu/3)^(1/3);
rHillPrim = ((1-mu)/3)^(1/3);
gam1=1;
gam2=1;

for c=1:itts
    gam1=gam1-((rHill/3)^c);
    gam2=gam2-((-1^c)*((rHill/3)^c));
end
gam1=rHill*gam1;
gam2=rHill*gam2;
xL1=(1-mu)-gam1;
xL2=(1-mu)+gam2;
xL3=-(1+(5*mu/12)); 
xL4=.5-mu;
xL5=xL4;

yL1=0;
yL4=sqrt(3)/2;

xPos=[xL1,xL2,xL3,xL4,xL5]';
yPos=[yL1,yL1,yL1,yL4,-yL4]';


end

