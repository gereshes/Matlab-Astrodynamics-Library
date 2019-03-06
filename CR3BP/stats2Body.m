function [mu,lStar,tStar,sec,prim] = stats2Body(secondaryNameOrIndex)
%This function takes in the name of a celestial object and returns a cell
%array of pertinent information. 
%
%Note: See https://gereshes.com/2018/11/12/dynamics-of-the-3-body-problem/
%   for eplanations/calculations of mu,lStar,tStar
%
%Input:
%   bodyNumberOrName - String or Int - Body name with first letter capitalized or index number
%
%Output:
%   mu- double - Nondimensional mass ratio
%   lStar - double - Charachtersitich length (Km)
%   tStar - double - Charachteristich time (s)
%   sec - 1x9 Cell array - It contains the following fields: Body
%   name, Body it's orbiting around, Axial Rotational rate (Rev/Day),
%   Equatorial Radius (Km),Gravitational Parameter (Km**3 s**-2),
%   Semi-major axis (Km), Orbital Period (days), Eccentricity, Inclination
%   of orbit to elliptic (deg), Mass of the planet (Kg)
%   prim - 1x9 Cell array - It contains the following fields: Body
%   name, Body it's orbiting around, Axial Rotational rate (Rev/Day),
%   Equatorial Radius (Km),Gravitational Parameter (Km**3 s**-2),
%   Semi-major axis (Km), Orbital Period (days), Eccentricity, Inclination
%   of orbit to elliptic (deg), Mass of the planet (Kg)
%                              
%   Ari Rubinsztejn
%   www.gereshes.com
%   2019.03.06

G=6.67408E-20; %Universal gravatational constant (Km**3 Kg**-1 s**-2)
sec=celestialLookUp(secondaryNameOrIndex);
prim=celestialLookUp(string(sec{2}));
mSec  = sec{5}/G;
mPrim = prim{5}/G;
prim{10}=mPrim;
sec{10}=mSec;
mStar=mSec+mPrim;
mu=mSec/mStar;
lStar=sec{6};
tStar=sqrt((lStar.^3)./(G.*mStar));
end

