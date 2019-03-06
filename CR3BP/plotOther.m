function plotOther(mu,jacobiValue,filled)
%Function for plotting other gerneal system points of interest.
%
%This function takes in system and trajectoy parameters, computes other
%objects of interest: 
%   1) Positions of the bodies
%   2) Zero Velocity Region (ZVC) in the X-Y plane
%   3) Locations of the Lagrange points
%and then plots them on the last graph
%   Inputs:
%       mu - double - Nondimenstional mass ratio
%       jacobiValue - double - Jacobi clue of desired ZVC (put 0 if no ZVC
%           desired)
%       filled - int - 0 if you want unfilled ZVC's, 1 if you want it
%           filled, 2 if you want a solid fill
%   Outputs:
%       None
%
%   Ari Rubinsztejn
%   www.gereshes.com
%   ari@gereshes.com
%   2019.03.06

gcf();
hold on


if (jacobiValue==0) %Check if you want a ZVC
   
else
    samples=1000;
    desJac= jacobiValue;
    x=linspace(-1.5,1.5,samples);
    y=linspace(0,1.5,samples);
    holder=NaN(samples,samples);
    holderPos=zeros(samples,samples);
    eps=.02;
    eps=.1;

    for c=1:samples
        for d=1:samples
            rVec=[x(c),y(d),0];
            tm=jacobiConstZVC(mu,rVec);
            if((desJac)>=tm && tm>(desJac*(1-eps)))
                
                holder(d,c)=tm;
                holderPos(d,c)=1;
                
            else
                
            end
        end
    end
    if(filled==1)
        contourf(x,y,holder,'DisplayName','ZVC','edgecolor','none')
        contourf(x,-y,holder,'HandleVisibility','off','edgecolor','none')
    elseif(filled==2)
        contourf(x,-y,holder*10,'HandleVisibility','off','edgecolor','none')
        contourf(x,y,holder,'DisplayName','ZVC','edgecolor','none')
        contourf(x,-y,holder,'HandleVisibility','off','edgecolor','none')
    else
        contour(x,y,holder,'DisplayName','ZVC')
        contour(x,-y,holder,'HandleVisibility','off')
    end
    
    
end


[xPos,yPos,rSec,rPrim]=lagrangePoints(mu,50);
scatter(xPos,yPos,'+','DisplayName','Lagrange Points')
hold on
% rSec
scatter(1-mu,0,rSec*100,'filled','DisplayName','Secondary','MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 1 0])
scatter(-mu,0,rPrim*100,'filled','DisplayName','Primary','MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[1 0 0])
grid on
grid minor
axis([-1.5,1.5,-1.5,1.5])
daspect([1,1,1])
end


function [jacobiConstVal] = jacobiConstZVC(mu,rVec)
%Calculates the Jacobi constant at a point
x=rVec(1);
y=rVec(2);
z=rVec(3);
d=sqrt(((x+mu)^2)+(y^2)+(z^2));
r=sqrt(((x-1+mu)^2)+(y^2)+(z^2));
jacobiConstVal=(x^2)+(y^2)+((2-(2*mu))/d)+(2*mu/r);
end

