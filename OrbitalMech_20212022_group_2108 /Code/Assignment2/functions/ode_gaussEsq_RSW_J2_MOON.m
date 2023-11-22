function d_kep = ode_gaussEsq_RSW_J2_MOON (t,kep,muP,aJ2_RSW,aMoon_RSW)
% 
% ode_gaussEsq_RSW_J2_MOON.m - computes planetary Gauss equations due to Moon's gravity and J2 perturbation
%                              in the RSW frame.
% 
% PROTOTYPE:
%  d_kep = ode_gaussEsq_RSW_J2_MOON (t,kep,muP,aJ2_RSW,aMoon_RSW)
% 
% INPUT:
% t [1]                      Time                                            [s]  
% kep [1x6]                  Kepler elements 
% muP [1]                    Earth's gravitational parameter                 [km^3/s^2]
% aJ2_RSW [1x3]              J2 Perturbed acceleration                       [km/s^2]
% aMoon_RSW [1x3]            Moon Perturbed acceleration                     [km/s^2]
% 
% OUTPUT:
% d_kep [6x1]                Derivatives of kepler elements
%
% AUTHORS:
%  Lanza Davide
%  Larocca Rocco
%  Mascelloni Matteo 
%  Shakeel Afaq
% 
% VERSION:
%  07-01-2022

% Keplerian parameters and orbit relations
a=kep(1);      e=kep(2);     i=kep(3);     
OM=kep(4);    om=kep(5);    th=kep(6);    

p=a*(1-e^2);                % Semi-lactus rectum
r=p/(1+e*cos(th));          % Position
b=a*sqrt(1-e^2);            % Semi-minor axis
n=sqrt(muP/a^3);            % Mean motion
h=n*a*b;                    % Angular momentum

% Components of AJ2_RSW:
aJ2_R= aJ2_RSW(1);
aJ2_S= aJ2_RSW(2);
aJ2_W= aJ2_RSW(3);

% Components of a_Moon_RSW:
aMoon_R= aMoon_RSW(1);
aMoon_S= aMoon_RSW(2);
aMoon_W= aMoon_RSW(3);

% Components of perturbed acceleration:
a_pert_R= aJ2_R + aMoon_R;
a_pert_S= aJ2_S + aMoon_S;
a_pert_W= aJ2_W + aMoon_W;

% Gauss equations in RSW frame:
da = (2*(a^2)/h) * (e*sin(th)*a_pert_R + (p/r)*a_pert_S);
de = (1/h) * ( p*sin(th)*a_pert_R + ( (p+r)*cos(th) + r*e )*a_pert_S );
di = ( r*cos(th+om)/h ) * a_pert_W;
dOM = ( r*sin(th+om)/(h*sin(i)) ) * a_pert_W;
dom = (1/(h*e)) * ( -p*cos(th)*a_pert_R + (p+r)*sin(th)*a_pert_S ) - ( (r*sin(th+om)*cos(i)) / (h*sin(i)) )*a_pert_W ;
dth =  (h/(r^2)) + (1/(e*h)) * ( p*cos(th)*a_pert_R - (p+r)*sin(th)*a_pert_S );

% Assemblying output
d_kep=[da;de;di;dOM;dom;dth];
