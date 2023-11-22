function aMoon_RSW = kepMoon_pert_RSW (t,kep,muMoon,muP)
% 
% kepMoon_pert_RSW.m - computes the pertubation accelaration, due to Moon's gravity, 
% in the RSW frame. 
%
% PROTOTYPE:
%  aMoon_RSW = kepMoon_pert_RSW (t,kep,muMoon,muP)
%
% INPUT:
%  t                      Time                                     [s]
%  kep [1x6]              Kepler elements                       
%  muMoon [1]             Moon's gravitational parameter           [km^3/s^2]
%  muP [1]                Earth's gravitational parameter          [km^3/s^2] 
%
% OUTPUT:
%  aMoon_RSW              Perturbed Moon acceleration in RSW frame [km/s^2]
%
% CONTRIBUTORS:
%  Lanza Davide
%  Larocca Rocco
%  Mascelloni Matteo 
%  Shakeel Afaq
% 
% VERSION:
%  07-01-2022


[r_sc,~] = kep2car(kep,muP);
              
t_days=t/86400; 
[r_moon, ~] = ephMoon(t_days); 
r_moon_norm=norm(r_moon); 


r_moon_sc=r_moon'-r_sc; 
r_moon_norm_sc=norm(r_moon_sc); %[km] modulus of Moon's position wrt the s/c

% Perturbation acceleration due to Moon gravity:
a_moon_XYZ = muMoon * ( (r_moon_sc/(r_moon_norm_sc^3) ) - (r_moon'/(r_moon_norm^3)) );

% Keplerian elements of s/c:
a_sc = kep(1); e_sc = kep(2); i_sc = kep(3);
OM_sc = kep(4); om_sc = kep(5); th_sc = kep(6);

% DCM from Cartesian to RSW  (pag 499 Curtis PDF)
w_sc=om_sc+th_sc; 
DCM_cartorsw=[-sin(OM_sc)*cos(i_sc)*sin(w_sc)+cos(OM_sc)*cos(w_sc), cos(OM_sc)*cos(i_sc)*sin(w_sc)+sin(OM_sc)*cos(w_sc), sin(i_sc)*sin(w_sc);
             -sin(OM_sc)*cos(i_sc)*cos(w_sc)-cos(OM_sc)*sin(w_sc), cos(OM_sc)*cos(i_sc)*cos(w_sc)-sin(OM_sc)*sin(w_sc), sin(i_sc)*cos(w_sc);
             sin(OM_sc)*sin(i_sc), -cos(OM_sc)*sin(i_sc),cos(i_sc)];


aMoon_RSW=DCM_cartorsw*a_moon_XYZ;      


end