function [dy] = tbp_ode_perturbed (t, y, muP,kep)
%
% tbp_ode_perturbed.m - Differential equation for the perturbed two body problem 
%
% INPUTS:
%  y [6x1]              State of the system              [km,km/s]
%  muP                  Gravitational parameter          [km^3/s^2]
%  kep                  keplerian parametre
%
% OUTPUTS:
%  dy                   State of the system derivative   [km/s^2,km/s^3]
%      
%
% First Edition: 01/07/2022
%
% Authors: Rocco Larocca, Matteo Mascelloni, Davide Lanza, Afaq Shakeel

r = y(1:3);
v = y(4:6);

r_norm = norm(r);

Re = astroConstants(23); % mean radius
J2 = astroConstants(9);

kJ2 = 1.5*J2*muP*Re^2/r_norm^4;

a_J2_i = kJ2 * r(1)/r_norm*(5*r(3)^2/r_norm^2-1);
a_J2_j = kJ2 * r(2)/r_norm*(5*r(3)^2/r_norm^2-1);
a_J2_k = kJ2 * r(3)/r_norm*(5*r(3)^2/r_norm^2-3);

muMoon = astroConstants(20); 
kep_sc=kep;
[rr_sc,~] = kep2car(kep_sc,muP);
              
% Set of Cartesian position of MOON (XYZ) from EPHEMERIDES OF MOON:
t_MJD2000=t/86400; %[day] time from second to day of MJD2000 time
[rm, ~] = ephMoon(t_MJD2000); %ATTENTION !!! DEFINE 't' as MJD2000 time.
nrm=norm(rm); %[km] modulus of Moon's position wrt the Earth
%where: rm,vm : are respectively position and velocity vector of Moon wrt the Earth.

% Set of Cartesian position of MOON wrt the S/C (XYZ):
rms=rm'-rr_sc; %[km] Position vector of Moon wrt the s/c
nrms=norm(rms); %[km] modulus of Moon's position wrt the s/c

%Perturbation acceleration due to Moon gravity:
a_moon_XYZ = muMoon * ( (rms/(nrms^3) ) - (rm'/(nrm^3)) );

%Keplerian elements OF SPACECRAFT:
    a_sc = kep_sc(1);
    e_sc = kep_sc(2);
    i_sc = kep_sc(3);
    OM_sc = kep_sc(4);
    om_sc = kep_sc(5);
    th_sc = kep_sc(6);

%Rotation matrix (direction cosines) from XYZ (cartesian) frame to RSW frame, as function of Keplerian elements of SPACECRAFT:
u_sc=om_sc+th_sc; %[rad] argoument of latitude OF SPACECRAFT;
ROT_xyz_rsw=[-sin(OM_sc)*cos(i_sc)*sin(u_sc)+cos(OM_sc)*cos(u_sc), cos(OM_sc)*cos(i_sc)*sin(u_sc)+sin(OM_sc)*cos(u_sc), sin(i_sc)*sin(u_sc);
             -sin(OM_sc)*cos(i_sc)*cos(u_sc)-cos(OM_sc)*sin(u_sc), cos(OM_sc)*cos(i_sc)*cos(u_sc)-sin(OM_sc)*sin(u_sc), sin(i_sc)*cos(u_sc);
                          sin(OM_sc)*sin(i_sc)                   ,              -cos(OM_sc)*sin(i_sc)                 ,      cos(i_sc)    ];
%Computation of aj2_XYZ in RSW frame:
aMoon_RSW=ROT_xyz_rsw*a_moon_XYZ;
                
% Rotation matrix in Curtis pag. 499

aMoon_R= aMoon_RSW(1);
aMoon_S= aMoon_RSW(2);
aMoon_W= aMoon_RSW(3);


dy = [v(1); v(2); v(3)                   
     -muP/r_norm^3*r(1)+a_J2_i + aMoon_R 
     -muP/r_norm^3*r(2)+a_J2_j + aMoon_S
     -muP/r_norm^3*r(3)+a_J2_k + aMoon_W ];
    
end
