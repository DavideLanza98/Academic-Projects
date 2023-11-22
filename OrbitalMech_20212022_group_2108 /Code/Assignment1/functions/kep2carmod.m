function [r,v]=kep2carmod(k,mu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% kep2carmod:  function useful to convert from Keplerian elements to Cartesian coordinates
%
% INPUTS:
%
%  a         [1x1] Semi-majopr axis        [km]
%  e         [1x1] Eccentricity            [-]
%  i         [1x1] Inclination             [rad]
%  OM        [1x1] RAAN                    [rad]
%  om        [1x1] Pericentre anomaly      [rad]
%  theta     [1x1] True anomaly            [rad]
%  mu        [1x1] Gravitational parameter [km^3/s^2]
%
% OUTPUTS:
%
%  r         [3x1] Position vector         [km]
%  v         [3x1] Velocity vector         [km/s]
%
% First Edition: 01/07/2022
%
% Authors: Rocco Larocca, Matteo Mascelloni, Davide Lanza, Afaq Shakeel

a=k(1);
e=k(2);
i=k(3);
OM=k(4);
om=k(5);
theta=k(6);


p=a*(1-(e^2));
r_norm=p/(1+(e*cos(theta)));
r_PF=r_norm*[cos(theta); sin(theta); 0];
v_PF=(sqrt(mu/p))*[-sin(theta); (e+cos(theta)); 0];
R_3_OM=[cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];
R_1_i=[1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
R_3_om=[cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];
T_PF=(R_3_OM)'*(R_1_i)'*(R_3_om)';
r=T_PF*r_PF;
v=T_PF*v_PF;