function [a, e, i, OM, om, th] = car2kep(r, v, mu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% car2kep: function useful to convert from Cartesian coordinates to
%          Keplerian elements
%
% INPUTS:
%
% r    [3x1]  Position vector          [km]
% v    [3x1]  Velocity vector          [km/s]
% mu   [1x1]  Gravitational parameter  [km^3/s^2]
%
% OUTPUT:
%
% a    [1x1]  Semi-major axis          [Km]
% e    [1x1]  Eccentricity             [-]
% i    [1x1]  Inclination              [rad]
% OM   [1x1]  RAAN                     [rad]
% om   [1x1]  Pericentre anomaly       [rad]
% th   [1x1]  True anomaly             [rad]
%
% First Edition: 01/07/2022
%
% Authors: Rocco Larocca, Matteo Mascelloni, Davide Lanza, Afaq Shakeel

r_norm=norm(r);
v_norm=norm(v);

h=cross(r,v);
h_norm=norm(h);

i=acos((h(3))/(h_norm));

e_v=(1/mu).*((((v_norm)^2 -(mu/(r_norm))).*r)-((dot(r,v))*v));
e=norm(e_v);

Epsilon= 0.5*((v_norm)^2)-(mu/(r_norm));  
a=-(mu/(2*Epsilon));

K=[0;0;1];
N=cross(K,h);
N_norm=norm(N);

if N(2)>=0
    OM=acos(( N(1) )/(N_norm));
else
    OM=2*pi - acos(( N(1) )/(N_norm));
end

if e_v(3)>=0
    om=acos((dot(N,e_v) )/((N_norm)*(e)));
else
    om=2*pi - acos((dot(N,e_v) )/((N_norm)*(e)));
end

v_r= (dot(r,v))/(r_norm);

if v_r >=0
    th=acos((dot(e_v,r) )/((e)*(r_norm)));
else
    th=2*pi - acos((dot(e_v,r) )/((e)*(r_norm)));
end

end