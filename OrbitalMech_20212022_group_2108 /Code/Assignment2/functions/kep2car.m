function [rr,vv] = kep2car(kep,mu)
% 
% kep2car.m - Computes the position and velocity vectors from keplerian
%             elements
% 
% PROTOTYPE:
%  [rr,vv] = kep2car(kep,mu)
% 
% INPUT:
%  kep [6]        Keplerian elements                                     [km, rad]
%  mu  [1]        Gravitational constant                                 [km^3/s^2]
% 
% OUTPUT:
%  rr [3,1]      Position vector                                        [km]
%  vv [3,1]      Velocity vector                                        [km/s]
% 
% CONTRIBUTORS:
%  Lanza Davide
%  Larocca Rocco
%  Mascelloni Matteo 
%  Shakeel Afaq
% 
% VERSION:
%  07-01-2022

% semi-latus reptum
p = kep(1)*(1-kep(2)^2);

% position 
r = p/(1+kep(2)*cos(kep(6)));

rr = r*[cos(kep(6));sin(kep(6));0];
vv = sqrt(mu/p)*[-sin(kep(6));kep(2)+cos(kep(6));0];

% transformation matrixes
R_OM = [ cos(kep(4)) sin(kep(4)) 0; -sin(kep(4)) cos(kep(4)) 0; 0 0 1];
R_i = [ 1 0 0; 0 cos(kep(3)) sin(kep(3)); 0 -sin(kep(3)) cos(kep(3))];
R_om = [ cos(kep(5)) sin(kep(5)) 0; -sin(kep(5)) cos(kep(5)) 0; 0 0 1];

T = R_om*R_i*R_OM;

% position and velocity vectors [3x1]
rr = T'*rr;
vv = T'*vv;
