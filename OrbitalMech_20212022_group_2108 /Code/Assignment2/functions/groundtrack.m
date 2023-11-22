function [lat, lon] = groundtrack (a, rr0,vv0, omE, muP, t_span,thetaG0 ,t0, perturbed,kep)
%
% groundtrack.m - Computes the latitude and longitude vector of an orbit during a given time period.
% 
% PROTOTYPE:
%  [lat, lon] = groundtrack (a, rr0,vv0, omE, muP, t_span,thetaG0 ,t0, perturbed, kep)
%
% INPUT:
%  a                                      semi-major axis                                [km]
%  rr0  [3x1]                             initial 3x1 position vector                    [km]
%  vv0  [3x1]                             initial 3x1 velocity vector                    [km/s]
%  omE                                    Earth's spin rate [rad/s]
%  muP                                    gravitational parameter of primary body        [km^3/s^2]
%  n_orb                                  number of orbits
%  thetaG0                                initial Greenwich siderial time                [rad]
%  t0                                     initial time
%  perturbed                              0 = non-perturbed motion ; 
%                                         1 = perturbed motion
% OUTPUT: 
%  lat [N,1]                              column vector of latitude at each time-step    [deg]
%  lon [N,1]                              column vector of longitude at each time-step   [deg]
%
% CONTRIBUTORS:
%  Lanza Davide
%  Larocca Rocco
%  Mascelloni Matteo 
%  Shakeel Afaq
% 
% VERSION:
%  07-01-2022

%% integration
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 ); %ODE solver tolerance
%ODE solver with Nx6 vector output (rrvv) providing x,y,z position coordinates and their velocity components at each time step:
    switch perturbed
        case 0
            [~, rrvv] = ode113( @(t,y) tbp_ode(t,y, muP), t_span, [rr0; vv0], options ); %non-perturbed motion
        case 1
            [~, rrvv] = ode113( @(t,y) tbp_ode_perturbed(t,y, muP,kep), t_span, [rr0; vv0], options ); %perturbed motion
    end

%% alpha, delta, latitude, longitude
[alpha,delta] = car2AlphaDelta (rrvv(:,1:3)); %right ascension alpha and declination delta vectors at each time step
 
thetaG = rad2deg(thetaG0 + omE * (t_span-t0)); %vector of Greenwich siderial time [deg]
alpha = rad2deg(alpha); %vector of right-ascension [deg]

lat = rad2deg(delta); %latitudes vector [deg]
lon = wrapTo180(alpha-thetaG'); %longitudes vector [deg] wrapped to [-180 180]

for i=2:length(lon)
    if abs(lon(i)-lon(i-1))>180 
        lon(i) = NaN;
    end

end

