function ds = kepl_orbit_J2_moon(s,t,muP,muMoon,J2,RE)
% 
% kepl_orbit_J2_moon.m - Computes the derivative of orbital state in 
% Cartesian coordinates considering J2 and Moon perturbing effects.
% 
% PROTOTYPE:
%  ds = kepl_orbit_J2_moon(x,t,muP,muMoon,J2,RE)
% 
% INPUT:
%  s [6,1]       Orbital state in cartesian coordinates                [km],[km/s]
%  t [1]         Time                                                  [s]
%  muP [1]       Gravitational parameter of the Earth                  [km^3/s^2]
%  muMoon [1]    Gravitational parameter of the Moon                   [km^3/s^2]
%  J2 [1]        Second zonal harmonic of the Earth   
%  RE [1]        Mean radius of the Earth                              [km]
% 
% OUTPUT:
%  ds [6,1]      Derivative of orbital state in cartesian coordinates  [km/s],[km/s^2]
% 
% CONTRIBUTORS:
%  Lanza Davide
%  Larocca Rocco
%  Mascelloni Matteo 
%  Shakeel Afaq
% 
% VERSION:
%  07-01-2022

r = norm(s(1:3));               
a_J2 = (3/2)*(J2/(r^4))*muP*(RE^2)*[(s(1)/r)*(5*((s(3)^2)/r^2)-1);(s(2)/r)*(5*((s(3)^2)/r^2)-1);(s(3)/r)*(5*((s(3)^2)/r^2)-3)]; 

t_days = t/86400;              
[r_moon, ~] = ephMoon(t_days);   % Moon position vector from ephemerides [km]
r_moon_norm=norm(r_moon);        % Norm of the Moon position vector [km]
r_s = s(1:3);                     % position vector of the s/c [km]
r_moon_s=r_moon'-r_s;              % position vector of the Moon wrt the s/c [km]
r_moon_s_norm=norm(r_moon_s);      % Norm of Moon position vector wrt the s/c [km]

% Perturbation acceleration due to the Moon's presence:
a_moon = muMoon .* ( (r_moon_s./(r_moon_s_norm^3) ) - (r_moon'./(r_moon_norm^3)) ); % Perturbing acceleration due to the Moon [km/s^2]

ds = [s(4);s(5);s(6);((-muP*s(1))/(r^3))+a_J2(1)+a_moon(1); ((-muP*s(2))/(r^3))+a_J2(2)+a_moon(2); ((-muP*s(3))/(r^3))+a_J2(3)+a_moon(3)];

end
