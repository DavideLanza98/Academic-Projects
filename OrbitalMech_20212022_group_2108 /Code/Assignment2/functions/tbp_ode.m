function ds=tbp_ode(~,s,muP)

%
% tbp_ode: Differential equation for the non-perturbed two body problem 
%
% INPUTS:
%  s    [6x1]   State of the system              [km,km/s]
%  muP  [1]     Gravitational parameter          [km^3/s^2]
%
% OUTPUTS:
%  ds   [1]     State of the system derivative   [km/s^2,km/s^3]
%
% First Edition: 01/07/2022
%
% Authors: Rocco Larocca, Matteo Mascelloni, Davide Lanza, Afaq Shakeel

r=s(1:3); % [km]
v=s(4:6); % [km/s]

r_norm=norm(r);
v_norm=norm(v);

ds= [v
    (-muP/(r_norm^3))*r];
    
   
end
    