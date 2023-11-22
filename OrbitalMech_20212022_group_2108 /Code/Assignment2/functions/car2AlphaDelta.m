function [alpha,delta] = car2AlphaDelta (xyz) 
%
% car2AlphaDelta.m - Computes the right ascension alpha and the declination delta of a S/C at each time step
%
% PROTOTYPE:
%  [alpha,delta] = car2AlphaDelta (xyz)
% INPUT:
%  xyz [Nx3]         vector of x,y,z position coordinates at each time step     [km]
%
% OUTPUT: 
%  alpha            column vector of right-ascension at each time-step          [rad]
%  delta            column vector of declination at each time-step              [rad]
%
% CONTRIBUTORS:
%  Lanza Davide
%  Larocca Rocco
%  Mascelloni Matteo 
%  Shakeel Afaq
% 
% VERSION:
%  07-01-2022

%% position coordinates and norm vector
x_vect = xyz(:,1);  y_vect = xyz(:,2);  z_vect = xyz(:,3); % creating three different vectors
r_vect = vecnorm(xyz'); % row vector containing norm of x,y,z coordinates at each time step

%% alpha, delta 
delta = asin(z_vect./r_vect'); %declination (latitude) [rad]

k = y_vect./r_vect'; 
alpha = acos((x_vect./r_vect')./cos(delta)); %right ascension [rad]

    for i = 1:1:size(k)
        if k(i) <= 0
            alpha(i) = 2*pi - alpha(i); 
        end
    end

end

