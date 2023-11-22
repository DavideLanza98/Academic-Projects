function [X,Y,Z] = plot3D(a,e,i,OM,om,thi,thf,dth,muP)
%
% plot3D.m - Computes the X,Y,Z coordinates by the keplerian elements in
% order to plot the orbit evolution
% 
% PROTOTYPE:
%  [X,Y,Z] = plot3D(a,e,i,OM,om,thi,thf,dth,muP)
%
% INPUT:
% a [1x1]   Semi-major axis        [km]
% e [1x1]   Eccentricity           [-]
% i [1x1]   Inclination            [rad]
% OM [1x1]  RAAN                   [rad]
% om [1x1]  Pericentre anomaly     [rad]
% thi [1x1] initial true anomaly   [rad]
% thf [1x1] final true anomaly     [rad]
%  muP      Gravitational parameter [km^3/s^2]
%
% OUTPUT:
%  X [N,1]      x coordinate      [km]
%  Y [N,1]      y coordinate      [km]
%  Z [N,1]      z coordinate      [km]

% CONTRIBUTORS:
%  Lanza Davide
%  Larocca Rocco
%  Mascelloni Matteo 
%  Shakeel Afaq
% 
% VERSION:
%  07-01-2022

% wrapping initial and final theta
thi=wrapTo2Pi(thi);
thf=wrapTo2Pi(thf);

if thi > thf
    thf = thf + 2*pi;
end

%building a theta vector
thvect = [thi:dth:thf];

rmat = zeros(length(thvect),3);

kep = [a e i OM om thi];

for j = 1:length(thvect)
    kep(6) = thvect(j);
    [rr,~] = kep2car(kep,muP);   
    rmat(j,:) = rr; 
end

% output
X = rmat(:,1);
Y = rmat(:,2);
Z = rmat(:,3);


end
