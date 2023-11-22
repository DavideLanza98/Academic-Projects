function aJ2_RSW = kep_pert_RSW (kep,muP,RE,J2)
%
% kep_pert_RSW.m - computes pertubation accelaration due to J2 effect in
% the RSW frame. 
%
% PROTOTYPE:
%  aJ2_RSW = kep_pert_RSW (kep,muP,RE,J2)
% 
% INPUT:
%  kep [1x6]             Keplerian elements        
%  muP [1]               Earth's gravitational parameter            [km^3/s^2]
%  RE  [1]               Earth's radius                             [km]
%  J2  [1]               Coefficient for the second zonal harmonic  [-] 
% 
% OUTPUT:
%  aJ2_RSW [1x3]         Perturbed acceleration in RSW frame
% 
% CONTRIBUTORS:
%  Lanza Davide
%  Larocca Rocco
%  Mascelloni Matteo 
%  Shakeel Afaq
% 
% VERSION:
%  07-01-2022

% Keplerian parameters
a=kep(1);   e=kep(2);   i=kep(3);   
OM=kep(4);  om=kep(5);  th=kep(6);  


[r,~] = kep2car([a,e,i,OM,om,th],muP); 

r_norm=norm(r);    
w=th+om;            % argument of latitude


% Perturbed acceleration due to J2:
const = -(3/2)*(J2*muP*(RE^2))/(r_norm^4);

M = [ 1-(3*(sin(i)^2)*(sin(w)^2)); (sin(i)^2)*sin(2*w); sin(2*i)*sin(w)];
               
aJ2_RSW = const * M; 

end