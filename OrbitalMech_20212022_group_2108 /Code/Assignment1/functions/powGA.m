function [dv_PGA,rp_s,v_min,v_plus,a_PGA_min,a_PGA_plus,e_PGA_min,e_PGA_plus]=powGA(V_pinf,V_minf,R_critic,muFB)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% powGA: function useful to estimate the fly-by hyperbola for an
%        interplanetary mission. 
%
% INPUTS:
% V_pinf     [3x1] Outcoming excess velocity vector        [km/s]
% V_minf     [3x1] Incoming excess velocity vector         [km/s]
% R_critic   [1]   Critical height from FB planet centre   [km]
% muFB       [1]   Gravitational parameter for FB planet   [km^3/s^2]
%
% OUTPUTS:
% dv_PGA     [1]   Cost of the fly-by mission              [km/s]
% rp_s       [1]   Radius of hyperbola perigee             [km]
% v_min      [1]   Incoming velocity magnitude at perigee  [km/s]
% v_plus     [1]   Outcoming velocity magnitude at perigee [km/s]
% a_PGA_min  [1]   Incoming hyperbola semi-major axis      [km]
% a_PGA_plus [1]   Outcoming hyperbola semi-major axis     [km]
% e_PGA_min  [1]   Incoming hyperbola eccentricity         [-]
% e_PGA_plus [1]   Outcoming hyperbola eccentricity        [-]
%
% First Edition: 01/07/2022
%
% Authors: Rocco Larocca, Matteo Mascelloni, Davide Lanza, Afaq Shakeel


delta_hint=acos((dot(V_pinf',V_minf'))/(norm(V_pinf)*norm(V_minf)));
rp_guess=6378.388;
options = optimset('Display','off');
fun=@(rp) asin(1/(1+rp/muFB*norm(V_minf)^2))+asin(1/(1+rp/muFB*norm(V_pinf)^2))-delta_hint;
rp_s=fsolve(fun,rp_guess,options);

dv_PGA=NaN;
v_min=NaN;
v_plus=NaN;
a_PGA_min=NaN;
a_PGA_plus=NaN;
e_PGA_min=NaN;
e_PGA_plus=NaN;

if rp_s>R_critic
    a_PGA_min=-muFB/norm(V_minf)^2;
    a_PGA_plus=-muFB/norm(V_pinf)^2;
    e_PGA_min=1+rp_s*norm(V_minf)^2/muFB;
    e_PGA_plus=1+rp_s*norm(V_pinf)^2/muFB;
    v_min=sqrt(2*muFB*(1./rp_s+1/2/abs(a_PGA_min)));
    v_plus=sqrt(2*muFB*(1./rp_s+1/2/abs(a_PGA_plus)));
    dv_PGA=abs(v_plus-v_min);
else
    rp_s=NaN;
end