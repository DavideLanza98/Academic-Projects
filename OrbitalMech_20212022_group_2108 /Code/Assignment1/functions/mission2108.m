function [dv,dv_dep,dv_PGA,dv_arr,rp_s,v_min,v_plus,V_L1_i,V_L1_f,V_L2_i,V_L2_f,r_FB,v_FB]=mission2108(dep,FB,arr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mission2108: function useful to compute all the quantities needed to
%              estimate the total cost dv for the Discovery One interplanetary mission
%
% INPUTS:
%
% dep     [1]    Departure date in days from 01/01/2000             [days]
% FB      [1]    Fly-By date in days from 01/01/2000                [days]
% arr     [1]    Arrival date in days from 01/01/2000               [days]
%
% OUTPUTS:
%
% dv      [1]    Total cost of the interplanetary mission           [km/s]
% dv_dep  [1]    Cost of the first leg (dep to FB)                  [km/s]
% dv_PGA  [1]    Cost of the fly-by maneuvre                        [km/s]
% dv_arr  [1]    Cost of the second leg (FB to arr)                 [km/s]
% rp_s    [1]    Radius of hyperbola perigee                        [km]
% v_min   [1]    Incoming velocity magnitude at perigee             [km/s]
% v_plus  [1]    Outcoming velocity magnitude at perigee            [km/s]
% V_L1_i  [3x1]  Initial velocity vector in the first Lambert leg   [km/s]
% V_L1_f  [3x1]  Final velocity vector in the first Lambert leg     [km/s]
% V_L2_i  [3x1]  Initial velocity vector in the second Lambert leg  [km/s]
% V_L2_f  [3x1]  Final velocity vector in the second Lambert leg    [km/s]
% r_FB    [3x1]  Position vector of FB planet at Fly-By             [km]
% v_FB    [3x1]  Velocity vector of FB planet at Fly-BY             [km/s]
%
% First Edition: 01/07/2022
%
% Authors: Rocco Larocca, Matteo Mascelloni, Davide Lanza, Afaq Shakeel


dv=NaN;
dv_dep=NaN;
dv_PGA=NaN;
dv_arr=NaN;
rp_s=NaN;
v_min=NaN;
v_plus=NaN;
V_L1_i=NaN;
V_L1_f=NaN;
V_L2_i=NaN;
V_L2_f=NaN;
r_FB=NaN;
v_FB=NaN;

if FB-dep>=0 && arr-FB>=0
    pl_dep=1;
    pl_FB=3;
    pl_arr=6;
    muS=astroConstants(4);
    muFB=astroConstants(13);
    
    kep_dep=uplanet(dep,pl_dep);
    [r_dep,v_dep]=kep2carmod(kep_dep,muS);
    kep_FB=uplanet(FB,pl_FB);
    [r_FB,v_FB]=kep2carmod(kep_FB,muS);
    kep_arr=uplanet(arr,pl_arr);
    [r_arr,v_arr]=kep2carmod(kep_arr,muS);
    
    options=1;
    TOF1=86400*(FB-dep);
    [~,~,~,ERROR1,V_L1_i,V_L1_f,~,~]=lambertMR(r_dep,r_FB,TOF1,muS,0,0,options);
    V_L1_i=V_L1_i';
    V_L1_f=V_L1_f';
    TOF2=86400*(arr-FB);
    [~,~,~,ERROR2,V_L2_i,V_L2_f,~,~]=lambertMR(r_FB,r_arr,TOF2,muS,0,0,options);
    V_L2_i=V_L2_i';
    V_L2_f=V_L2_f';
    
    if ERROR1==0 && ERROR2==0
        V_minf=V_L1_f-v_FB;
        V_pinf=V_L2_i-v_FB;
        R_Earth=6378.388;
        h_atm=200;
        R_critic=R_Earth+h_atm;
        [dv_PGA,rp_s,v_min,v_plus]=powGA(V_pinf,V_minf,R_critic,muFB);
        if not(isnan(rp_s))
            dv_dep=norm(V_L1_i-v_dep);
            dv_arr=norm(V_L2_f-v_arr);
            dv=dv_dep+dv_arr+dv_PGA;
        end
    end
end
        
    
    
    