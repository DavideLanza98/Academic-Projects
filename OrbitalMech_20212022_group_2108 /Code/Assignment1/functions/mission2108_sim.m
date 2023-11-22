function [dv,dv_dep,dv_arr]=mission2108_sim(dep,arr,pl_depID,pl_arrID)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mission2108_sim: simplified version of mission2108.m useful to estimate
%                  the total cost of a single Lambert leg
%
% INPUTS:
%
% dep      [1]    Departure date in days from 01/01/2000              [days]
% FB       [1]    Fly-By date in days from 01/01/2000                 [days]
% arr      [1]    Arrival date in days from 01/01/2000                [days]
% pl_depID [1]   ID value for departure planet                        [-]
% pl_arrID [1]   ID value for arrival planet                          [-]
%
% OUTPUTS:
%
% dv       [1]    Total cost of the interplanetary mission            [km/s]
% dv_dep   [1]    Cost of the first leg (dep to FB)                   [km/s]
% dv_PGA   [1]    Cost of the fly-by maneuvre                         [km/s]
% dv_arr   [1]    Cost of the second leg (FB to arr)                  [km/s]
% rp_s     [1]    Radius of hyperbola perigee                         [km]
% v_min    [1]    Incoming velocity magnitude at perigee              [km/s]
% v_plus   [1]    Outcoming velocity magnitude at perigee             [km/s]
% V_L1_i   [3x1]  Initial velocity vector in the first Lambert leg    [km/s]
% V_L1_f   [3x1]  Final velocity vector in the first Lambert leg      [km/s]
% V_L2_i   [3x1]  Initial velocity vector in the second Lambert leg   [km/s]
% V_L2_f   [3x1]  Final velocity vector in the second Lambert leg     [km/s]
% r_FB     [3x1]  Position vector of FB planet at Fly-By              [km]
% v_FB     [3x1]  Velocity vector of FB planet at Fly-BY              [km/s]
%
% First Edition: 01/07/2022
%
% Authors: Rocco Larocca, Matteo Mascelloni, Davide Lanza, Afaq Shakeel


dv=NaN;
dv_dep=NaN;
dv_arr=NaN;

if arr-dep>=0 
    pl_dep=pl_depID;
    pl_arr=pl_arrID;
    muS=astroConstants(4);
    
    kep_dep=uplanet(dep,pl_dep);
    [r_dep,v_dep]=kep2carmod(kep_dep,muS);
    kep_arr=uplanet(arr,pl_arr);
    [r_arr,v_arr]=kep2carmod(kep_arr,muS);
    
    options=1;
    TOF1=86400*(arr-dep);
    [~,~,~,ERROR1,V_L1_i,V_L1_f,~,~]=lambertMR(r_dep,r_arr,TOF1,muS,0,0,options);
    V_L1_i=V_L1_i';
    V_L1_f=V_L1_f';
    
    if ERROR1==0
        dv_dep=norm(V_L1_i-v_dep);
        dv_arr=norm(-V_L1_f+v_arr);
        dv=dv_dep+dv_arr;
    end
end
       