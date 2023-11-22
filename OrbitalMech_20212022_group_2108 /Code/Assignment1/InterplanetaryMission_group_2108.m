%% MAIN SCRIPT - ASSIGNIMENT 2: INTERPLANETARY EXPLORER MISSION
%
% CONTRIBUTORS:
%  Lanza Davide
%  Larocca Rocco
%  Mascelloni Matteo 
%  Shakeel Afaq
% 
% VERSION:
%  07/01/2022 
%

%% Data
clc; clear; close all

addpath 'functions'
addpath 'Assigned Functions'/

pl_dep=1;     % uplanet ID for Mercury
pl_FB=3;      % uplanet ID for Earth
pl_arr=6;     % uplanet ID for Saturn

muS=astroConstants(4);              % Gravitational Parameter of Sun
mu_dep=astroConstants(11);          % Gravitational Parameter of Mercury
mu_FB=astroConstants(13);           % Gravitational Parameter of Earth
mu_arr=astroConstants(16);          % Gravitational Parameter of Saturn
R_FB=6378.388;                      % Radius of FlyBy planet: Earth
h_atm=200;                          % Critical atmospheric heigth

dep_earliest=[2026 02 01 00 00 00];    % Earliest departure date
arr_latest=[2066 02 01 23 59 59];      % Latest arrival date
dep_min=date2mjd2000(dep_earliest);
arr_max=date2mjd2000(arr_latest);

[kep_Mercury,~] = uplanet(0,pl_dep);
T_Mercury=1/(86400)*2*pi*sqrt(kep_Mercury(1)^3/muS);
[r_Mercury,v_Mercury]=kep2carmod(kep_Mercury,muS);
[kep_Earth,~] = uplanet(0,pl_FB);
T_Earth=1/(86400)*2*pi*sqrt(kep_Earth(1)^3/muS);
[r_Earth,v_Earth]=kep2carmod(kep_Earth,muS);
[kep_Saturn,~] = uplanet(0,pl_arr);
T_Saturn=1/(86400)*2*pi*sqrt(kep_Saturn(1)^3/muS);
[r_Saturn,v_Saturn]=kep2carmod(kep_Saturn,muS);

T_synME=(T_Mercury*T_Earth)/abs(T_Mercury-T_Earth);  % Syndoic period of Earth WRT Mercury
T_synESA=(T_Earth*T_Saturn)/abs(T_Earth-T_Saturn);   % Synodic period of Saturn WRT Earth

c1=abs(norm(r_Mercury)-norm(r_Earth));    % Chord of the first parabola (dep to FB)
s1=(c1+norm(r_Mercury)+norm(r_Earth))/2;  % Semi-perimeter of the first space triangle

t_p1=sqrt(2)/3*sqrt(s1^3/muS)*(1-((s1-c1)/s1)^(3/2));

t_p1=t_p1/86400;       % Parabolic TPAR in days for the first leg

c2=abs(norm(r_Earth)-norm(r_Saturn));     % Chord of the second parabola (FB to arr)
s2 =(c2+norm(r_Earth)+norm(r_Saturn))/2;  % Semi-perimeter of the second space triangle 

t_p2=sqrt(2)/3*sqrt(s2^3/muS)*(1-((s2-c2)/s2)^(3/2));

t_p2=t_p2/86400;       % Parabolic TPAR in days for the second leg

% The synodic periods appear too short to be individuated as helpful to
% determine the departure window: let's impose the latest departure date as
% the sum of dep_min plus half of one Saturn revolution period T_Saturn

dep=[dep_min dep_min+0.5*T_Saturn];     % Departure window determination

% In order to estimate the FlyBy and arrival time windows, we approximate
% our trajectories to two Hohmann transfers, calculating the relative
% geometries and TOFs, as well as the transfer velocities

a_L1=(kep_Mercury(1)+kep_Earth(1))/2;    % Semi-major axis of first Hohmann leg
dt_L1=(1/86400)*pi*sqrt((a_L1^3)/muS);   % TOF of first Hohmann leg
a_L2=(kep_Saturn(1)+kep_Earth(1))/2;     % Semi-major axis of second Hohmann leg
dt_L2=(1/86400)*pi*sqrt((a_L2^3)/muS);   % TOF of second Hohmann leg
dvA=sqrt(muS/kep_Mercury(1))*(sqrt((2*kep_Earth(1))/(2*a_L1))-1);
dvB=sqrt(muS/kep_Earth(1))*(-sqrt((2*kep_Mercury(1))/(2*a_L1))+1);
dvC=sqrt(muS/kep_Earth(1))*(sqrt((2*kep_Saturn(1))/(2*a_L2))-1);
dvD=sqrt(muS/kep_Saturn(1))*(-sqrt((2*kep_Earth(1))/(2*a_L2))+1);

% In order to compute the Fly-By and arrival windows we consider as lower
% time limits the parabolic times for leg 1 and leg 2, while as upper time
% limits the 125% of the Hohmann times of flights.

FB=[dep(1)+t_p1 dep(2)+1.25*dt_L1];      % FlyBy time window
arr=[FB(1)+t_p2 FB(2)+1.25*dt_L2];       % Arrival time window
FB_earliest = mjd20002date(FB(1));
FB_latest = mjd20002date(FB(2));
arr_earliest=mjd20002date(arr(1));
arr_latest1=mjd20002date(arr(2));
dep_latest=mjd20002date(dep(2));
arr_latest=mjd20002date(arr(2));



%% Genetic Alhorithm "ga.m"

rng default                     % For reproducibility
lb=[dep(1) FB(1) arr(1)];       % Lower boundaries for genetic algorithm
ub=[dep(2) FB(2) arr(2)];       % Upper boundaries for genetic algorithm
N_it=5;                   % To check the algorithm convergence
t_ga=NaN(N_it,3);         % Matrix of possible times combinations
dvtot_ga=NaN(N_it,1);     % Vector of mission costs

tic
options=optimoptions('ga','PlotFcn',@gaplotbestf,'FunctionTolerance',1e-04,'PopulationSize',1000,'MaxGeneration',100,'Display','off','HybridFcn','patternsearch');

for i=1:N_it
    [t_ga(i,:),dvtot_ga(i),~,~,~]=ga(@(x) mission2108(x(1),x(2),x(3)),3,[],[],[],[],lb,ub,[],options);
end

comp_time_ga=toc;        % Computational time for genetic algorithm

[dvtot_ga_min,j]=min(dvtot_ga);   % Minimum cost calculated by genetic algorithm
t_ga_min=t_ga(j,:);               
t_ga_final=[mjd20002date(t_ga_min(1));mjd20002date(t_ga_min(2));mjd20002date(t_ga_min(3))]; 

fprintf('Genetic Algorithm Results:\n')
fprintf('Departure date:\n');
t_ga_final(1,:)
fprintf('FlyBy date:\n');
t_ga_final(2,:)
fprintf('Arrival date\n');
t_ga_final(3,:)
fprintf('Cost of the mission:\n')
dvtot_ga_min
fprintf('Computational time:\n')
comp_time_ga



%% Second optimization: grid search

% We have to define three new time windows using as "middle values" the
% t_ga_final obtained with the Genetic algorithm. We choose to add and
% subtract 365 days to each value in t_ga_min.
% The linspace interval has been choosen in order to consider two dates for
% month in dep, FB and arr windows.

dep_w=linspace(t_ga_min(1)-365, t_ga_min(1)+365,48);  
FB_w=linspace(t_ga_min(2)-365, t_ga_min(2)+365,48);
arr_w=linspace(t_ga_min(3)-365,t_ga_min(3)+365,48);

dvtot_gs=NaN(length(dep_w),length(FB_w),length(arr_w));

tic
for i = 1:length(dep_w)
    for j = 1:length(FB_w)
        for k = 1:length(arr_w)
            [dvtot_gs(i,j,k)] = mission2108(dep_w(i),FB_w(j),arr_w(k));
        end
    end
end
comp_time_gs=toc;   % Computational time of grid search algorithm

[dvtot_gsmin,place] = min(dvtot_gs(:));  % Minimum cost calculated by grid search algorithm
[iii,jjj,kkk] = ind2sub(size(dvtot_gs),place);
t_gs = [dep_w(iii) FB_w(jjj) arr_w(kkk)];  
t_gs_final=[mjd20002date(t_gs(1));mjd20002date(t_gs(2));mjd20002date(t_gs(3))];

fprintf('Grid Search Results:\n')
fprintf('Departure date:\n');
t_gs_final(1,:)
fprintf('FlyBy date:\n');
t_gs_final(2,:)
fprintf('Arrival date\n');
t_gs_final(3,:)
fprintf('Cost of the mission:\n')
dvtot_gsmin
fprintf('Computational time:\n')
comp_time_gs


%% Third optimization:

% As shown, the computational time of ga.m and grid search algorithms are too
% elevated: to go over this problem, we apply a fmincon algorithm searching
% for the minimum dvtot_f and for the associated dates. We also impose a
% inequality condition to specify that the dates must be such that FB>dep+TOF1,
% arr>FB+TOF2 and arr>dep+TOF1+TOF2. 

rng default                       % For reproducibility
lbf=[dep(1); FB(1); arr(1)];      % Lower boundaries for fmincon algorithm
ubf=[dep(2); FB(2); arr(2)];      % Upper boundaries for fmincon algorithm

A=[1 -1 0; 0 1 -1; 1 0 -1];           % Inequality matrix
b=[-dt_L1; -dt_L2; -(dt_L1+dt_L2)];   % Inequality vector
options=optimoptions('fmincon','OptimalityTolerance',1e-12,'StepTolerance',1e-12,'PlotFcn',@optimplotfval,'Display','off');
x0=[t_gs(1);t_gs(2);t_gs(3)];         % Initial conditions taken from grid search

tic
[t_f,dvtot_f]=fmincon(@(x) mission2108(x(1),x(2),x(3)),x0,A,b,[],[],lbf,ubf,[],options);
comp_time_f=toc;        % Computational time of fmincon algorithm

t_f_final=[mjd20002date(t_f(1));mjd20002date(t_f(2));mjd20002date(t_f(3))];

fprintf('Fmincon Algorithm Results:\n')
fprintf('Departure date:\n');
t_f_final(1,:)
fprintf('FlyBy date:\n');
t_f_final(2,:)
fprintf('Arrival date\n');
t_f_final(3,:)
fprintf('Cost of the mission:\n')
dvtot_f
fprintf('Computational time:\n')
comp_time_f

%% Final results

if dvtot_ga_min<=dvtot_gsmin && dvtot_ga_min<=dvtot_f
    t_s=t_ga_min;
    deltatot_s=dvtot_ga_min;
    fprintf('The genetic algorithm shows the minimum cost of the mission:\n')
    fprintf('Chosen dates:\n')
    fprintf('Departure date:\n');
    t_ga_final(1,:)
    fprintf('FlyBy date:\n');
    t_ga_final(2,:)
    fprintf('Arrival date\n');
    t_ga_final(3,:)
    fprintf('Cost of the mission:\n')
    dvtot_ga_min
elseif dvtot_gsmin<=dvtot_ga_min && dvtot_gsmin<=dvtot_f
    t_s=t_gs;
    deltatot_s=dvtot_gsmin;
    fprintf('The grid search algorithm shows the minimum cost of the mission:\n')
    fprintf('Chosen dates:\n')
    fprintf('Departure date:\n');
    t_gs_final(1,:)
    fprintf('FlyBy date:\n');
    t_gs_final(2,:)
    fprintf('Arrival date\n');
    t_gs_final(3,:)
    fprintf('Cost of the mission:\n')
    dvtot_gsmin
elseif dvtot_f<=dvtot_ga_min && dvtot_f<=dvtot_gsmin
    t_s=t_f;
    deltatot_s=dvtot_f;
    fprintf('The fmincon algorithm shows the minimum cost of the mission:\n')
    fprintf('Chosen dates:\n')
    fprintf('Departure date:\n');
    t_f_final(1,:)
    fprintf('FlyBy date:\n');
    t_f_final(2,:)
    fprintf('Arrival date\n');
    t_f_final(3,:)
    fprintf('Cost of the mission:\n')
    dvtot_f
end

%% Plot of the mission

[kep_dep,~]=uplanet(t_s(1),pl_dep);
[kep_FB,~]=uplanet(t_s(2),pl_FB);
[kep_arr,~]=uplanet(t_s(3),pl_arr);

[r_dep,v_dep]=kep2carmod(kep_dep,muS);
[r_FB,v_FB]=kep2carmod(kep_FB,muS);
[r_arr,v_arr]=kep2carmod(kep_arr,muS);

% Trajectory representation 

[kep_Mercury,~] = uplanet(t_s(1),pl_dep);
T_Mercurys=2*pi*sqrt(kep_Mercury(1)^3/muS);
[kep_Earth,~] = uplanet(t_s(2),pl_FB);
T_Earths=2*pi*sqrt(kep_Earth(1)^3/muS);
[kep_Saturn,~] = uplanet(t_s(3),pl_arr);
T_Saturns=2*pi*sqrt(kep_Saturn(1)^3/muS);

[~,dv_dep,~,dv_arr,~,~,~,V_L1_i,V_L1_f,V_L2_i,V_L2_f,~,~]=mission2108(t_s(1),t_s(2),t_s(3));

[a_L1_i,~,e_L1_i,i_L1_i,OM_L1_i,om_L1_i,th_L1_i]=car2kep(r_dep,V_L1_i,muS);  
[a_L1_f,~,e_L1_f,i_L1_f,OM_L1_f,om_L1_f,th_L1_f]=car2kep(r_FB,V_L1_f,muS);

[a_L2_i,~,e_L2_i,i_L2_i,OM_L2_i,om_L2_i,th_L2_i]=car2kep(r_FB,V_L2_i,muS);
[a_L2_f,~,e_L2_f,i_L2_f,OM_L2_f,om_L2_f,th_L2_f]=car2kep(r_arr,V_L2_f,muS);

tspan_dep=linspace(0,T_Mercurys,1000);
tspan_FB=linspace(0,T_Earths,1000);
tspan_arr=linspace(0,T_Saturns,1000);
tspan_1=linspace(0,86400*(t_s(2)-t_s(1)),1000);
tspan_2=linspace(0,86400*(t_s(3)-t_s(2)),1000);

options=odeset('RelTol',1e-13, 'AbsTol',1e-13);
[time_dep,state_dep]=ode113(@(t,s) tbp_ode(t,s,muS), tspan_dep, [r_dep;v_dep],options);
[time_FB,state_FB]=ode113(@(t,s) tbp_ode(t,s,muS), tspan_FB, [r_FB;v_FB],options);
[time_arr,state_arr]=ode113(@(t,s) tbp_ode(t,s,muS), tspan_arr, [r_arr;v_arr],options);
[time_1,state_1]=ode113(@(t,s) tbp_ode(t,s,muS), tspan_1, [r_dep;V_L1_i],options);
[time_2,state_2]=ode113(@(t,s) tbp_ode(t,s,muS), tspan_2, [r_FB;V_L2_i],options);

f=figure;

hold on
plot3(state_dep(:,1),state_dep(:,2),state_dep(:,3),'-.','LineWidth',2)
plot3(state_FB(:,1),state_FB(:,2),state_FB(:,3),'-.','LineWidth',2)
plot3(state_arr(:,1),state_arr(:,2),state_arr(:,3),'-.k','LineWidth',2)
plot3(state_1(:,1),state_1(:,2),state_1(:,3),'LineWidth',2)
plot3(state_2(:,1),state_2(:,2),state_2(:,3),'LineWidth',2)
pl3d(1,r_dep,f,6000)
pl3d(3,r_FB,f,5000)
pl3d(6,r_arr,f,1000)
pl3d(0,[0;0;0],f,25)
grid on
axis equal
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
title('Discovery One: Interplanetary Trajectory')
legend('Mercury orbit','Earth orbit','Saturn orbit','Leg 1: dep to FB','Leg 2: FB to arr')
hold off

%% Plot of the fly-by hyperbola

V_minf=V_L1_f-v_FB;
V_pinf=V_L2_i-v_FB;
R_critic=R_FB+h_atm;
[dv_PGA,rp,v_min,v_plus,a_PGA_min,a_PGA_plus,e_PGA_min,e_PGA_plus]=powGA(V_pinf,V_minf,R_critic,mu_FB);

tspan_hyp = 0:100:9000;
f=figure;
hold on
y01 = [[rp;0;0],[0;-v_min;0]]; 
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-13 );
[~, Y1] = ode113( @(t,y) tbp_ode(t,y,mu_FB), tspan_hyp, y01, options );
y02 = [[rp;0;0],[0;v_plus;0]]; 
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-13 );
[~, Y2] = ode113( @(t,y) tbp_ode(t,y,mu_FB), tspan_hyp, y02, options );
plot3(Y2(:,2),Y2(:,1),Y2(:,3), 'LineWidth', 2)
plot3(Y1(:,2),Y1(:,1),Y1(:,3), 'LineWidth', 2)
pl3d(3, [0 0 0], f, 1)
grid on
axis equal
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
title('Discovery One: Fly-By Hyperbola')
legend('Outcoming hyperbola','Incoming hyperbola')
hold off

%% Countour-Plots and Pork-Chop plots

int_dep_contour = NaN(1,length(dep_w));
int_FB_contour = NaN(1,length(FB_w));
int_arr_contour = NaN(1,length(arr_w));

for i = 1:length(dep_w)
    for j=1:length(FB_w)
        for k=1:length(arr_w)
         int_dep_contour(i) = datenum(mjd20002date(dep_w(i)));
         int_FB_contour(j) = datenum(mjd20002date(FB_w(j)));
         int_arr_contour(k)=datenum(mjd20002date(arr_w(k)));
        end
    end
end

dv_leg1=NaN(length(dep_w),length(FB_w));
dv_leg2=NaN(length(FB_w),length(arr_w));


for i=1:length(dep_w)
    for j=1:length(FB_w)
        [~,dv_leg1(i,j),~] = mission2108_sim(dep_w(i),FB_w(j),pl_dep,pl_FB);
    end
end

for j=1:length(FB_w)
    for k=1:length(arr_w)
        [~,~,dv_leg2(j,k)] = mission2108_sim(FB_w(j),arr_w(k),pl_FB,pl_arr);
    end
end

% Pork-Chop plot for Mercury-to-Earth transfer

hold on
plot(datenum(t_ga_final(1,:)),datenum(t_ga_final(2,:)),'kd','LineWidth',3)
[C_leg1,h_leg1]=contour(int_dep_contour,int_FB_contour,(dv_leg1)',min(min(dv_leg1))+(0:1:10));
caxis([min(min(dv_leg1)) min(min(dv_leg1))+10])
c = colorbar;
c.Label.String = 'Cost of the first maneuvre [km/s]';
clabel(C_leg1,h_leg1,min(min(dv_leg1))+(0:1:10))
datetick('x', 'yyyy mm dd', 'keeplimits','keepticks')
datetick('y', 'yyyy mm dd', 'keeplimits')
xtickangle(45)
ytickangle(45)
xlabel('Departure window')
ylabel('Fly-By window')
title('Pork-Chop plot for the dep-FB leg')
legend('Chosen dates')
hold off

% Pork-Chop plot for Earth-to-Saturn transfer

plot(datenum(t_ga_final(2,:)),datenum(t_ga_final(3,:)),'kd','LineWidth',3)
hold on
[C_leg2,h_leg2]=contour(int_FB_contour,int_arr_contour,(dv_leg2)',min(min(dv_leg2))+(0:1:10));
caxis([min(min(dv_leg2)) min(min(dv_leg2))+10])
c = colorbar;
c.Label.String = 'Cost of the second maneuvre [km/s]';
clabel(C_leg2,h_leg2,min(min(dv_leg2))+(0:1:10))
datetick('x', 'yyyy mm dd', 'keeplimits','keepticks')
datetick('y', 'yyyy mm dd', 'keeplimits')
xtickangle(45)
ytickangle(45)
xlabel('Fly-By window')
ylabel('Arrival window')
title('Pork-Chop plot for the FB-arr leg')
legend('Chosen dates')
hold off





