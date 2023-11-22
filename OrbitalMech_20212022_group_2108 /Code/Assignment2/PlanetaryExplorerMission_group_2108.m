%% MAIN SCRIPT - ASSIGNIMENT 2: PLANETARY EXPLORER MISSION
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
% INSTRUCTIONS: 
%  In ORBIT PROPAGATION (PERTURBED) section remember to run from switch 1 
%  to 3 in order to avoid errors.

%% DATA 
clc; clear; close all

addpath Provided-scripts; 
addpath functions;

% Nominal orbit
a = 2.6039e4;           % semi-major axis [km]
e = 0.5792;             % eccentricity [-]
i = deg2rad(26.1613);   % inclination [rad]          
OM = 0;                 % right ascension of ascending node [deg] 
om = 0;                 % argument of perigee [deg]
th = 0;                 % initial true anomaly [deg]
h = 4586.201;           % altitude from Earth [km]

% Constants
muP = astroConstants(13);       % gravitational parameter of the Earth [km^3/s^2]
muMoon = astroConstants(20);    % gravitational parameter of the Moon [km^3/s^2]
J2 = astroConstants(9);         % second zonal harmonic of the Earth
RE = astroConstants(23);        % mean equatorial radius of the Earth [km]
omE = 2*pi/(23*3600+56*60);   % Earth angular velocity [rad/s]

% Repeating GT ratio
k = 2;  % S/C revoultions
m = 1;  % Earth revoulution

% semi-major axis for the repeated unperturbed case
a_rep = ((muP*m^2)/(omE*k)^2)^(1/3); % modified semi-major axis to obtain a repeating ground track [km]

c = (-3/2)*((sqrt(muP)*J2*(RE^2))/((1-e^2)^2)); %constant coefficient
% semi-major axis for the repeated perturbed case
% non linear equation with a_rep as initial guess
a_rep_pert = fzero(@(a) (m/k) - ((omE - (c/a^(7/2))*cos(i))/((sqrt(muP/a^3))+(c/a^(7/2))*((5/2)*(sin(i))^2-2)-((c*sqrt(1-e^2))/a^(7/2))*(1-(3/2)*(sin(i))^2))),a_rep);

% Keplerin elements vectors
kep = [a e i OM om th];         % Keplerian elements vector [1x6]
kep0 = kep;                     % Initial condition for propagation
kep_rep = [a_rep e i OM om th]; % Keplerian modified elements vector [1x6]
kep_rep_pert = [a_rep_pert e i OM om th]; % keplerian elements in repetuted and perturbed problem [1x6]

% Orbital periods
T = 2*pi*sqrt((a^3)/muP);     % orbital period [s]
T_rep = 2*pi*sqrt((a_rep^3)/muP);    % modified orbital period [s]
T_rep_pert = 2*pi*sqrt((a_rep_pert^3)/muP); % modified orbital period in rep and perturbed problem [s]

% Initial conditions 
t0 = 0;         % initial time [s]
thetaG0 = 0;    % right ascension of Greenwich at t0 [rad]

% time span vectors
step = 10; 
tspan_1orbit = t0:step:T;            % 1 orbit [s]
tspan_1day = t0:step:3600*24;      % 1 day [s]
tspan_10days = t0:step:3600*24*10;   % 10 days [s]
tspan_rep = t0:step:T_rep;     % 1 modified orbit [s]   
t_span_rep_pert = t0:step:T_rep_pert; % 1 rep and perturbed orbit [s]

date_in = [2020,12,01,00,00,00]; % real satellite start date [date]
t0_prop = date2mjd2000(date_in) * 86400; % initial date for propagation [s]

n_orbs = 300; % number of orbits 
tspan_prop = linspace(t0_prop, t0_prop + n_orbs * T, 10000);

options = odeset('RelTol',1e-13,'AbsTol',1e-14); % ode options for integrations

%% GROUND TRACK SELECTION:

% 1: non repeated and unperturbed ground tracks
% 2: repeated and unperturbed ground tracks
% 3: non repeated and perturbed ground tracks
% 4: repeated and perturbed ground tracks

GT_selection = 1; % Choose one of the previous options

if GT_selection ~= 1 & GT_selection ~= 2 & GT_selection ~= 3 & GT_selection ~= 4
    error('ERROR: GT_selection must be a number between 1 and 4!');
end

switch GT_selection

    case 1 % Non repeated and unperturbed ground track

        [rr0, vv0] = kep2car(kep, muP); % initial 3x1 position and velocity vectors 
        
        % Ground track of the nominal orbit over 1 orbit
        [lat, lon] = groundtrack (a, rr0, vv0, omE, muP, tspan_1orbit, thetaG0, t0, 0, kep); %latitude and longitude colummn vectors
        plot_Earth_2D();
        hold on;
        plot(lon,lat,'g','LineWidth',2)
        plot(lon(1),lat(1),'s','MarkerSize',10,'LineWidth',2); % starting point
        plot(lon(end),lat(end), 'x','MarkerSize',10,'LineWidth',2); % ending point
        title('1 ORBIT non repeated and unperturbed GT');
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');
        legend('Ground track', 'start', 'end');

        % Ground track of the nominal orbit over 1 day
        [lat, lon] = groundtrack (a, rr0,vv0, omE, muP, tspan_1day, thetaG0, t0, 0,kep); %latitude and longitude colummn vectors
        plot_Earth_2D();
        hold on;
        plot(lon,lat,'g','LineWidth',2)
        plot(lon(1),lat(1),'s','MarkerSize',10,'LineWidth',2); % starting point
        plot(lon(end),lat(end), 'x','MarkerSize',10,'LineWidth',2); % ending point
        title('1 DAY non repeated and unperturbed GT');
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');
        legend('Ground track', 'start', 'end')
        
        % Ground track of the nominal orbit over 10 days
        [lat, lon] = groundtrack (a, rr0,vv0, omE, muP, tspan_10days, thetaG0, t0, 0,kep); %latitude and longitude colummn vectors
        plot_Earth_2D();
        hold on;
        plot(lon,lat,'g','LineWidth',2)
        plot(lon(1),lat(1),'s','MarkerSize',10,'LineWidth',2); % starting point
        plot(lon(end),lat(end), 'x','MarkerSize',10,'LineWidth',2); % ending point
        title('10 DAY non repeated and unperturbed GT');
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');
        legend('Ground track', 'start', 'end')
        
    case 2 % Repeated and unperturbed ground tracks
        
        [rr0, vv0] = kep2car(kep_rep, muP); % initial 3x1 position and velocity vectors 
        
        % Repeating ground track over 1 modified orbit
        [lat, lon] = groundtrack (a_rep, rr0,vv0, omE, muP, tspan_rep, thetaG0, t0, 0,kep_rep); % latitude and longitude colummn vectors
        plot_Earth_2D();
        plot(lon,lat,'g','LineWidth',2)
        plot(lon(1),lat(1),'s','MarkerSize',10,'LineWidth',2); % starting point
        plot(lon(end),lat(end), 'x','MarkerSize',10,'LineWidth',2); % ending point
        title('1 ORBIT repeated and unperturbed GT');
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');
        legend('Ground track', 'start', 'end')
        
        % Repeating ground track over 1 day
        [lat, lon] = groundtrack (a_rep, rr0,vv0, omE, muP, tspan_1day, thetaG0, t0, 0,kep_rep); 
        plot_Earth_2D();
        hold on;
        plot(lon,lat,'g','LineWidth',2)
        plot(lon(1),lat(1),'s','MarkerSize',10,'LineWidth',2); % starting point
        plot(lon(end),lat(end), 'x','MarkerSize',10,'LineWidth',2); % ending point
        title('1 DAY repeated and unperturbed GT');
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');
        legend('Ground track', 'start', 'end')
        
        % Repeating ground track over 10 days
        [lat, lon] = groundtrack (a_rep, rr0,vv0, omE, muP, tspan_10days, thetaG0, t0, 0,kep_rep); % latitude and longitude colummn vectors
        plot_Earth_2D();
        hold on;
        plot(lon,lat,'g','LineWidth',2)
        plot(lon(1),lat(1),'s','MarkerSize',10,'LineWidth',2); % starting point
        plot(lon(end),lat(end), 'x','MarkerSize',10,'LineWidth',2); % ending point
        title('10 DAY repeated and unperturbed GT');
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');
        legend('Ground track', 'start', 'end')

    case 3 % Non repeated and perturbed ground tracks

        [rr0, vv0] = kep2car(kep, muP); % initial 3x1 position and velocity vectors 
        
        % Non repeated and perturbed ground track over 1 orbit
        [lat, lon] = groundtrack (a, rr0, vv0, omE, muP, tspan_1orbit, thetaG0, t0, 1,kep); % latitude and longitude colummn vectors
        plot_Earth_2D();
        hold on;
        plot(lon,lat,'g','LineWidth',2)
        plot(lon(1),lat(1),'s','MarkerSize',10,'LineWidth',2); % starting point
        plot(lon(end),lat(end), 'x','MarkerSize',10,'LineWidth',2); % ending point
        title('1 ORBIT non repeated and perturbed GT');
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');
        legend('Ground track', 'start', 'end')
        
        % Non repeated and perturbed ground track over 1 day
        [lat, lon] = groundtrack (a, rr0, vv0, omE, muP, tspan_1day, thetaG0, t0, 1,kep); % latitude and longitude colummn vectors
        plot_Earth_2D();
        plot(lon,lat,'g','LineWidth',2)
        plot(lon(1),lat(1),'s','MarkerSize',10,'LineWidth',2); % starting point
        plot(lon(end),lat(end), 'x','MarkerSize',10,'LineWidth',2); % ending point
        title('1 DAY non repeated and perturbed GT');
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');
        legend('Ground track', 'start', 'end')
        
        % Non repeated and perturbed ground track over 10 days
        [lat, lon] = groundtrack (a, rr0, vv0, omE, muP, tspan_10days, thetaG0, t0, 1,kep); % latitude and longitude colummn vectors 
        plot_Earth_2D();
        hold on;
        plot(lon,lat,'g','LineWidth',2)
        plot(lon(1),lat(1),'s','MarkerSize',10,'LineWidth',2); % starting point
        plot(lon(end),lat(end), 'x','MarkerSize',10,'LineWidth',2); % ending point
        title('10 DAYS non repeated and perturbed GT');
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');
        legend('Ground track', 'start', 'end')

    case 4 % Repeated and perturbed ground tracks
        
        [rr0, vv0] = kep2car(kep_rep_pert, muP); % initial [3x1] position and velocity vectors

        % Repeated and perturbed ground track over 1 modified orbit
        [lat, lon] = groundtrack (a_rep_pert, rr0, vv0, omE, muP, tspan_rep, thetaG0, t0, 1,kep_rep_pert); % latitude and longitude colummn vectors
        plot_Earth_2D();
        hold on;
        plot(lon,lat,'g','LineWidth',2)
        plot(lon(1),lat(1),'s','MarkerSize',10,'LineWidth',2); % starting point
        plot(lon(end),lat(end), 'x','MarkerSize',10,'LineWidth',2); % ending point
        title('1 ORBIT repeated and perturbed GT');
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');
        legend('Ground track', 'start', 'end')

        % Repeated and perturbed ground track over 1 day
        [lat, lon] = groundtrack (a_rep_pert, rr0, vv0, omE, muP, tspan_1day, thetaG0, t0, 1,kep_rep_pert); % latitude and longitude colummn vectors
        plot_Earth_2D();
        hold on;
        plot(lon,lat,'g','LineWidth',2)
        plot(lon(1),lat(1),'s','MarkerSize',10,'LineWidth',2); % starting point
        plot(lon(end),lat(end), 'x','MarkerSize',10,'LineWidth',2); % ending point
        title('1 DAY repeated and perturbed GT');
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');
        legend('Ground track', 'start', 'end')
        
        % Repeated and perturbed ground track over 10 days
        [lat, lon] = groundtrack (a_rep_pert, rr0, vv0, omE, muP, tspan_10days, thetaG0, t0, 1,kep_rep_pert); % latitude and longitude colummn vectors
        plot_Earth_2D();
        plot(lon,lat,'g','LineWidth',2)
        plot(lon(1),lat(1),'s','MarkerSize',10,'LineWidth',2); % starting point
        plot(lon(end),lat(end), 'x','MarkerSize',10,'LineWidth',2); % ending point
        title('10 DAYS repeated and perturbed GT');
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');
        legend('Ground track', 'start', 'end')

end

%%  ORBIT PROPAGATION (PERTURBED): Gauss and Cartesian method 

% METHOD SELECTION and ERRORS:
% INSTRUCTIONS: run 1 and 2 before running 3!!! ******IMPORTANT******
% 1: Gauss propagation
% 2: Cartesian propagation
% 3: Errors between Guass and Cartesian propagation

prop_method = 1; % Choose one of the previous options

if prop_method ~= 1 & prop_method ~= 2 & prop_method ~= 3 
    error('ERROR: pert_method must be a number between 1 and 3!')
end

switch prop_method

    case 1 % Gauss propagation

        Gauss_eqs = @(t,kep) ode_gaussEsq_RSW_J2_MOON (t,kep,muP,kep_pert_RSW (kep,muP,RE,J2),kepMoon_pert_RSW (t,kep,muMoon,muP)); 
        [T_Gauss,kep_gauss] = ode113(Gauss_eqs,tspan_prop,kep0,options); % orbit Gauss propagation
        
        kep_gauss(:,3:6) = rad2deg(kep_gauss(:,3:6)); % converting keplerian elements in deg from radians
        
        % plotting the six graphs in a figure as subplots
        figure()
        sgtitle('Time evolution with Gauss propagation')
        
        subplot(2,3,1)
        plot((T_Gauss-t0_prop)/T, (kep_gauss(:,1)), 'LineWidth', 2);
        xlabel('time [T]');
        ylabel('a [km]');
        title('a_{Gauss}');
        
        subplot(2,3,2)
        plot((T_Gauss-t0_prop)/T, (kep_gauss(:,2)), 'LineWidth', 2);
        xlabel('time [T]');
        ylabel('e [-]');
        title('e_{Gauss}');
        
        subplot(2,3,3)
        plot((T_Gauss-t0_prop)/T, (kep_gauss(:,3)),'LineWidth', 2);
        xlabel('time [T]');
        ylabel('i [deg]');
        title('i_{Gauss}');
        
        subplot(2,3,4)
        plot((T_Gauss-t0_prop)/T, (kep_gauss(:,4)),'LineWidth', 2);
        xlabel('time [T]');
        ylabel('\Omega [deg]');
        title('\Omega_{Gauss}');
        
        subplot(2,3,5)
        plot((T_Gauss-t0_prop)/T, (kep_gauss(:,5)),'LineWidth', 2);
        xlabel('time [T]');
        ylabel('\omega [deg]');
        title('\omega_{Gauss}');
        
        subplot(2,3,6)
        plot((T_Gauss-t0_prop)/T, ((kep_gauss(:,6))),'LineWidth', 2);
        xlabel('time [T]');
        ylabel('\theta [deg]');
        title('\theta_{Gauss}');

    case 2 % Cartesian propagation

        [rr,vv] = kep2car(kep0,muP); % initial position and velocity [3x1] vectors
        y0 = [rr(1:3); vv(1:3)]; % initial conditions
        
        kepl_eqs = @(t,s) kepl_orbit_J2_moon(s,t,muP,muMoon,J2,RE); 
        [T_Car,x_car] = ode113(kepl_eqs, tspan_prop,y0,options); % cartesian propagation
        
        % building the kep_car vector
        
        for k = 1:size(x_car,1)
            [a,e,i,OM,om,th] = car2kep(x_car(k,1:3),x_car(k,4:6), muP);
             kep_car(k,1) = a;
             kep_car(k,2) = e;
             kep_car(k,3) = i;
             kep_car(k,4) = OM;
             kep_car(k,5) = om;
             kep_car(k,6) = th;
        end
        
        kep_car(:,3:6) = unwrap(kep_car(:,3:6)); % unwrapping keplerian angles
        kep_car(:,3:6) = rad2deg(kep_car(:,3:6)); % converting from rad to deg
        
        % plotting the six graphs in a figure as subplots
        figure()
        sgtitle('Time evolution with Cartesian propagation')
        
        subplot(2,3,1)
        plot((T_Car-t0_prop)/T, (kep_car(:,1)),'LineWidth', 2);
        xlabel('time [T]');
        ylabel('a [km]');
        title('a_{car}');
        
        subplot(2,3,2)
        plot((T_Car-t0_prop)/T, (kep_car(:,2)),'LineWidth', 2);
        xlabel('time [T]');
        ylabel('e [-]');
        title('e_{car}');
        
        subplot(2,3,3)
        plot((T_Car-t0_prop)/T, (kep_car(:,3)),'LineWidth', 2);
        xlabel('time [T]');
        ylabel('i [deg]');
        title('i_{car}');
        
        subplot(2,3,4)
        plot((T_Car-t0_prop)/T, (kep_car(:,4)),'LineWidth', 2);
        xlabel('time [T]');
        ylabel('\Omega [deg]');
        title('\Omega_{car}');
        
        subplot(2,3,5)
        plot((T_Car-t0_prop)/T, (kep_car(:,5)),'LineWidth', 2);
        xlabel('time [T]');
        ylabel('\omega [deg]');
        title('\omega_{car}');
        
        subplot(2,3,6)
        plot((T_Car-t0_prop)/T, (kep_car(:,6)),'LineWidth', 2);
        xlabel('time [T]');
        ylabel('\theta [deg]');
        title('\theta_{car}');

    case 3 % Errors between Guass and Cartesian propagation
        
        % computing errors between Cartesian and Gauss methods
        err_a_plot = abs(kep_car(:,1)-kep_gauss(:,1))./abs(kep0(1));            % relative error (a)   
        err_e_plot = abs(kep_car(:,2)-kep_gauss(:,2));                          % absolute error (e)                              
        err_i_plot = abs(kep_car(:,3)-kep_gauss(:,3))/(360);                    % relative error (i)
        err_OM_plot = abs(kep_car(:,4)-kep_gauss(:,4))/(360);                   % relative error (OM)
        err_om_plot = abs(kep_car(:,5)-kep_gauss(:,5))/(360);                   % relative error (om)
        err_th_plot = abs(kep_car(:,6)-kep_gauss(:,6))./(abs(kep_gauss(:,6)));  % relative error (th)
        
        % plotting the errors
        figure();
        sgtitle('Relative errors between Gauss and Cartesian')
        
        subplot(2,3,1)
        semilogy((T_Gauss - t0_prop)/T,err_a_plot);
        xlabel('time [T]');
        ylabel('|a_{car} - a_{gauss}|/|a_0| ');
        title('err_{a}');
        
        subplot(2,3,2)
        semilogy((T_Gauss - t0_prop)/T,err_e_plot);
        xlabel('time [T]');
        ylabel('|e_{car} - e_{gauss}|');
        title('err_{e}');
        
        subplot(2,3,3)
        semilogy((T_Gauss - t0_prop)/T, err_i_plot);
        xlabel('time [T]');
        ylabel('|i_{car} - i_{gauss}|/360°');
        title('err_{i}')
        
        subplot(2,3,4)
        semilogy((T_Gauss - t0_prop)/T,err_OM_plot);
        xlabel('time [T]');
        ylabel('|\Omega_{car} - \Omega_{gauss}|/360°');
        title('err_{\Omega}')
        
        subplot(2,3,5)
        semilogy((T_Gauss - t0_prop)/T,err_om_plot)
        xlabel('time [T]');
        ylabel('\omega_{car} - \omega_{gauss}|/360°');
        title('err_{\omega}')
        
        subplot(2,3,6)
        semilogy((T_Gauss - t0_prop)/T,err_th_plot)
        xlabel('time [T]');
        ylabel('\theta_{car} - \theta_{gauss}|/360°');
        title('err_{\theta}')
end

%% Computaional efficiency comparison

[rr,vv] = kep2car(kep0,muP);
y0 = [rr(1:3); vv(1:3)]; 

n_orbs_comp = 50; % number of orbits
step_vec = [2000:300:5000]; % time step vector for each orbit

Gauss_eqs = @(t,kep) ode_gaussEsq_RSW_J2_MOON (t,kep,muP,kep_pert_RSW (kep,muP,RE,J2),kepMoon_pert_RSW (t,kep,muMoon,muP));
kepl_eqs = @(t,s) kepl_orbit_J2_moon(s,t,muP,muMoon,J2,RE); % Gauss propagation

% allocating space
comput_gauss = [];
comput_cart = [];

% computing cpu computational cost for each time step 
for k = 1:length(step_vec)
    tspan_comp = linspace(0, n_orbs_comp* T, n_orbs_comp*step_vec(k));
    
    t1 = cputime;
    [T_Gauss,kep_gauss] = ode113(Gauss_eqs,tspan_comp,kep0,options);
    t2 = cputime;
    comput_gauss = [comput_gauss; (t2-t1)];
        
    t1 = cputime;
    [T_Car,x_car] = ode113(kepl_eqs, tspan_comp,y0,options);
    t2 = cputime;
    comput_cart = [comput_cart; (t2-t1)];
end

% plotting results
plot(step_vec, comput_gauss, 'r','LineWidth',1.5)
hold on
grid on
plot(step_vec, comput_cart, 'b','LineWidth',1.5)
xlabel('n_{steps} per orbit')
ylabel('t_{CPU} [s]')
legend('Gauss ', 'Cartesian ');

%% ORBIT EVOLUTION REPRESENTATION 

N_points = 50;                                              % number of points for the timespan discretization
t_span_evol=linspace(0,365*2*86400,N_points);                          % timespan for propagation [s]
Gauss_eqs= @(t,kep) ode_gaussEsq_RSW_J2_MOON (t,kep,muP,kep_pert_RSW (kep,muP,RE,J2),kepMoon_pert_RSW (t,kep,muMoon,muP));
[T_Gauss,kep_gauss] = ode113(Gauss_eqs,t_span_evol,kep0,options); % propagation of the orbit through Gauss equations

% plots and settings
f=figure;
hold on
pl3d(3, [0,0,0], f, 1);
colmap=winter(length(T_Gauss));
set(gcf,'DefaultAxesColorOrder',colmap);
colormap(colmap);
hcb=colorbar;
set(get(hcb,'Title'),'String','Time in a year [days]');
caxis([0 (t_span_evol(end)/86400)]);

for k=1:N_points
    plotcolour=colmap(k,:);
    [X,Y,Z]=plot3D(kep_gauss(k,1),kep_gauss(k,2),kep_gauss(k,3),kep_gauss(k,4),kep_gauss(k,5),0,2*pi,0.01,muP);
    plot3(X,Y,Z, 'HandleVisibility', 'off', 'LineStyle', '-', 'LineWidth', 1.4,  'Color', plotcolour);
    xlabel('X [km]')
    ylabel('Y [km]')
    zlabel('Z [km]')
    axis equal
    title('Evolution of the orbit in a year')
end
hold off

%% FILTERING OF HIGH FREQUENCIES 

Gauss_eqs = @(t,kep) ode_gaussEsq_RSW_J2_MOON (t,kep,muP,kep_pert_RSW (kep,muP,RE,J2),kepMoon_pert_RSW (t,kep,muMoon,muP));
[T_Gauss,kep_gauss] = ode113(Gauss_eqs,tspan_prop,kep0,options); % Gauss propagation

kep_gauss(:,3:6) = rad2deg(kep_gauss(:,3:6)); % keplerian elements [deg]

window_size = 1:2:5000; % vector of different window sizes
SAD = []; % initialize Sum of Absolute Difference vector
for j = 1:6
    for l = 1:length(window_size)
        smoothedSignal(:,j) = movmean(kep_gauss(:,j),window_size(l)); % filtered signal for each windows size
        SAD(l,j) = sum(abs(smoothedSignal(:,j) - kep_gauss(:,j))); % Sum of Absolute Difference vector
    end
end

% Subplotting the original signal and SAD graph
keplabel = {'a [km]', 'e', 'i [deg]', '\Omega [deg]', '\omega [deg]', '\theta [deg]'}; % string vector with labels

for j = 1:6
    figure(j)
    subplot(2,1,1);
    plot((T_Gauss- t0_prop)/T, kep_gauss(:,j));
    xlabel('time [T]');
    ylabel(sprintf('%s', keplabel{j}))
    grid on;
    subplot(2,1,2);
    plot(window_size, SAD(:,j));
    xlabel('window size');
    ylabel('SAD');
    grid on;
end

% choosing the window size when the SAD value is more than 99% of the final value for each keplerian element
for j = 1:6
    index(j) = find(SAD(:,j) > 0.99 * SAD(end,j), 1, 'first');  
    window_size_secular(j) = window_size(index(j));
end

window_size_long = window_size_secular./1.5; % window size narrowing for long-term effects filtering

% movmean low-pass filtering function for each keplerian element
kep_filtered_long_a = movmean(kep_gauss(:,1),window_size_long(1),1); 
kep_filtered_long_e = movmean(kep_gauss(:,2),window_size_long(2),1); 
kep_filtered_long_i = movmean(kep_gauss(:,3),window_size_long(3),1);
kep_filtered_long_OM = movmean(kep_gauss(:,4),window_size_long(4),1);
kep_filtered_long_om = movmean(kep_gauss(:,5),window_size_long(5),1);
kep_filtered_long_th = movmean(kep_gauss(:,6),window_size_long(6),1);

kep_filtered_secular_a = movmean(kep_gauss(:,1),window_size_secular(1),1);
kep_filtered_secular_e = movmean(kep_gauss(:,2),window_size_secular(2),1);
kep_filtered_secular_i = movmean(kep_gauss(:,3),window_size_secular(3),1);
kep_filtered_secular_OM = movmean(kep_gauss(:,4),window_size_secular(4),1);
kep_filtered_secular_om = movmean(kep_gauss(:,5),window_size_secular(5),1);
kep_filtered_secular_th = movmean(kep_gauss(:,6),window_size_secular(6),1);

% plotting unfiltered and filtered keplerian elements in both long and secular cases
figure()
subplot(3,2,1)
plot((T_Gauss-t0_prop)/T,kep_gauss(:,1),'LineWidth',2);
hold on
plot((T_Gauss-t0_prop)/T,kep_filtered_secular_a,'g','LineWidth',2);
plot((T_Gauss-t0_prop)/T,(kep_filtered_long_a),'r','LineWidth',2);
xlabel('time [T]');
ylabel(' a [km]');
legend('Complete','Secular','long');
title('Filtered semi-major axis');
grid on

subplot(3,2,2)
plot((T_Gauss-t0_prop)/T,kep_gauss(:,2),'LineWidth',2);
hold on
plot((T_Gauss-t0_prop)/T,kep_filtered_secular_e,'g','LineWidth',2);
plot((T_Gauss-t0_prop)/T,(kep_filtered_long_e),'r','LineWidth',2);
xlabel('time [T]');
ylabel(' e [-]');
legend('Complete','Secular','long');
title('Filtered eccentricity');
grid on

subplot(3,2,3)
plot((T_Gauss-t0_prop)/T,kep_gauss(:,3),'LineWidth',2);
hold on
plot((T_Gauss-t0_prop)/T,kep_filtered_secular_i,'g','LineWidth',2);
plot((T_Gauss-t0_prop)/T,(kep_filtered_long_i),'r','LineWidth',2);
xlabel('time [T]');
ylabel(' i [deg]');
legend('Complete','Secular','long');
title('Filtered inclination');
grid on

subplot(3,2,4)
plot((T_Gauss-t0_prop)/T,kep_gauss(:,4),'LineWidth',2);
hold on
plot((T_Gauss-t0_prop)/T,kep_filtered_secular_OM,'g','LineWidth',2);
plot((T_Gauss-t0_prop)/T,(kep_filtered_long_OM),'r','LineWidth',2);
xlabel('time [T]');
ylabel(' \Omega [deg]');
legend('Complete','Secular','long');
title('Filtered RAAN');
grid on

subplot(3,2,5)
plot((T_Gauss-t0_prop)/T,kep_gauss(:,5),'LineWidth',2);
hold on
plot((T_Gauss-t0_prop)/T,kep_filtered_secular_om,'g','LineWidth',2);
plot((T_Gauss-t0_prop)/T,(kep_filtered_long_om),'r','LineWidth',2);
xlabel('time [T]');
ylabel(' \omega [deg]');
legend('Complete','Secular','long');
title('Filtered argument of periapsis');
grid on

subplot(3,2,6)
plot((T_Gauss-t0_prop)/T,kep_gauss(:,6),'LineWidth',2);
hold on
plot((T_Gauss-t0_prop)/T,kep_filtered_secular_th,'g','LineWidth',2);
plot((T_Gauss-t0_prop)/T,(kep_filtered_long_th),'r','LineWidth',2);
xlabel('time [T]');
ylabel(' \theta [deg]');
legend('Complete','Secular','long');
title('Filtered true anomaly');
grid on

%% COMPARISON WITH REAL DATA
%Satellite SL-12 R/B - NORAD Catalog number: 27504

kep0_real = [2.6447e+04, 0.5946, deg2rad(23.0470), deg2rad(94.9916), deg2rad(29.2472), deg2rad(87.4165)]; %initial keplerian elements of the satellite
T_real= 2*pi*sqrt(kep0_real(1)^3/muP);    % orbital period of the satellite [s]                   
    
initial_date = [2020,12,01,00,00,00];   % inizial S/C date   
final_date = [2021,12,01,00,00,00];     % final S/C date          
t_in = date2mjd2000(initial_date)*86400;    % [s]         
t_fin = date2mjd2000(final_date)*86400;     % [s]                                  
step = 86400;              
tspan_real = [t_in:step:t_fin];  

Gauss_eqs = @(t,kep) ode_gaussEsq_RSW_J2_MOON (t,kep,muP,kep_pert_RSW (kep,muP,RE,J2),kepMoon_pert_RSW (t,kep,muMoon,muP)); 
[T_real_gauss,kep_real_gauss] = ode113(Gauss_eqs,tspan_real,kep0_real,options); % Gauss propagation

comparisonData = readmatrix('comparisonData.xls'); % reading the Excell file and converting in a matlab matrix

for k=1:1:length(comparisonData(:,1))
    mjd2000_date(k) = jd2mjd2000(comparisonData(k,1)); % converting date
end

T_real_eph = (mjd2000_date')*86400; %[sec]

% building the kep vector from the columns of the excel file
kep_real_eph(:,1) = comparisonData(:,11); % a [km]
kep_real_eph(:,2) = comparisonData(:,2); % e [-]
kep_real_eph(:,3) = deg2rad(comparisonData(:,4)); % i [rad]
kep_real_eph(:,4) = deg2rad(comparisonData(:,5)); % OM [rad]
kep_real_eph(:,5) = deg2rad(comparisonData(:,6)); % om [rad]
kep_real_eph(:,6) = deg2rad(comparisonData(:,10)); % th [rad]

kep_real_gauss(:,3:6) = rad2deg(unwrap(kep_real_gauss(:,3:6))); % unwrapping and converting in deg kep elements
kep_real_eph(:,3:6) = rad2deg(unwrap(kep_real_eph(:,3:6)));  

% plotting comparison graphs between keplerian elements obtained by gauss
% propagation and ephemerides from space track and nasa horizon
keplabel = {'a [km]', 'e', 'i [deg]', '\Omega [deg]', '\omega [deg]', '\theta [deg]'}; % string vector with labels
figure('Name','Gauss propagation and NASA Horizon ephemerides comparison');
for k = 1:6
    subplot(3,2,k)
    plot((T_real_gauss-t_in)/T_real, kep_real_gauss(:,k),'LineWidth',1.5);
    hold on
    plot((T_real_eph-t_in)/T_real, kep_real_eph(:,k),'LineWidth', 1.5);
    xlabel('periods [T]')
    ylabel(sprintf('%s', keplabel{k}))
    legend('Gauss propagation','NASA Horizon Ephemerides');
    grid on
end
