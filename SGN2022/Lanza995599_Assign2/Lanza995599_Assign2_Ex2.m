% Spacecraft Guidance and Navigation (2022/2023)
% Assignment # 2
% Exercise # 2
% Author: Davide Lanza

%% EXERCISE 2.1: Visibility Windows

% Clear workspace, close all figures, and clear command window
clearvars; close all; clc;

% Clear kernels
cspice_kclear(); 

% loading SPICE kernels:
cspice_furnsh('assignment02.tm') ; 

% Reference epoch in UTC 
t0 = '2022-11-11T19:08:49.824';
et0 = cspice_str2et(t0) ; % Ephemeris time [s]

tf = '2022-11-12T04:30:00.000' ;
etf = cspice_str2et(tf) ; % Ephemeris time [s]

% Initial mean state from ex 01 (ECI J2000)
rr0 = [6054.30795817484; -3072.03883303992; -133.115352431876]; % Position [km]
vv0 = [4.64750094824087;  9.18608475681236; -0.62056520749034]; % Velocity [km/s]

xx0 = [rr0; vv0] ; % Initial mean state vector 
tspan = [et0, etf] ; 

% Station names 
station_Ko = 'KOUROU';
station_Pe = 'PERTH';

% Propagation 
options = odeset('AbsTol', 1e-13, 'RelTol', 1e-13);
[~, xx] = ode113(@TBP_ode, tspan, xx0, options) ;

% Updating the initial state for exercise 2 (ECI J2000)
xx0 = xx(end,:)' ; % [6x1] 

% Initial time 
t0 = '2022-11-12T04:30:00.000' ; 
et0 = cspice_str2et(tf) ; % Ephemeris time [s]

% Final time 
tf = '2022-11-14T16:30:00.000' ;
etf = cspice_str2et(tf) ; % Ephemeris time [s]

% Time grid definition 
timegrid = et0 : 60 : etf ; 
timegrid(end) = etf ; 

% Orbit propagation over the et0 and etf
[~, xx] = ode113(@TBP_ode, timegrid, xx0, options) ;

% Compute Azimuth & Elevation 
[Az_Ko, El_Ko, ~] = SC_LatCordinates(station_Ko, timegrid, xx);
[Az_Pe, El_Pe, ~] = SC_LatCordinates(station_Pe, timegrid, xx);

%%%%%%%%%%%%%%%%%%%%%%%%% Time Visibility Windws %%%%%%%%%%%%%%%%%%%%%%%%%

% Create the logical vector for visible elevations @KOUROU
i_visible_Ko = El_Ko > 10 * cspice_rpd;

% Extract the time and visible elevations @KOUROU
time_visible_Ko = timegrid(i_visible_Ko);
elevation_visible_Ko = El_Ko(i_visible_Ko);

% Create the logical vector for visible elevations @PERTH
i_visible_Pe = El_Pe > 5 * cspice_rpd;

% Extract the time and visible elevations @PERTH
time_visible_Pe = timegrid(i_visible_Pe);
elevation_visible_Pe = El_Pe(i_visible_Pe);

% Dates and times to be labeled on the x-axis
x_dates = [
cspice_str2et('2022-11-12T04:30:00') ;
cspice_str2et('2022-11-12T12:30:00') ;
cspice_str2et('2022-11-12T20:30:00') ;
cspice_str2et('2022-11-13T04:30:00') ;
cspice_str2et('2022-11-13T12:30:00') ;
cspice_str2et('2022-11-13T20:30:00') ;
cspice_str2et('2022-11-14T04:30:00') ;
cspice_str2et('2022-11-14T12:30:00') ;
cspice_str2et('2022-11-14T20:30:00')];

% Elevation(t) @KOUROU
figure(1)
plot(timegrid, El_Ko * cspice_dpr, 'LineWidth', 1.5);
hold on;
plot(time_visible_Ko, elevation_visible_Ko * cspice_dpr, '.', ...
     'MarkerSize', 8, 'Color', '#FF6600');
grid on;

ylabel('Elevation [deg]','Interpreter', 'latex','FontSize', 20)
xlabel('Date', 'Interpreter', 'latex','FontSize', 20)
title('Satellite Elevation (@KOUROU-Topocentric)', 'Interpreter', 'latex')
legend('Elevation', 'Visible Elevation', 'Location', 'Best', 'Interpreter', 'latex')
xticks(x_dates);
xticklabels(cspice_et2utc(x_dates', 'C', 0));
xtickangle(30);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12)

% Azimuth(t) @KOUROU
figure(2)
plot(timegrid, Az_Ko * cspice_dpr, 'LineWidth', 1.5);
grid on;

ylabel('Azimuth [deg]','Interpreter', 'latex','FontSize', 20)
xlabel('Date', 'Interpreter', 'latex','FontSize', 20)
title('Satellite Azimuth (@KOUROU-Topocentric)', 'Interpreter', 'latex')
legend('Azimuth', 'Location', 'Best', 'Interpreter', 'latex')
xticks(x_dates);
xticklabels(cspice_et2utc(x_dates', 'C', 0));
xtickangle(30);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12)

% Elevation(t) @PERTH
figure(3)
plot(timegrid, El_Pe * cspice_dpr, 'LineWidth', 1.5);
hold on;
plot(time_visible_Pe, elevation_visible_Pe * cspice_dpr, '.', 'MarkerSize', 8, 'Color', '#800000');
grid on;

ylabel('Elevation [deg]','Interpreter', 'latex','FontSize', 20)
xlabel('Date', 'Interpreter', 'latex','FontSize', 20)
title('Satellite Elevation (@PERTH-Topocentric)', 'Interpreter', 'latex')
legend('Elevation', 'Visible Elevation', 'Location', 'Best', 'Interpreter', 'latex')
xticks(x_dates);
xticklabels(cspice_et2utc(x_dates', 'C', 0));
xtickangle(30);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12)

% Azimuth(t) @PERTH
figure(4)
plot(timegrid, Az_Pe * cspice_dpr, 'LineWidth', 1.5);
grid on;

ylabel('Azimuth [deg]','Interpreter', 'latex','FontSize', 20)
xlabel('Date', 'Interpreter', 'latex','FontSize', 20)
title('Satellite Azimuth (@PERTH-Topocentric)', 'Interpreter', 'latex')
legend('Azimuth', 'Location', 'Best', 'Interpreter', 'latex')
xticks(x_dates);
xticklabels(cspice_et2utc(x_dates', 'C', 0));
xtickangle(30);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12)

% Combined visibility windows
figure(5)
plot(time_visible_Ko, elevation_visible_Ko * cspice_dpr, '.', 'MarkerSize', 8, 'Color', '#FF6600');
hold on ; grid on ; 
plot(time_visible_Pe, elevation_visible_Pe * cspice_dpr, '.', 'MarkerSize', 8, 'Color', '#800000');

ylabel('Elevation [deg]','Interpreter', 'latex','FontSize', 20)
xlabel('Date', 'Interpreter', 'latex','FontSize', 20)
title('Combined Visibility Windows (@Topocentric)', 'Interpreter', 'latex')
legend('@KOUROU', '@PERTH', ...
       'Location', 'Best', 'Interpreter', 'latex')
xticks(x_dates);
xticklabels(cspice_et2utc(x_dates', 'C', 0));
xtickangle(30);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12)

% Combined visibile Altitude and Azimuth 
figure(6)
plot(Az_Ko(i_visible_Ko)*cspice_dpr, elevation_visible_Ko * cspice_dpr, '.', 'MarkerSize', 8, 'Color', '#FF6600');
grid on ; hold on ; 
plot(Az_Pe(i_visible_Pe)*cspice_dpr, elevation_visible_Pe * cspice_dpr, '.', 'MarkerSize', 8, 'Color', '#800000');

ylabel('Elevation [deg]','Interpreter', 'latex','FontSize', 20)
xlabel('Azimuth [deg]', 'Interpreter', 'latex','FontSize', 20)
title('Visibility Windows: Azimuth \& Elevation (@Topocentric)', 'Interpreter', 'latex')
legend('@KOUROU', '@PERTH', 'Location', 'Best', 'Interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12)

% Printing Table with visibility windows
fprintf('Station: %s\n', 'KOUROU');
fprintf('   Pass #1:\n');
fprintf('       Start time: %s\n', cspice_et2utc(timegrid(65), 'ISOC', 0));
fprintf('       End time: %s\n', cspice_et2utc(timegrid(556), 'ISOC', 0));
fprintf('   Pass #2:\n');
fprintf('       Start time: %s\n', cspice_et2utc(timegrid(1771), 'ISOC', 0));
fprintf('       End time: %s\n', cspice_et2utc(timegrid(2006), 'ISOC', 0));
fprintf('   Pass #3:\n');
fprintf('       Start time: %s\n', cspice_et2utc(timegrid(2871), 'ISOC', 0));
fprintf('       End time: %s\n', cspice_et2utc(timegrid(3545), 'ISOC', 0));
fprintf('----------------------------------\n')

fprintf('\nStation: %s\n', 'PERTH');
fprintf('   Pass #1:\n');
fprintf('       Start time: %s\n', cspice_et2utc(timegrid(588), 'ISOC', 0));
fprintf('       End time: %s\n', cspice_et2utc(timegrid(1504), 'ISOC', 0));
fprintf('   Pass #2:\n');
fprintf('       Start time: %s\n', cspice_et2utc(timegrid(2229), 'ISOC', 0));
fprintf('       End time: %s\n', cspice_et2utc(timegrid(2837), 'ISOC', 0));
fprintf('----------------------------------\n')


%% EXERCISE 2.2: Simulate Measurements 

% Close all figures, and clear command window
close all; clc;

addpath('sgp4');

typerun = 'u'   ; % User-provided inputs to SGP4 Matlab function
opsmode = 'a'   ; % afspc approach ('air force space command')
whichconst = 72 ; % WGS72 constants (radius, gravitational parameter)

% Two-Line elements 
TLE_1 = '1 87654U 22110B   22316.00967942  .00000002  00000-0  32024-3 0  9990';
TLE_2 = '2 87654   3.6309 137.1541 8138191 196.1632  96.6141  1.26411866   834';

% Initialize the satrec structure, using the function twoline2rv
satrec = twoline2rv(TLE_1, TLE_2, typerun, 'e', opsmode, whichconst) ; 

% TLE reference epoch 
[year, mon, day, hr, min, sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f',[year,mon,day,hr,min,sec]);
et0 = cspice_str2et(sat_epoch_str);

% Conversion factor from arcseconds to radians
arcsec2rad = pi / (180 * 3600);

% Nutation correction coefficients
ddpsi = -0.113603 * arcsec2rad; % Nutation correction in radians
ddeps = -0.007049 * arcsec2rad; % Nutation correction in radians

% Allocation
rr_ECI = zeros(3, length(timegrid)) ; % Position @ECI 
vv_ECI = zeros(3, length(timegrid)) ; % Velocity @ECI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2.2.a %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(timegrid)

    % SGP4 propagator 
    tsince = (timegrid(i) - et0) / 60; % Convert time to minutes from TLE epoch
    [~, rteme, vteme] = sgp4(satrec, tsince); % Calculate position and velocity in TEME frame
    
    % Centuries from 2000-01-01T00:00:00.00 TDT
    ttt = cspice_unitim(timegrid(i), 'ET', 'TDT') / (cspice_jyear() * 100);
    
    % Convert from TEME to ECI frame
    [rr_ECI(:, i), vv_ECI(:, i)] = teme2eci(rteme, vteme, [0; 0; 0], ttt, ddpsi, ddeps);
end

% Compute Azimuth & Elevation 
[Az_Ko_SGP4, El_Ko_SGP4, rho_Ko_SGP4] = SC_LatCordinates(station_Ko, timegrid, [rr_ECI', vv_ECI']) ; 
[Az_Pe_SGP4, El_Pe_SGP4, rho_Pe_SGP4]= SC_LatCordinates(station_Pe, timegrid, [rr_ECI', vv_ECI']) ; 


% Create the logical vector for visible elevations @KOUROU
i_visible_Ko_SGP4 = El_Ko_SGP4 > 10 * cspice_rpd;

% Extract the time and visible elevations @KOUROU
time_visible_Ko_SGP4 = timegrid(i_visible_Ko_SGP4);
elevation_visible_Ko_SGP4 = El_Ko_SGP4(i_visible_Ko_SGP4);

% Create the logical vector for visible elevations @PERTH
i_visible_Pe_SGP4 = El_Pe_SGP4 > 5 * cspice_rpd;

% Extract the time and visible elevations @PERTH
time_visible_Pe_SGP4 = timegrid(i_visible_Pe_SGP4);
elevation_visible_Pe_SGP4 = El_Pe_SGP4(i_visible_Pe_SGP4);

% Plots of the results 
figure(1)
plot(time_visible_Ko_SGP4, elevation_visible_Ko_SGP4 * cspice_dpr, '.', ...
     'MarkerSize', 8, 'Color', '#FF6600');
hold on ; grid on ;
plot(time_visible_Pe_SGP4, elevation_visible_Pe_SGP4 * cspice_dpr, '.', ...
     'MarkerSize', 8, 'Color', '#800000');

ylabel('Elevation [deg]','Interpreter', 'latex','FontSize', 20)
xlabel('Date', 'Interpreter', 'latex','FontSize', 20)
title('Combined Visibility Windows (with SGP4, @Topocentric)', 'Interpreter', 'latex')
legend('@KOUROU', '@PERTH', ...
       'Location', 'Best', 'Interpreter', 'latex')
xticks(x_dates);
xticklabels(cspice_et2utc(x_dates', 'C', 0));
xtickangle(30);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12)

% Combined visibile Altitude and Azimuth 
figure(2)
plot(Az_Ko_SGP4(i_visible_Ko_SGP4)*cspice_dpr, elevation_visible_Ko_SGP4 * cspice_dpr, '.', 'MarkerSize', 8, 'Color', '#FF6600');
grid on ; hold on ; 
plot(Az_Pe_SGP4(i_visible_Pe_SGP4)*cspice_dpr, elevation_visible_Pe_SGP4 * cspice_dpr, '.', 'MarkerSize', 8, 'Color', '#800000');

ylabel('Elevation [deg]','Interpreter', 'latex','FontSize', 20)
xlabel('Azimuth [deg]', 'Interpreter', 'latex','FontSize', 20)
title('Visibility Windows: Azimuth \& Elevation (with SGP4, @Topocentric)', 'Interpreter', 'latex')
legend('@KOUROU', '@PERTH', 'Location', 'Best', 'Interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12)


figure(3)
plot(time_visible_Ko_SGP4, rho_Ko_SGP4(i_visible_Ko_SGP4), '.', ...
     'MarkerSize', 8, 'Color', '#FF6600');
hold on ; grid on ;
plot(time_visible_Pe_SGP4, rho_Pe_SGP4(i_visible_Pe_SGP4), '.', ...
     'MarkerSize', 8, 'Color', '#800000');

ylabel('Range [km]','Interpreter', 'latex','FontSize', 20)
xlabel('Date', 'Interpreter', 'latex','FontSize', 20)
title('RANGE: Combined Visibility Windows (with SGP4, @Topocentric)', 'Interpreter', 'latex')
legend('@KOUROU', '@PERTH', ...
       'Location', 'Best', 'Interpreter', 'latex')
xticks(x_dates);
xticklabels(cspice_et2utc(x_dates', 'C', 0));
xtickangle(30);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2.2.b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Building covariance matrix and mean 
sigma_Az_El = 100e-3 ; % [deg]
sigma_range = 0.01 ;   % [km]  

R = diag([sigma_Az_El^2, sigma_Az_El^2, sigma_range^2]) ; 
mu_Ko = [Az_Ko_SGP4(i_visible_Ko_SGP4)*cspice_dpr, El_Ko_SGP4(i_visible_Ko_SGP4)* cspice_dpr,...
         rho_Ko_SGP4(i_visible_Ko_SGP4)] ; 

mu_Pe = [Az_Pe_SGP4(i_visible_Pe_SGP4)* cspice_dpr, El_Pe_SGP4(i_visible_Pe_SGP4)* cspice_dpr,...
          rho_Pe_SGP4(i_visible_Pe_SGP4)] ; 

% Random error simultation
meas_noise_Ko = mvnrnd(mu_Ko, R);
meas_noise_Pe = mvnrnd(mu_Pe, R);

figure(4)
plot(time_visible_Ko_SGP4, meas_noise_Ko(:,2), '.', ...
     'MarkerSize', 8, 'Color', '#FF6600');
hold on ; grid on ;
plot(time_visible_Pe_SGP4, meas_noise_Pe(:,2), '.', ...
     'MarkerSize', 8, 'Color', '#800000');

ylabel('Elevation [deg]','Interpreter', 'latex','FontSize', 20)
xlabel('Date', 'Interpreter', 'latex','FontSize', 20)
title('Combined Visibility Windows (with SGP4 \& noise, @Topocentric)', 'Interpreter', 'latex')
legend('@KOUROU', '@PERTH', ...
       'Location', 'Best', 'Interpreter', 'latex')
xticks(x_dates);
xticklabels(cspice_et2utc(x_dates', 'C', 0));
xtickangle(30);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12)

figure(5)
plot(meas_noise_Ko(:,1), meas_noise_Ko(:,2), '.', 'MarkerSize', 8, 'Color', '#FF6600');
grid on ; hold on ; 
plot(meas_noise_Pe(:,1), meas_noise_Pe(:,2), '.', 'MarkerSize', 8, 'Color', '#800000');

ylabel('Elevation [deg]','Interpreter', 'latex','FontSize', 20)
xlabel('Azimuth [deg]', 'Interpreter', 'latex','FontSize', 20)
title('Visibility Windows: Azimuth \& Elevation (with SGP4 \& noise, @Topocentric)', 'Interpreter', 'latex')
legend('@KOUROU', '@PERTH', 'Location', 'Best', 'Interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12)

figure(6)
plot(time_visible_Ko_SGP4, meas_noise_Ko(:,3), '.', ...
     'MarkerSize', 8, 'Color', '#FF6600');
hold on ; grid on ;
plot(time_visible_Pe_SGP4, meas_noise_Pe(:,3), '.', ...
     'MarkerSize', 8, 'Color', '#800000');

ylabel('Range [km]','Interpreter', 'latex','FontSize', 20)
xlabel('Date', 'Interpreter', 'latex','FontSize', 20)
title('RANGE: Combined Visibility Windows (with SGP4 \& noise, @Topocentric)', 'Interpreter', 'latex')
legend('@KOUROU', '@PERTH', ...
       'Location', 'Best', 'Interpreter', 'latex')
xticks(x_dates);
xticklabels(cspice_et2utc(x_dates', 'C', 0));
xtickangle(30);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12)

%% EXERCISE 2.3: Navigation Problem

% Close all figures, and clear command window
close all; clc;

% Optimization options
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'off', ...
                       'StepTolerance', 1e-11, 'OptimalityTolerance',1e-11) ; 

xx0_ECI = [rr_ECI(:,1) ; vv_ECI(:,1)] ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CASE (a) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute least squares
[xx_ECI_a, RESNORM_a, RESIDUAL_a, FLAG_a, OUT_a, LAMBDA_a, JACOBIAN_a] = lsqnonlin(@(xx) ...
     costfunction_a(xx, timegrid, meas_noise_Pe, station_Pe, i_visible_Pe_SGP4), ...
     xx0_ECI, [], [], options) ; 

fprintf('CASE (a): pure Keplerian dynamics || @PERTH\n');
fprintf('       Initial State [km, km/s] = \n');
fprintf('           %12.6f\n', xx_ECI_a);
fprintf('-------------------------------------------\n')

% Covariance matrix 
Jacobian_a = full(JACOBIAN_a);
P_a = RESNORM_a/(length(RESIDUAL_a)-length(xx0_ECI)).* inv(Jacobian_a'*Jacobian_a);
Trace_rr_a = sqrt(trace(P_a(1:3,1:3))) ;
Trace_vv_a = sqrt(trace(P_a(4:6,4:6))) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CASE (b) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xx_ECI_b, RESNORM_b, RESIDUAL_b, FLAG_b, OUT_b, LAMBDA_b, JACOBIAN_b] = lsqnonlin(@(xx) ...
     costfunction_b(xx, timegrid, meas_noise_Ko, meas_noise_Pe, i_visible_Ko_SGP4, i_visible_Pe_SGP4), ...
     xx0_ECI, [], [], options) ; 

fprintf('CASE (b): pure Keplerian dynamics || @KOUROU & @PERTH\n');
fprintf('       Initial State [km, km/s] = \n');
fprintf('           %12.6f\n', xx_ECI_b);
fprintf('-------------------------------------------\n')

% Covariance matrix 
Jacobian_b = full(JACOBIAN_b);
P_b = RESNORM_b/(length(RESIDUAL_b)-length(xx0_ECI)).* inv(Jacobian_b'*Jacobian_b);
Trace_rr_b = sqrt(trace(P_b(1:3,1:3))) ;
Trace_vv_b = sqrt(trace(P_b(4:6,4:6))) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CASE (c) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xx_ECI_c, RESNORM_c, RESIDUAL_c, FLAG_c, OUT_c, LAMBDA_c, JACOBIAN_c] = lsqnonlin(@(xx) ...
     costfunction_c(xx, timegrid, meas_noise_Ko, meas_noise_Pe, i_visible_Ko_SGP4, i_visible_Pe_SGP4), ...
     xx0_ECI, [], [], options) ; 

fprintf('CASE (c): J2-perturbed motion || @KOUROU & @PERTH\n');
fprintf('       Initial State [km, km/s] = \n');
fprintf('           %12.6f\n', xx_ECI_c);
fprintf('-------------------------------------------\n')

% Covariance matrix 
Jacobian_c = full(JACOBIAN_c);
P_c = RESNORM_c/(length(RESIDUAL_c)-length(xx0_ECI)).* inv(Jacobian_c'*Jacobian_c);
Trace_rr_c = sqrt(trace(P_c(1:3,1:3))) ;
Trace_vv_c = sqrt(trace(P_c(4:6,4:6))) ;

%% FUNCTIONS

function xx_dot = TBP_ode(~, xx)
%
% xx_dot = TBP_ode(~, xx)
%
% Two-Body Problem: Provide the odefun for the dynamic integration
%
% INPUT:
%   xx                  State               [6x1]           
%                         - Position        [3x1]           [km]
%                         - Velocity        [3x1]           [km/s]
% 
% OUTPUT:
%   xx_dot              State derivative    [6x1]           [km/s]
%                                                           [km/s^2]
% AUTHOR:
%   Davide Lanza
%

mu = 3.986004354360959e+05 ; % Gravitational parameter 

% Unpacking 
rr = xx(1:3) ; % Position 
vv = xx(4:6) ; % Velocity 

r = norm(rr) ; % Magnitude 

% State vector derivatives
xx_dot = [vv ; 
          (-mu/(r^3)) * rr] ; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xx_dot = TBP_J2_ode(tt, xx)
%
% xx_dot = TBP_J2_ode(~, xx)
%
% Two-Body Problem with J2 Perturbation: Provide the odefun for the dynamic integration
%
% INPUT:
%   xx                  State          [6x1]
%                           - Position [3x1]    [km]
%                           - Velocity [3x1]    [km/s]
%
% OUTPUT:
%   xx_dot             State derivative [6x1]    [km/s], [km/s^2]
% 
% AUTHOR:
%  Davide Lanza
%

% Gravitational parameter of Earth
mu = 3.986004354360959e+05; 

% Earth data
Re = cspice_bodvrd('EARTH', 'RADII', 3);
R_e = Re(1);

% Unpack state variables
x = xx(1);                      % Position 
y = xx(2);
z = xx(3);
norm_r = sqrt(x^2 + y^2 + z^2); % Norm of position 

x_dot = xx(4);                  % Velocity 
y_dot = xx(5);
z_dot = xx(6);

% Compute the matrix to convert position from ECI to ECEF
rotm = cspice_pxform('J2000', 'ITRF93', tt);

% Convert the ECI position in ECEF
rr_ECEF = rotm * [x; y; z];
norm_rr_ECEF = norm(rr_ECEF);

% Compute the J2 acceleration in ECEF:
J2 = 0.0010826269 ; 
a_J2 = 3/2 * mu * J2 * rr_ECEF / norm_rr_ECEF^3 * (R_e / norm_rr_ECEF)^2 .* ...
        (5 * (rr_ECEF(3) / norm_rr_ECEF)^2 - [1; 1; 3]);

% Convert to ECI
a_J2 = rotm \ a_J2;

% Add A_J2 to the Keplerian acceleration

xx_dot = [x_dot; 
          y_dot; 
          z_dot; 
          -mu/norm_r^3 * [x; y; z] + a_J2 ] ;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Az, El, rho] = SC_LatCordinates(stat_name, timegrid, xx_SC_ECI)
%
% [Az, El, rho] = SC_LatCordinates(stat_name, timegrid, xx_SC_ECI)
% 
% SC_LatCoordinates.m - Compute the sensor measurements 
%
% INPUT:
% stat_name              Name of the ground station
% timegrid               Vector of times                    [s]
% xx_SC_ECI              Matrix of states                   [km], [km/s]
%
% Outputs:
% Az                    Azimuth angles                      [rad]
% El                    Elevation angles                    [rad]
% rho                   Ranges                              [km]
%
% AUTHOR:
% Davide Lanza
%


rho = zeros(length(timegrid), 1) ;      % Range             [km]
Az = zeros(length(timegrid), 1) ;       % Azimuth angle     [rad]
El = zeros(length(timegrid), 1) ;       % Elevation angle   [rad]

xx_SC_ECI = xx_SC_ECI' ;

for i = 1 : length(timegrid)

    % State rotation matrix from ECI J2000 to TOPOCENTRIC reference frames
    ROT_ECI2TOPO = cspice_sxform('J2000', [stat_name, '_TOPO'], timegrid(i)) ;
    
    % Station position in ECI
    xx_station_ECI = cspice_spkezr(stat_name, timegrid(i), 'J2000', 'NONE', 'EARTH') ;
    
    % SC state wrt station in ECI
    xx_station_SC_ECI = xx_SC_ECI(:, i) - xx_station_ECI ;
    
    % SC state wrt station in TOPO
    xx_station_SC_TOPO = ROT_ECI2TOPO * xx_station_SC_ECI ;
    
    % SC state from rectangular to latitudinal coordinates
    xx_station_SC_lat = cspice_xfmsta(xx_station_SC_TOPO, 'RECTANGULAR', 'LATITUDINAL', 'EARTH') ;
    
    % Extract azimuth, elevation, range, and range rate
    rho(i) = xx_station_SC_lat(1) ;      % Range             [km]
    Az(i) = xx_station_SC_lat(2) ;       % Azimuth angle     [rad]
    El(i) = xx_station_SC_lat(3) ;       % Elevation angle   [rad]
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function residual = costfunction_a(xx0, timegrid, meas_real, station_name, i_visible)
%
% residual = costfunction_a(xx0, timegrid, meas_real, station_name, i_visible)
%
% costfunction_a: Compute the cost function for case (a) 
%
% INPUT:
% xx0            Initial state vector                   [6x1]
% timegrid       Vector of time points for measurements [Nx1]
% meas_real      Matrix of real measurements            [Nx3]
% station_name   Name of the ground station
% i_visible      Indices of visible satellites
%
% OUTPUT:
% residual      Cost function residual
%
% AUTHOR:
% Davide Lanza
% 

% Initialize output variable
residual = [] ;

% Compute the weights 
sigma_Az_El = 100e-3 ; % [deg]
sigma_range = 0.01 ;   % [km]  

W_m = diag(1./[sigma_Az_El^2, sigma_Az_El^2, sigma_range^2]) ; 

% Propagate x to the epochs of the measurements 
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13) ;
[~, xx_prop] = ode113(@(t,xx) TBP_ode(t,xx), timegrid, xx0, options) ;

% Compute predicted measurements
[Az, El, rho] = SC_LatCordinates(station_name, timegrid, xx_prop) ; 

% Visible
meas_pred(:,1)=Az(i_visible)*cspice_dpr ;
meas_pred(:,2)=El(i_visible)*cspice_dpr ;
meas_pred(:,3)=rho(i_visible) ;

% Compute the residual of the measurements and append it to the output
diff_meas_weighted = W_m*(meas_pred-meas_real)' ;

% Ooutput
residual = [ residual ;
             diff_meas_weighted ] ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function residual = costfunction_b(xx0, timegrid, meas_real_Ko, meas_real_Pe,i_visible_Ko, i_visible_Pe)  
%
% residual = costfunction_b(xx0, timegrid, meas_real_Ko, meas_real_Pe,i_visible_Ko, i_visible_Pe)  
%
% costfunction_b - Compute the cost function for case (b) 
%
% INPUT:
% xx0               Initial state vector 
% timegrid          Vector of time points for measurements
% meas_real_Ko      Matrix of real measurements for KOUROU ground station 
% meas_real_Pe      Matrix of real measurements for PERTH ground station 
% i_visible_Ko      Indices of visible satellites for KOUROU ground station
% i_visible_Pe      Indices of visible satellites for PERTH ground station
%
% OUTPUT:
% residual - Cost function residual
%
% AUTHOR:
% Davide Lanza
%

% Initialize output variable
residual = [] ;

% Compute the weights 
sigma_Az_El = 100e-3 ; % [deg]
sigma_range = 0.01 ;   % [km]  

W_m = diag(1./[sigma_Az_El^2, sigma_Az_El^2, sigma_range^2]) ; 

% Station names 
station_Ko = 'KOUROU';
station_Pe = 'PERTH';

% Propagate x to the epochs of the measurements 
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13) ;
[~, xx_prop] = ode113(@(t,xx) TBP_ode(t,xx), timegrid, xx0, options) ;

% Compute predicted measurements
[Az_Ko, El_Ko, rho_Ko] = SC_LatCordinates(station_Ko, timegrid, xx_prop) ; 
[Az_Pe, El_Pe, rho_Pe] = SC_LatCordinates(station_Pe, timegrid, xx_prop) ; 

% Visible
meas_pred_Ko(:,1) = Az_Ko(i_visible_Ko)*cspice_dpr ;
meas_pred_Ko(:,2) = El_Ko(i_visible_Ko)*cspice_dpr ;
meas_pred_Ko(:,3) = rho_Ko(i_visible_Ko);

meas_pred_Pe(:,1) = Az_Pe(i_visible_Pe) * cspice_dpr ;
meas_pred_Pe(:,2) = El_Pe(i_visible_Pe) * cspice_dpr ;
meas_pred_Pe(:,3) = rho_Pe(i_visible_Pe) ;


%Compute the residual
diff_meas_weighted_Ko = W_m*(meas_pred_Ko - meas_real_Ko)' ;
diff_meas_weighted_Pe = W_m*(meas_pred_Pe - meas_real_Pe)' ;

% Ooutput
residual = [residual ; 
            diff_meas_weighted_Ko(:) ; 
            diff_meas_weighted_Pe(:) ];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function residual = costfunction_c(xx0, timegrid, meas_real_Ko, meas_real_Pe,i_visible_Ko, i_visible_Pe)  
%
% residual = costfunction_c(xx0, timegrid, meas_real_Ko, meas_real_Pe,i_visible_Ko, i_visible_Pe)  
%
% costfunction_b - Compute the cost function for case (c) 
%
% INPUT:
% xx0               Initial state vector 
% timegrid          Vector of time points for measurements
% meas_real_Ko      Matrix of real measurements for KOUROU ground station 
% meas_real_Pe      Matrix of real measurements for PERTH ground station 
% i_visible_Ko      Indices of visible satellites for KOUROU ground station
% i_visible_Pe      Indices of visible satellites for PERTH ground station
%
% OUTPUT:
% residual - Cost function residual
%
% AUTHOR:
% Davide Lanza
%

% Initialize output variable
residual = [] ;

% Compute the weights 
sigma_Az_El = 100e-3 ; % [deg]
sigma_range = 0.01 ;   % [km]  

W_m = diag(1./[sigma_Az_El^2, sigma_Az_El^2, sigma_range^2]) ; 

% Station names 
station_Ko = 'KOUROU';
station_Pe = 'PERTH';

% Propagate x to the epochs of the measurements 
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13) ;
[~, xx_prop] = ode113(@(t,xx) TBP_J2_ode(t,xx), timegrid, xx0, options) ;

% Compute predicted measurements
[Az_Ko, El_Ko, rho_Ko] = SC_LatCordinates(station_Ko, timegrid, xx_prop) ; 
[Az_Pe, El_Pe, rho_Pe] = SC_LatCordinates(station_Pe, timegrid, xx_prop) ; 

% Visible
meas_pred_Ko(:,1) = Az_Ko(i_visible_Ko)*cspice_dpr ;
meas_pred_Ko(:,2) = El_Ko(i_visible_Ko)*cspice_dpr ;
meas_pred_Ko(:,3) = rho_Ko(i_visible_Ko) ;

meas_pred_Pe(:,1) = Az_Pe(i_visible_Pe)*cspice_dpr ;
meas_pred_Pe(:,2) = El_Pe(i_visible_Pe)*cspice_dpr ;
meas_pred_Pe(:,3) = rho_Pe(i_visible_Pe);

%Compute the residual
diff_meas_weighted_Ko = W_m*(meas_pred_Ko - meas_real_Ko)' ;
diff_meas_weighted_Pe = W_m*(meas_pred_Pe - meas_real_Pe)' ;

residual = [ residual ; 
             diff_meas_weighted_Ko(:); 
             diff_meas_weighted_Pe(:) ] ;

end

