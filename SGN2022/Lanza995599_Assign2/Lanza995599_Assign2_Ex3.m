% Spacecraft Guidance and Navigation (2022/2023)
% Assignment # 2
% Exercise # 3
% Author: Davide Lanza

%% EXERCISE 3.1: Measurements Acquisition 

% Clear workspace, close all figures, and clear command window
clearvars; close all; clc;

% Clear kernels
cspice_kclear(); 

% loading SPICE kernels:
cspice_furnsh('assignment02.tm') ; 
addpath('sgp4');
addpath('simulator_ex3');

% Reference epoch 
t0_ref = '2023-04-01T14:55:12.023'; % [UTC]
et0 = cspice_str2et(t0_ref) ; %  [s]

% Initial states (quaternion and angular velocity) @target in BODY -> INERTIAL 
q0 = [ 0.674156352338764  ; % Initial quaternion [-]
       0.223585877389611  ; 
       0.465489474399161  ; 
       0.528055032413102 ]; 

omega0 = [ -0.001262427155865  ; % Initial angular velocity [-]
            0.001204540074343  ; 
           -0.000039180139156 ]; 


% Gravitational parameter 
mu = cspice_bodvrd( 'EARTH', 'GM', 1 );  % [km^3/s^2]

% GEO orbit
radius = 42164 ;                 % Radius of the orbit        [km]
T = 2*pi * sqrt(radius^3/mu) ;   % Orbital period             [s]
n = 2*pi/T ;                     % Mean motion of the target  [rad/s]

% Initial state nominal trajectory 
rr0 = [12.0; -60.0; 0.0] ;              % Position 
vv0 = [1e-4; -2*n*rr0(1); -1.2e-3] ;    % Velocity 
xx0_nom_traj = [rr0; vv0] ;                      % State 

etf = et0 + cspice_spd ; % 1 day [s]
tspan = et0 : 1 : etf ; % 1 Hz of frequency 
tspan_hr = (tspan-et0)/3600 ; 
options = odeset('AbsTol', 1e-13, 'RelTol', 1e-13);

% Propagate the nominal trajctory in LVLH frame
[~, xx_nom_traj] = ode113(@CW_ode, tspan, xx0_nom_traj, options);

% Plotting the nominal trajectory
figure(1)
plot3(xx_nom_traj(:,1), xx_nom_traj(:,2), xx_nom_traj(:,3), 'k', 'LineWidth', 2);
hold on;
grid on;
plot3(0, 0, 0, 'b.', 'MarkerSize', 35);

xlabel('x [m]', 'Interpreter', 'latex')
ylabel('y [m]', 'Interpreter', 'latex')
zlabel('z [m]', 'Interpreter', 'latex')
legend('OSR nominal trajectory', 'SGN-I', 'Location', 'Best', 'Interpreter', ...
        'latex', 'Fontsize', 15);
title('Relative nominal trajectory of OSR wrt SGN (@SGN-I LVLH)', 'Interpreter', ...
        'Latex', 'fontsize', 20)
set(gca, 'FontSize', 15);

% Attitude determination 
om_q0 = [omega0; q0] ; 
[~, om_q] = ode113(@AD_ode, tspan, om_q0, options); 

figure(2)
% Angular velocity plot
subplot(2,1,1)
plot(tspan_hr, om_q(:,1:3), 'LineWidth', 1)
grid on
xlabel('Time [hours]', 'Interpreter', 'latex')
ylabel('$\omega$ $[\frac{rad}{s}]$', 'Interpreter', 'latex')
title('SGN-I Angular Velocity (@SGN-I Inertial)', 'Interpreter', 'latex')
legend('$\omega_{x}$', '$\omega_{y}$', '$\omega_{z}$', 'Interpreter', 'latex')
set(gca, 'FontSize',15)
xlim([0 24])

% Quaternion plot 
subplot(2,1,2)
plot(tspan_hr, om_q(:,4:7), 'LineWidth', 1)
grid on
xlabel('Time [hours]', 'Interpreter', 'latex')
ylabel('$q$ $[-]$', 'Interpreter', 'latex')
title('Quaternion (@SGN-I Inertial)', 'Interpreter', 'latex')
legend('${q}_{1}$', '${q}_{2}$', '${q}_{3}$', '${q}_{4}$', 'Interpreter', 'latex')
set(gca, 'FontSize', 15)
xlim([0 24])

% Parameters of the stereo camera
Camera.f = 30;                           % Focal Length [mm]                                    
Camera.d = 54;                           % Pixel density [pix/mm]    
Camera.p0 = [960 600];                   % Central pixel coordinates [pix] 
Camera.b = 1;                            % Baseline [m]                       
Camera.Cframe = [1 0 0; 0 0 -1; 0 1 0] ; % Direction cosine matric
Camera.R = 10;                           % Variance                                  
Camera.SensorSize = [1920 1200];         % Sensor size [pix]                    

% Allocation
visible_meas_vect = zeros(length(tspan), 8) ;
meas_vis_coord = zeros(3, 8, length(tspan)) ; 
n_visible_vert = zeros(length(tspan), 1) ; 

% Iterating for each time istant of the time window 
for i = 1 : length(tspan)
    % Measurement simulation   
    measurement = meas_sim(n, xx_nom_traj(i,1:3)', om_q(i,4:7)', tspan(i) - et0, et0, Camera);

    % Each i-th entry identifies the vertex that generated the corresponding measurement 
    visible_meas_vect(i,1:length(measurement.visible)) = measurement.visible ; 
   
    % Matrix containing the triplet of horizontal pixel, vertical pixel, and baseline
    meas_vis_coord(:, 1 : length(measurement.visible), i) = measurement.y ; 

    % Number of visible vertices 
    n_visible_vert(i) = length(measurement.visible) ; 
end

% PLot of visible vertices 
figure(3)
plot(tspan_hr, n_visible_vert, '.k', 'LineWidth', 1.5)
grid on
xlabel('Time [hours]', 'Interpreter', 'latex')
ylabel('Number of vertices [-]', 'Interpreter', 'Latex')
title('Number of visibile vertices at each measurement instant', 'Interpreter', 'Latex')
set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 15) ; 
xlim([0 24])


%% EXERCISE 3.2: FoV verification 
figure(4);
% Reshape the measurement coordinate matrix to a 2xN matrix
meas_vis_coord_reshaped = reshape(meas_vis_coord(1:2, :, :), 2, []) ;

% Plot the scattered points
scatter(meas_vis_coord_reshaped(1, :), meas_vis_coord_reshaped(2, :), ...
        2, 'MarkerEdgeColor', '#8B0000', 'MarkerFaceColor', '#8B0000');

hold on; grid on;
% Sensor FoV
rectangle('Position', [0, 0, Camera.SensorSize(1), Camera.SensorSize(2)], ...
           'EdgeColor', 'k', 'LineWidth', 0.3) ;
xlim([-100 2000])
ylim([-100 1300])

% Set the legend, labels, and title
legend('coordinates vertex', 'Interpreter', 'latex', 'FontSize', 12) ;
text(1600, 30, 'Sensor FOV', 'Color', 'k', 'LineWidth', 1, 'FontSize', 12, 'Interpreter', 'latex') ;
xlabel('X [pix]', 'Interpreter', 'latex', 'FontSize', 14) ;
ylabel('Y [pix]', 'Interpreter', 'latex', 'FontSize', 14) ;
title(['Vertex Pixel Coordinates of the Measurements (@chaser and fixed with' ...
       'camera)'], 'Interpreter', 'latex', 'FontSize', 14) ;

%% EXERCISE 3.3: EKF & UKF

% Close all figures, and clear command window
close all ; clc ;

% Estimate of the initial relative state state at t0 @LVLH 
rr0_mean = [15.792658268071492 ; -59.044939772661586 ;  3.227106250277039] ; % Mean position [m] 
vv0_mean = [-0.053960274403210 ;  -0.053969644762889 ; -0.089140748762173] ; % Mean velocity [m/s]
xx0_mean = [rr0_mean ; vv0_mean] ; % Mean state

P0 = diag([10, 10, 10, 0.1, 0.1, 0.1]) ; % Covariance matrix [m^2, m^2/s, m^2/s^2]

q = om_q(:, 4:7) ; % quaternion 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EKF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic ; 
[xx_EKF, P_EKF] = EKF(xx0_mean, P0, q, meas_vis_coord, visible_meas_vect, tspan) ; 
CPU_time_EKF = toc ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UKF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic ; 
[xx_UKF, P_UKF] = UKF(xx0_mean, P0, q, meas_vis_coord, visible_meas_vect, tspan) ; 
CPU_time_UKF = toc ; 

%% EXERCISE 3.4: EKF & UKF results and comparison 

close all ; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EKF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocation 
err_EKF_x = zeros(length(tspan), 1) ; 
err_EKF_y = zeros(length(tspan), 1) ; 
err_EKF_z = zeros(length(tspan), 1) ;
err_EKF_vx = zeros(length(tspan), 1) ; 
err_EKF_vy = zeros(length(tspan), 1) ; 
err_EKF_vz = zeros(length(tspan), 1) ; 

sigma_EKF_x = zeros(length(tspan), 1) ; 
sigma_EKF_y = zeros(length(tspan), 1) ; 
sigma_EKF_z = zeros(length(tspan), 1) ; 
sigma_EKF_vx = zeros(length(tspan), 1) ; 
sigma_EKF_vy = zeros(length(tspan), 1) ; 
sigma_EKF_vz = zeros(length(tspan), 1) ; 

% Errors and Standard deviations computation 
for i=1:length(tspan)
    % State errors
    err_EKF_x(i)=norm(xx_EKF(1,i)-xx_nom_traj(i,1)') ; % Position 
    err_EKF_y(i)=norm(xx_EKF(2,i)-xx_nom_traj(i,2)') ;
    err_EKF_z(i)=norm(xx_EKF(3,i)-xx_nom_traj(i,3)') ;

    err_EKF_vx(i)=norm(xx_EKF(4,i)-xx_nom_traj(i,4)') ; % Velocity
    err_EKF_vy(i)=norm(xx_EKF(5,i)-xx_nom_traj(i,5)') ;
    err_EKF_vz(i)=norm(xx_EKF(6,i)-xx_nom_traj(i,6)') ;

    % Sqare root of covariance matrix diagonal elements(3*sigma)
    sigma_EKF_x(i)=3*sqrt(P_EKF(1,1,i)) ; % Position 
    sigma_EKF_y(i)=3*sqrt(P_EKF(2,2,i)) ;
    sigma_EKF_z(i)=3*sqrt(P_EKF(3,3,i)) ;

    sigma_EKF_vx(i)=3*sqrt(P_EKF(4,4,i)) ; % Velocity 
    sigma_EKF_vy(i)=3*sqrt(P_EKF(5,5,i)) ;
    sigma_EKF_vz(i)=3*sqrt(P_EKF(6,6,i)) ;
end

% EKF: Plot results 
figure(1)
subplot(3,1,1)
semilogy(tspan_hr,err_EKF_x,'LineWidth',1.5, 'Color','#8B0000');
hold on ; grid on ;  
semilogy(tspan_hr, sigma_EKF_x,'LineWidth',1.5, 'Color','#000000');
ylabel('$\epsilon_x$ [m]', 'Interpreter', 'Latex', 'fontsize', 18)
legend('error', '$3\sigma$', 'Interpreter', 'Latex', 'FontSize', 15, 'Location', 'best');

subplot(3,1,2)
semilogy(tspan_hr, err_EKF_y,'LineWidth',1.5, 'Color','#8B0000');
hold on ; grid on ; 
semilogy(tspan_hr, sigma_EKF_y,'LineWidth',1.5, 'Color','#000000');
ylabel('$\epsilon_y$ [m]', 'Interpreter', 'Latex', 'fontsize', 18)

subplot(3,1,3)
semilogy(tspan_hr,err_EKF_z,'LineWidth',1.5, 'Color','#8B0000');
hold on ; grid on ; 
semilogy(tspan_hr,sigma_EKF_z,'LineWidth',1.5, 'Color','#000000');
xlabel('Time [hours]', 'Interpreter', 'Latex', 'fontsize', 18)
ylabel('$\epsilon_z$ [m]', 'Interpreter', 'Latex', 'fontsize', 18)
sgtitle('EKF: Position Error Analysis \& Standard Deviation (@SGN-I LVLH)', ...
        'FontSize', 15, 'Interpreter', 'Latex');

figure(2)
subplot(3,1,1)
semilogy(tspan_hr,err_EKF_vx,'LineWidth',1.5, 'Color','#8B0000');
hold on ; grid on ; 
semilogy(tspan_hr,sigma_EKF_vx,'LineWidth',1.5, 'Color','#000000');
ylabel('$\epsilon_{v_x}$ [m/s]', 'Interpreter', 'Latex', 'fontsize', 18)
legend('error', '$3\sigma$', 'Interpreter', 'Latex', 'FontSize', 15, 'Location', 'best');

subplot(3,1,2)
semilogy(tspan_hr,err_EKF_vy,'LineWidth',1.5, 'Color','#8B0000');
hold on ; grid on ; 
semilogy(tspan_hr,sigma_EKF_vy,'LineWidth',1.5, 'Color','#000000');
ylabel('$\epsilon_{v_y}$ [m/s]', 'Interpreter', 'Latex', 'fontsize', 18)

subplot(3,1,3)
semilogy(tspan_hr,err_EKF_vz,'LineWidth',1.5, 'Color','#8B0000');
hold on ; grid on ; 
semilogy(tspan_hr,sigma_EKF_vz,'LineWidth',1.5, 'Color','#000000');
ylabel('$\epsilon_{v_z}$ [m/s]', 'Interpreter', 'Latex', 'fontsize', 18)
xlabel('Time [hours]', 'Interpreter', 'Latex', 'fontsize', 18)
sgtitle('EKF: Velocity Error Analysis \& Standard Deviation (@SGN-I LVLH)', ...
        'FontSize', 15, 'Interpreter', 'Latex');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UKF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocation 
err_UKF_x = zeros(length(tspan), 1) ; 
err_UKF_y = zeros(length(tspan), 1) ; 
err_UKF_z = zeros(length(tspan), 1) ;
err_UKF_vx = zeros(length(tspan), 1) ; 
err_UKF_vy = zeros(length(tspan), 1) ; 
err_UKF_vz = zeros(length(tspan), 1) ; 

sigma_UKF_x = zeros(length(tspan), 1) ; 
sigma_UKF_y = zeros(length(tspan), 1) ; 
sigma_UKF_z = zeros(length(tspan), 1) ; 
sigma_UKF_vx = zeros(length(tspan), 1) ; 
sigma_UKF_vy = zeros(length(tspan), 1) ; 
sigma_UKF_vz = zeros(length(tspan), 1) ; 

% Errors and Standard deviations computation 
for i=1:length(tspan)
    % State errors
    err_UKF_x(i)=norm(xx_UKF(1,i)-xx_nom_traj(i,1)') ; % Position 
    err_UKF_y(i)=norm(xx_UKF(2,i)-xx_nom_traj(i,2)') ;
    err_UKF_z(i)=norm(xx_UKF(3,i)-xx_nom_traj(i,3)') ;

    err_UKF_vx(i)=norm(xx_UKF(4,i)-xx_nom_traj(i,4)') ; % Velocity
    err_UKF_vy(i)=norm(xx_UKF(5,i)-xx_nom_traj(i,5)') ;
    err_UKF_vz(i)=norm(xx_UKF(6,i)-xx_nom_traj(i,6)') ;

    % Sqare root of covariance matrix diagonal elements(3*sigma)
    sigma_UKF_x(i)=3*sqrt(P_UKF(1,1,i)) ; % Position 
    sigma_UKF_y(i)=3*sqrt(P_UKF(2,2,i)) ;
    sigma_UKF_z(i)=3*sqrt(P_UKF(3,3,i)) ;

    sigma_UKF_vx(i)=3*sqrt(P_UKF(4,4,i)) ; % Velocity 
    sigma_UKF_vy(i)=3*sqrt(P_UKF(5,5,i)) ; 
    sigma_UKF_vz(i)=3*sqrt(P_UKF(6,6,i)) ;
end

% UKF: Plot results 
figure(3)
subplot(3,1,1)
semilogy(tspan_hr,err_UKF_x,'LineWidth',1.5, 'Color','#8B0000');
hold on ; grid on ;  
semilogy(tspan_hr,sigma_UKF_x,'LineWidth',1.5, 'Color','#000000');
ylabel('$\epsilon_x$ [m]', 'Interpreter', 'Latex', 'fontsize', 18)
legend('error', '$3\sigma$', 'Interpreter', 'Latex', 'FontSize', 15, 'Location', 'best');

subplot(3,1,2)
semilogy(tspan_hr,err_UKF_y,'LineWidth',1.5, 'Color','#8B0000');
hold on ; grid on ; 
semilogy(tspan_hr,sigma_UKF_y,'LineWidth',1.5, 'Color','#000000');
ylabel('$\epsilon_y$ [m]', 'Interpreter', 'Latex', 'fontsize', 18)

subplot(3,1,3)
semilogy(tspan_hr,err_UKF_z,'LineWidth',1.5, 'Color','#8B0000');
hold on ; grid on ; 
semilogy(tspan_hr,sigma_UKF_z,'LineWidth',1.5, 'Color','#000000');
xlabel('Time [hours]', 'Interpreter', 'Latex', 'fontsize', 18)
ylabel('$\epsilon_z$ [m]', 'Interpreter', 'Latex', 'fontsize', 18)
sgtitle('UKF: Position Error Analysis \& Standard Deviation (@SGN-I LVLH)', ...
        'FontSize', 15, 'Interpreter', 'Latex');

figure(4)
subplot(3,1,1)
semilogy(tspan_hr,err_UKF_vx,'LineWidth',1.5, 'Color','#8B0000');
hold on ; grid on ; 
semilogy(tspan_hr,sigma_UKF_vx,'LineWidth',1.5, 'Color','#000000');
ylabel('$\epsilon_{v_x}$ [m/s]', 'Interpreter', 'Latex', 'fontsize', 18)
legend('error', '$3\sigma$', 'Interpreter', 'Latex', 'FontSize', 15, 'Location', 'best');

subplot(3,1,2)
semilogy(tspan_hr,err_UKF_vy,'LineWidth',1.5, 'Color','#8B0000');
hold on ; grid on ; 
semilogy(tspan_hr,sigma_UKF_vy,'LineWidth',1.5, 'Color','#000000');
ylabel('$\epsilon_{v_y}$ [m/s]', 'Interpreter', 'Latex', 'fontsize', 18)

subplot(3,1,3)
semilogy(tspan_hr,err_UKF_vz,'LineWidth',1.5, 'Color','#8B0000');
hold on ; grid on ; 
semilogy(tspan_hr,sigma_UKF_vz,'LineWidth',1.5, 'Color','#000000');
ylabel('$\epsilon_{v_z}$ [m/s]', 'Interpreter', 'Latex', 'fontsize', 18)
xlabel('Time [hours]', 'Interpreter', 'Latex', 'fontsize', 18)
sgtitle('UKF: Velocity Error Analysis \& Standard Deviation (@SGN-I LVLH)', ...
        'FontSize', 15, 'Interpreter', 'Latex');

%%%%%%%%%%%%%%%%%%%%%%%%% Comparison of results %%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
subplot(3,1,1)
semilogy(tspan_hr, abs(err_UKF_x-err_EKF_x),'LineWidth',1.5, 'Color','#8B0000');
hold on ; grid on ;  
semilogy(tspan_hr,abs(sigma_UKF_x-sigma_EKF_x),'LineWidth',1.5, 'Color','#000000');

ylabel('$\epsilon_x$ [m]', 'Interpreter', 'Latex', 'fontsize', 18)
legend('$|\epsilon_{UKF}-\epsilon_{EKF}|$ [m]', '$|3\sigma_{UKF}-3\sigma_{EKF}|$ [m]', ...
        'Interpreter', 'Latex', 'FontSize', 15, 'Location', 'best');

subplot(3,1,2)
semilogy(tspan_hr,abs(err_UKF_y-err_EKF_y),'LineWidth',1.5, 'Color','#8B0000');
hold on ; grid on ;  
semilogy(tspan_hr,abs(sigma_UKF_y-sigma_EKF_y),'LineWidth',1.5, 'Color','#000000');
ylabel('$\epsilon_y$ [m]', 'Interpreter', 'Latex', 'fontsize', 18)

subplot(3,1,3)
semilogy(tspan_hr,abs(err_UKF_z-err_EKF_z),'LineWidth',1.5, 'Color','#8B0000');
hold on ; grid on ;  
semilogy(tspan_hr,abs(sigma_UKF_z-sigma_EKF_z),'LineWidth',1.5, 'Color','#000000');
ylabel('$\epsilon_z$ [m]', 'Interpreter', 'Latex', 'fontsize', 18)
xlabel('Time [hours]', 'Interpreter', 'Latex', 'fontsize', 18)

sgtitle('UKF \& EKF: Difference in Position  (@SGN-I LVLH)', 'FontSize', 15, 'Interpreter', 'Latex');

figure(6)
subplot(3,1,1)
semilogy(tspan_hr,abs(err_UKF_vx-err_EKF_vx),'LineWidth',1.5, 'Color','#8B0000');
hold on ; grid on ;  
semilogy(tspan_hr,abs(sigma_UKF_vx-sigma_EKF_vx),'LineWidth',1.5, 'Color','#000000');

ylabel('$\epsilon_x$ [m/s]', 'Interpreter', 'Latex', 'fontsize', 18)
legend('$|\epsilon_{UKF}-\epsilon_{EKF}|$ [m/s]', '$|3\sigma_{UKF}-3\sigma_{EKF}|$ [m/s]', ...
        'Interpreter', 'Latex', 'FontSize', 15, 'Location', 'best');

subplot(3,1,2)
semilogy(tspan_hr,abs(err_UKF_vy-err_EKF_vy),'LineWidth',1.5, 'Color','#8B0000');
hold on ; grid on ;  
semilogy(tspan_hr,abs(sigma_UKF_vy-sigma_EKF_vy),'LineWidth',1.5, 'Color','#000000');
ylabel('$\epsilon_y$ [m/s]', 'Interpreter', 'Latex', 'fontsize', 18)

subplot(3,1,3)
semilogy(tspan_hr,abs(err_UKF_vz-err_EKF_vz),'LineWidth',1.5, 'Color','#8B0000');
hold on ; grid on ;  
semilogy(tspan_hr,abs(sigma_UKF_vz-sigma_EKF_vz),'LineWidth',1.5, 'Color','#000000');
ylabel('$\epsilon_z$ [m/s]', 'Interpreter', 'Latex', 'fontsize', 18)
xlabel('Time [hours]', 'Interpreter', 'Latex', 'fontsize', 18)

sgtitle('UKF \& EKF: Difference in Velocity  (@SGN-I LVLH)', 'FontSize', 15, 'Interpreter', 'Latex');


%% FUNCTIONS

function xx_dot = CW_ode(~, xx)
%
% xx_dot = CW_ode(~, xx)
%
% Clohessy-Wiltshire equations: Provide the odefun for the dynamic
%                               integration in the LVLH reference frame 
%                               centred @SGN-I
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

% Gravitational parameter 
mu = 3.986004354360959e+05 ;  

% Geostationary orbit (GEO)
radius = 42164 ;                 % Radius of the orbit        [km]
T = 2*pi * sqrt(radius^3/mu) ;   % Orbital period             [s]
n = 2*pi/T ;                     % Mean motion of the target  [rad/s]

% Unpacking the state 
rr = xx(1:3) ; % position vector 
vv = xx(4:6) ; % velocity vector 

% State vector derivatives
xx_dot = [                               vv ; 
            3 * n^2 * rr(1) + 2 * n * vv(2) ;
                             -2 * n * vv(1) ;
                              -n^2 * rr(3) ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xxSTM_dot = CW_STM_ode(~,xxSTM)
% 
% xxSTM_dot = CW_STM_ode(~, xxSTM)
% 
% Clohessy-Wiltshire equations: Provide the odefun for dynamic integration
% in the LVLH reference frame centered at SGN-I.
%
% INPUT:
%   xxSTM       State and State Transition Matrix       [42x1]
%                      - Position                        [3x1]  [km]
%                      - Velocity                        [3x1]  [km/s]
%                      - STM                             [36x1] 
%
% OUTPUT:
%   xxSTM_dot   State and STM derivative                 [42x1]
%                - Velocity derivative                    [3x1] [km/s]
%                - Acceleration                           [3x1] [km/s^2]   
%                - STM derivative                        [36x1]
%
% AUTHOR:
%   Davide Lanza
% 

% Gravitational parameter 
mu = 3.986004354360959e+05 ;  

% Geostationary orbit (GEO)
radius = 42164 ;                 % Radius of the orbit        [km]
T = 2*pi * sqrt(radius^3/mu) ;   % Orbital period             [s]
n = 2*pi/T ;                     % Mean motion of the target  [rad/s]

% Unpacking the state
rr = xxSTM(1:3) ; % Position 
vv = xxSTM(4:6) ; % Velocity 

% State vector derivatives
xx_dot = [                               vv ; 
            3 * n^2 * rr(1) + 2 * n * vv(2) ;
                             -2 * n * vv(1) ;
                              -n^2 * rr(3)] ;
% State transition matrix 
A = [0           , 0           , 0           , 1    , 0    , 0;
     0           , 0           , 0           , 0    , 1    , 0;
     0           , 0           , 0           , 0    , 0    , 1;
     3*n^2       , 0           , 0           , 0    , 2*n  , 0;
     0           , 0           , 0           ,-2*n  , 0    , 0;
     0           , 0           ,-n^2         , 0    , 0    , 0];

% Computing STM derivative
STM = (reshape(xxSTM(7:42),6,6)) ; % Reshaping into a matrix
STM_dot = A * STM ; 

STM_dot = reshape(STM_dot, [36 1]); % Reshaped STM derivatives 

% Building the ouput 
xxSTM_dot = [xx_dot ; STM_dot] ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function om_q_dot = AD_ode (~, om_q)
%
% om_q_dot = AD_ode (~, om_q)
% 
% Angular dynamics equations: Provides the odefun for angular dynamics 
% integration in quaternion representation.
%
% INPUT:
%   om_q     : Angular velocity and quaternion            [7x1]
%              - Angular velocity                         [3x1]    [rad/s]
%              - Quaternion                               [4x1]       [-]
%
% OUTPUT:
%   om_q_dot : Angular velocity and quaternion derivative [7x1]
%              - Angular acceleration                     [3x1]   [rad/s^2]
%              - Quaternion derivative                    [4x1]
%
% AUTHOR:
%   Davide Lanza
%
        
% Parameters of SGN-I
l = 10 ; h = 5 ; d = 3 ;                             % Size [m]
rho = 1420 ;                                         % Density [kg/m^3]
m = l * h * d * rho ;                                % Mass [kg]
J = m/12 * diag([d^2 + h^2, l^2 + h^2, l^2 + d^2]) ; % Inertia matrix [kg/m^2]

% Unpack
om = om_q(1:3) ;    % angular velocity [rad/s]
q = om_q(4:7) ;     % quaternion       [-]

% Angular velocity derivatives (from Euler equations)
om_dot = J \ (-cross(om, J*om)) ; 

% quaternio derivative 
qdot = (1/2) * [    0   -om(1)  -om(2)  -om(3);
                 om(1)      0    om(3)  -om(2);
                 om(2)  -om(3)      0    om(1);
                 om(3)   om(2)  -om(1)      0 ] * q ;

% Packing into a unique vector the output 
om_q_dot = [om_dot; qdot] ; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xx, P] = EKF(xx0, P0, q, meas_vis_coord, visible_meas_vect, tspan)  
%
% [xx, P] = EKF(xx0, P0, q, meas_vis_coord, visible_meas_vect, tspan)  
%
% Extended Kalman Filter: Implements the EKF algorithm for state estimation.
%
% INPUT:
%   xx0                 : Initial state estimate [6x1]
%                         - Position [3x1] [km]
%                         - Velocity [3x1] [km/s]
%   P0                  : Initial covariance matrix [6x6]
%   q                   : Quaternion [Nx4]
%   meas_vis_coord      : Measured visible coordinates [3xN]
%   visible_meas_vect   : Visibility matrix [NxN]
%   tspan               : Time vector [1xM] [s]
%
% OUTPUT:
%   xx                  : Estimated state [6xM]
%   P                   : Covariance matrix [6x6xM]
%
% AUTHOR:
%   Davide Lanza
%

options = odeset('RelTol',1e-13,'AbsTol',1e-13) ; 
xx(:,1) = xx0 ; 
P(:,:,1) = P0 ; 

% Reference epoch 
t0_ref = '2023-04-01T14:55:12.023'; % [UTC]
et0=cspice_str2et(t0_ref) ; % [s]


% Iteration from time t_(k-1) to t_k
for k = 2 : length(tspan)

    % State and Covariance at step k-1
    xx_old = xx(:,k-1) ; 
    P_old = P(:,:,k-1) ;

    % Propagation for t_(k-1) to t_k 
    [~, xxSTM] = ode113(@CW_STM_ode, [tspan(k-1) tspan(k)], [xx_old; reshape(eye(6),[], 1)] , options) ;
    
    % Perform the prediction step 
    xx_k_min  = xxSTM(end,1:6)' ; % Mean 

    STM_k = reshape(xxSTM(end,7:42), [6 6]) ;
    P_k_min = STM_k * P_old * STM_k' ; % Covariance 

    %Real measurements
    vis_vert = nonzeros(visible_meas_vect(k,:)) ;
    
    % Measurement computation coming from model 
    [xx_vert, H_k] = EKF_measurement(tspan(k)-et0, tspan(1), xx_k_min(1:3), q(k,:)') ; 
    
    % Measurement from simulation 
    est_meas = nonzeros(xx_vert(:,vis_vert)) ;

    % Compute the Jacobian matrix  
    H = zeros(3*length(vis_vert),6) ;
    for j = 1 : length(vis_vert)
        H((3*j-2): 3*j, :) = H_k(:, :, vis_vert(j)) ;
    end
    
    % Compute noise matrix 
    meas_nois = diag([10, 10, 10]);
    R_k = [] ;
    for i = 1 : length(vis_vert)
        R_k = blkdiag(R_k, meas_nois) ;
    end

    % Compute Kalman gain 
    K_k = P_k_min * H'/ (H * P_k_min * H'+ R_k) ;
    
    % Compute a posteriori mean and covariance 
    meas_sim = nonzeros(meas_vis_coord(:,:,k)) ; % measurement from simulation 
    
    Id = eye(6) ; % Identity matrix 
    xx_k_plus = xx_k_min + K_k*(meas_sim-est_meas) ;
    P_plus = (Id - K_k * H) * P_k_min ;

    % Building the output 
    xx(:,k) = xx_k_plus ;
    P(:,:,k) = P_plus ;

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xx, P] = UKF(xx0, P0, q, meas_vis_coord, visible_meas_vect, tspan)  
%
% [xx, P] = UKF(xx0, P0, q, meas_vis_coord, visible_meas_vect, tspan)  
%
% Unscented Kalman Filter: Implements the UKF algorithm for state estimation.
%
% INPUT:
%   xx0                 : Initial state estimate [6x1]
%                         - Position [3x1] [km]
%                         - Velocity [3x1] [km/s]
%   P0                  : Initial covariance matrix [6x6]
%   q                   : Quaternion [Nx4]
%   meas_vis_coord      : Measured visible coordinates [3xN]
%   visible_meas_vect   : Visibility matrix [NxN]
%   tspan               : Time vector [1xM] [s]
%
% OUTPUT:
%   xx                  : Estimated state [6xM]
%   P                   : Covariance matrix [6x6xM]
%
% AUTHOR:
%   Davide Lanza
%

% Reference epoch 
t0_ref = '2023-04-01T14:55:12.023'; % [UTC]
et0=cspice_str2et(t0_ref) ; % [s]

% UT parameters  
n = 6 ; % size of P
alpha = 1e-3 ;  % Determines the spread of sigma points
k = 0 ; 
lambda = alpha^2 * (n+k) - n ; % Scaling paramenter 
beta = 2 ; 

% Weights computation
W_mean = zeros(2*n +1, 1) ; 
W_cov = zeros(2*n +1, 1) ; 

W_mean(1) = lambda/(n+lambda) ; 
W_cov(1) = lambda/(n+lambda) + (1 - alpha^2 + beta) ; 

for i = 1 : 2*n 
    W_mean(i+1) = 1/(2*(n+lambda)) ;
    W_cov(i+1) = 1/(2*(n+lambda)) ;
 end

% Initial conditions 
xx(:,1) = xx0 ; 
P(:,:,1) = P0 ; 

% Iteration from time t_(k-1) to t_k
for k = 2 : length(tspan)
    % State and Covariance at step k-1
    xx_old = xx(:, k-1) ; 
    P_old = P(:,:, k-1);

    P_sqr_root = sqrtm((n + lambda)*P_old) ; 

    % Initial Sigma points composition 
    Sigma_pnt0(:, 1) = xx_old ; % Sigma_0 
    for i = 1 : n
        Sigma_pnt0(:, i+1) = xx_old + P_sqr_root(:,i) ; 
    end
    for i = n+1 : 2*n 
        Sigma_pnt0(:, i+1) = xx_old - P_sqr_root(:,i-n) ; 
    end

    % Sigma points propagation
    options = odeset('AbsTol', 1e-13, 'RelTol', 1e-13) ;
    Y = zeros(n,(2*n)+1) ; 
    for i = 1 : (2*n)+1
        [~, Y_iesimo] = ode113(@(t, xx) CW_ode(t, xx), [tspan(k-1) tspan(k)], Sigma_pnt0(:,i), options) ; 
        Y(:,i) = Y_iesimo(end, :) ; 
    end

    % Retrieving visible measurements 
    vis_vert = nonzeros(visible_meas_vect(k,:));
    
    est_meas = [] ; % Variable allocatoion 
    % Computing sigma points in the measurement space 
    for i = 1 : (2*n)+1
        [xx_vert] = UKF_measurement(tspan(k) - et0, tspan(1), Y(1:3 , i), q(k,:)') ; 
        est_meas(:,i) = nonzeros(xx_vert(:,vis_vert));
    end
    
    % Compute a priori mean state  
    xx_k_min = zeros(n, 1) ; 
    for i = 1 : (2*n + 1)
        xx_k_min = W_mean(i) *  Y(:,i) + xx_k_min ; 
    end

    % Compute a priori covariance matrix 
    % Weighted sample covariance computation 
    P_k_min = zeros(n,n) ;   
    for i = 1 : (2*n + 1) 
        P_k_min = W_cov(i) * ((Y(:,i) - xx_k_min)*(Y(:,i) - xx_k_min)') + P_k_min ; 
    end

    % Compute the a priori mean of the measurements
    yy_k_min = zeros(size(est_meas,1),1) ; 
    for i = 1 : (2*n + 1)
        yy_k_min = W_mean(i) * est_meas(:,i) + yy_k_min ; 
    end

    % Perform the update step
    % Compute measurements covariance and cross covariance
    P_yy_k = zeros(size(est_meas,1),size(est_meas,1)) ; 
    for i = 1 : (2*n + 1) 
        P_yy_k = W_cov(i) * ((est_meas(:,i) - yy_k_min)*(est_meas(:,i) - yy_k_min)') + P_yy_k ; 
    end

    meas_nois = diag([10, 10, 10]) ;
    R_k = [] ;
    for i = 1 : length(vis_vert)
        R_k = blkdiag(R_k, meas_nois) ;
    end

    P_yy_k = P_yy_k + R_k ; 

    % Compute measurements covariance and cross covariance
    P_xy_k = zeros(6,size(est_meas,1)) ; 
    for i = 1 : (2*n + 1) 
        P_xy_k = W_cov(i) * ((Y(:,i) - xx_k_min)*(est_meas(:,i) - yy_k_min)') + P_xy_k ; 
    end
    
    % Compute Kalman gain
    K_k = P_xy_k / P_yy_k ;

    % Compute a posteriori mean and covariance 
    meas_sim = nonzeros(meas_vis_coord(:,:,k)); % measurement from simulation 

    xx_k_plus = xx_k_min + K_k * (meas_sim - yy_k_min) ;
    P_plus = P_k_min - K_k * P_yy_k*K_k' ;

    % Building the output 
    xx(:,k) = xx_k_plus ;
    P(:,:,k) = P_plus ;

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xx_vert, H]=EKF_measurement(t,~ ,xx_k, q)
% 
% [xx_vert, H]=EKF_measurement(t,~ ,xx_k, q)
%
% Measurement model for Extended Kalman Filter: Computes the measured 
% visible coordinates and the measurement Jacobian matrix.
%
% INPUT:
%   t             : Current time [s]
%   xx_k          : Estimated state [3x1] [km]
%   q             : Quaternion [-]
%
% OUTPUT:
%   xx_vert       : Measured visible coordinates [3x8] [pix]
%   H             : Measurement Jacobian matrix [3x6x8]
%
% AUTHOR:
%   Davide Lanza
%

% Parameters of the stereo camera
foc = 30 ;                         % Focal Length [mm]                                    
D = 54 ;                           % Pixel density [pix/mm]    
p0 = [960 600] ;                   % Central pixel coordinates [pix] 
b = 1 ;                            % Baseline [m]                       
Cframe = [1 0  0;
          0 0 -1; 
          0 1  0] ;  % Direction cosine matric

% Gravitational parameter 
mu = cspice_bodvrd( 'EARTH', 'GM', 1 );  % [km^3/s^2]

% GEO orbit
radius = 42164 ;                 % Radius of the orbit        [km]
T = 2*pi * sqrt(radius^3/mu) ;   % Orbital period             [s]
n = 2*pi/T ;                     % Mean motion of the target  [rad/s]

% Parameters of SGN-I
l = 10 ; h = 5 ; d = 3 ;                             % Size [m]

% Vertices 
p =       [ l/2 -d/2 -h/2 ; 
           l/2  d/2 -h/2 ; 
           l/2  d/2  h/2 ; 
           l/2 -d/2  h/2 ; 
          -l/2 -d/2 -h/2 ; 
          -l/2  d/2 -h/2 ; 
          -l/2  d/2  h/2 ; 
          -l/2 -d/2  h/2 ] ; 

% Director cosine matrix  LVLH to Inertial 
C_L_I=[cos(n*t), sin(n*t), 0 ;
      -sin(n*t), cos(n*t), 0 ;
              0,        0, 1];

% Direction cosine matrix from quaternion 
A_B_I = quat2dcm(q') ;

% Direction cosine matrix from body to LVLH
A_B_L = A_B_I / C_L_I ;

%Vertices coordinates in camera reference frame 
vert_coord=zeros(3,8) ;
p = p' ; 

for i=1:8
    vert_coord(:,i) = -Cframe * xx_k(1:3) + Cframe * (A_B_L \ p(:,i)) ;
end

% Camera measurement model in camera fixed frame 
yy = @(X,Y,Z) [p0(1)-D*foc*Y/Z ;
               p0(2)+D*foc*X/Z ; 
                    b*D*foc/Z] ;

xx_vert=zeros(3,8) ;

for i=1:8
    xx_vert(:,i) = yy(vert_coord(1,i), vert_coord(2,i), vert_coord(3,i)) ;
end

% dh/dX Jacobian matrix 
H=zeros(3,6,8);

H_til=@(X,Y,Z) [      0, -D*foc/Z,  D*foc*Y/Z^2 ;
                D*foc/Z,        0, -D*foc*X/Z^2 ;
                      0,        0, -b*D*foc/Z^2];

% Building the final Jacobian matrix
for i=1:8
    H(:,:,i)=[-H_til(vert_coord(1,i), vert_coord(2,i), ...
               vert_coord(3,i)) * Cframe, zeros(3,3)] ;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xx_vert] = UKF_measurement(t,~ ,xx_k, q)
% 
% [xx_vert, H]=UKF_measurement(t,~ ,xx_k, q)
%
% Measurement model for Unscented Kalman Filter: Computes the measured 
% visible coordinates.
%
% INPUT:
%   t             : Current time [s]
%   xx_k          : Estimated state [3x1] [km]
%   q             : Quaternion [-]
%
% OUTPUT:
%   xx_vert       : Measured visible coordinates [3x8] [pix]
%   
% AUTHOR:
%   Davide Lanza
%

% Parameters of the stereo camera
foc = 30 ;                           % Focal Length [mm]                                    
D = 54 ;                           % Pixel density [pix/mm]    
p0 = [960 600] ;                   % Central pixel coordinates [pix] 
b = 1 ;                            % Baseline [m]                       
Cframe = [1 0  0; 
          0 0 -1; 
          0 1  0] ; % Direction cosine matrix

% Gravitational parameter 
mu = cspice_bodvrd( 'EARTH', 'GM', 1 ) ;  % [km^3/s^2]

% GEO orbit
radius = 42164 ;                 % Radius of the orbit        [km]
T = 2*pi * sqrt(radius^3/mu) ;   % Orbital period             [s]
n = 2*pi/T ;                     % Mean motion of the target  [rad/s]

% Parameters of SGN-I
l = 10 ; h = 5 ; d = 3 ;      % Size [m]

% Vertices position  
p =       [ l/2 -d/2 -h/2 ; 
           l/2  d/2 -h/2 ; 
           l/2  d/2  h/2 ; 
           l/2 -d/2  h/2 ; 
          -l/2 -d/2 -h/2 ; 
          -l/2  d/2 -h/2 ; 
          -l/2  d/2  h/2 ; 
          -l/2 -d/2  h/2 ] ; 

% Director cosine matrix  LVLH to Inertial 
C_L_I=[cos(n*t), sin(n*t), 0 ;
      -sin(n*t), cos(n*t), 0 ;
              0,        0, 1];

% Direction cosine matrix from quaternion 
A_B_I = quat2dcm(q') ;

% Direction cosine matrix from body to LVLH
A_B_L = A_B_I / C_L_I ;

%Vertices coordinates in camera reference frame 
vert_coord = zeros(3,8) ;
p = p' ; 

for i=1:8
    vert_coord(:,i) = -Cframe * xx_k(1:3) + Cframe * (A_B_L \ p(:,i)) ;
end

% Camera measurement model in camera fixed frame 
yy = @(X,Y,Z) [p0(1) - D*foc*Y/Z ;
               p0(2) + D*foc*X/Z ; 
                       b*D*foc/Z] ;

% Allocation 
xx_vert = zeros(3,8) ;

% Evaluating the output 
for i=1:8
    xx_vert(:,i) = yy(vert_coord(1,i), vert_coord(2,i), vert_coord(3,i)) ;
end

end

