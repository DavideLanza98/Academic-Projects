% Spacecraft Guidance and Navigation (2022/2023)
% Assignment # 2
% Exercise # 1
% Author: Davide Lanza

%% EXERCISE 1.1: LinCov & UT

% Clear workspace, close all figures, and clear command window
clearvars; close all; clc;

% Clear kernels
cspice_kclear(); 

% loading SPICE kernels:
cspice_furnsh('assignment02.tm') ; 

% Reference initial epoch in UTC 
t0 = '2022-11-11T19:08:49.824';
et0 = cspice_str2et(t0) ; % Ephemeris time [s]

% Initial mean state (ECI J2000)
rr0 = [6054.30795817484; -3072.03883303992; -133.115352431876]; % Position [km]
vv0 = [4.64750094824087;  9.18608475681236; -0.62056520749034]; % Velocity [km/s]

xx0 = [rr0; vv0] ; % Initial mean state vector 

% Initial covariance (ECI J2000) [km^2, km^2/s, km^2/s^2]
P0 = [ 5.6e-3  3.5e-3 -7.1e-4       0       0      0  ;
       3.5e-3  9.7e-3  7.6e-4       0       0      0  ;    
      -7.1e-4  7.6e-4  8.1e-4       0       0      0  ;     
            0       0       0  2.8e-7       0      0  ;    
            0       0       0       0  2.7e-7      0  ;      
            0       0       0       0       0 9.6e-8 ];

mu = cspice_bodvrd('EARTH', 'GM', 1); % Gravitational parameter [km^3/s^2]
a = norm(rr0)/(2 - norm(rr0)*dot(vv0,vv0)/mu); % Semi-major axis [km]
T = 2*pi*sqrt(a^3/mu); % Period of the orbit [s]

tt = linspace(et0, et0+4*T, 9) ;
PP0 = P0 ; 
xx00 = xx0 ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LinCov %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocation 
xx_mean_LinCov = zeros(6,8) ; 
PP_LinCov = zeros(6,6, 8) ; 

% Mean state and Covariance matrix computation 
for i = 1 : 8 
   [x_mean_LinCov, P_LinCov] = LinCov(xx00, PP0, tt(i), tt(i+1)) ; 
    xx_mean_LinCov(:,i) = x_mean_LinCov ; 
    PP_LinCov(:,:,i) = P_LinCov ; 
    xx00 = x_mean_LinCov ; 
    PP0 = P_LinCov ; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters 
n = 6 ; % size of P
alpha = 1e-3 ;  % Determines the spread of sigma points
k = 0 ; 
lambda = alpha^2 * (n+k) - n ; % Scaling paramenter 

% Matrix square root 
P_sqr_root = sqrtm((n + lambda)*P0) ; 

% Initial Sigma points composition 
Sigma_pnt0(:, 1) = xx0 ; % Sigma_0 
for i = 1 : n
    Sigma_pnt0(:, i+1) = xx0 + P_sqr_root(:,i) ; 
end
for i = n+1 : 2*n 
    Sigma_pnt0(:, i+1) = xx0 - P_sqr_root(:,i-n) ; 
end

% Allocation 
xx_mean_UT = zeros(6,8) ; 
PP_UT = zeros(6,6,8) ; 

% Mean state and Covariance matrix computation 
for i = 1 : 8 
    [x_mean_UT, P_UT, Sigma_pnt] = UT(Sigma_pnt0, tt(i),  tt(i+1)) ; 
    xx_mean_UT(:,i) = x_mean_UT ; 
    PP_UT(:,:,i) = P_UT ; 
    Sigma_pnt0 = Sigma_pnt ; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Total Dispersion %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocation 
LinCov_trace_pos = zeros(8,1) ; 
UT_trace_pos = zeros(8,1) ; 
LinCov_trace_vel = zeros(8, 1) ; 
UT_trace_vel = zeros(8,1) ; 

% Computing the total variance 
for i = 1 : 8 
    LinCov_trace_pos(i) = sqrt(trace(PP_LinCov(1:3,1:3,i))) ; 
    UT_trace_pos(i) = sqrt(trace(PP_UT(1:3,1:3,i))) ; 
    LinCov_trace_vel(i) = sqrt(trace(PP_LinCov(4:6,4:6,i))) ; 
    UT_trace_vel(i) = sqrt(trace(PP_UT(4:6,4:6,i))) ;
end


%% EXERCISE 1.2: Monte-Carlo Simulation 

%%%%%%%%%%%%%%%%%%%%%%%%% Monte-Carlo Simulation %%%%%%%%%%%%%%%%%%%%%%%%%

% Allocation
xx_mean_MC = zeros(6,8) ; 
PP_MC = zeros(6,6, 8) ; 
n = 200 ; % Nuber of samples
R0 = mvnrnd(xx0, P0, n); % Samples generation 
for i = 1 : 8 
    [x_mean_MC, P_MC, R] = MC_sim(R0, tt(i), tt(i+1)) ;  
    xx_mean_MC(:,i) = x_mean_MC ; 
    PP_MC(:,:,i) = P_MC ; 
    R0 = R ; 
end

% Computing the total deviation 
MC_trace_pos = zeros(8,1) ; 
MC_trace_vel = zeros(8,1) ; 
for i = 1 : 8 
    MC_trace_pos(i) = sqrt(trace(PP_MC(1:3,1:3,i))) ; 
    MC_trace_vel(i) = sqrt(trace(PP_MC(4:6,4:6,i))) ; 
end

%%%%%%%%%%%%%%%%%%%%%%%% Plots of total deviations %%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
% P_rr @apogee
subplot(2,1,1)
hold on ; grid on ; 
plot(1:4 , LinCov_trace_pos(1:2:8),'o- b', 'LineWidth', 1.5)
plot(1:4 , UT_trace_pos(1:2:8),'s-- r', 'LineWidth', 1.5)
plot(1:4 , MC_trace_pos(1:2:8),'^-- k', 'LineWidth', 1.5)

title('Total Position Deviation @Apogee', 'Interpreter', 'latex')
xticks(1:4) 
xticklabels({'$\frac{T}{2}$', '$\frac{3}{2}T$', '$\frac{5}{2}T$', '$\frac{7}{2}T$'})
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 15)

% P_rr @perigee
subplot(2,1,2)
hold on; grid on ; 
plot(1:4 , LinCov_trace_pos(2:2:8),'o-- b', 'LineWidth', 1.5)
plot(1:4 , UT_trace_pos(2:2:8),'s-- r', 'LineWidth', 1.5)
plot(1:4 , MC_trace_pos(2:2:8),'^-- k', 'LineWidth', 1.5)

title('Total Position Deviation @Perigee', 'Interpreter', 'latex')
ylabel('Total Deviation [km]','Interpreter', 'latex','FontSize', 20, 'Position', [0.7161 5.3695e+03 -1])
xlabel('N. of Periods [-]', 'Interpreter', 'latex','FontSize', 20)
xticks(1:4)
xticklabels({'$T$', '$2T$', '$3T$', '$4T$'})
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 15)
legend('LinCov', 'UT', 'MC', 'Location', 'Best', 'Interpreter', 'Latex', 'fontsize', 15)

figure(2)
% P_vv @apogee
subplot(2,1,1)
hold on ; grid on ; 
plot(1:4 , LinCov_trace_vel(1:2:8),'o- b', 'LineWidth', 1.5)
plot(1:4 , UT_trace_vel(1:2:8),'s-- r', 'LineWidth', 1.5)
plot(1:4 , MC_trace_vel(1:2:8),'^-- k', 'LineWidth', 1.5)

title('Total Velocity Deviation @Apogee', 'Interpreter', 'latex')
xticks(1:4) 
xticklabels({'$\frac{T}{2}$', '$\frac{3}{2}T$', '$\frac{5}{2}T$', '$\frac{7}{2}T$'})
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 15)

% P_vv @perigee
subplot(2,1,2)
hold on; grid on ; 
plot(1:4 , LinCov_trace_vel(2:2:8),'o-- b', 'LineWidth', 1.5)
plot(1:4 , UT_trace_vel(2:2:8),'s-- r', 'LineWidth', 1.5)
plot(1:4 , MC_trace_vel(2:2:8),'^-- k', 'LineWidth', 1.5)

title('Total Velocity Deviation @Perigee', 'Interpreter', 'latex')
ylabel('Total Deviation [km/s]','Interpreter', 'latex','FontSize', 20, 'Position', [0.7161 4.5 -1])
xlabel('N. of Periods [-]', 'Interpreter', 'latex','FontSize', 20)
xticks(1:4)
xticklabels({'$T$', '$2T$', '$3T$', '$4T$'})
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 15)
legend('LinCov', 'UT', 'MC', 'Location', 'Best', 'Interpreter', 'Latex', 'fontsize', 15)

%%%%%%%%%%%%%%%%%%%%%%%% Plots of cloud of points %%%%%%%%%%%%%%%%%%%%%%%%

% APOGEE - Cloud of Points Sampled with Propagated Covariance
figure(3)
subplot(2,2,1)
[tt, xx] = ode113(@TBP_ode, linspace(et0, et0 + T, 1000), xx0) ; 
plot3(xx(493:507, 1), xx(493:507,2), xx(493:507,3), 'Color',[0.8 0.8 0.8], 'LineWidth', 0.8)
hold on ; grid on ; 
samples_LinCov=mvnrnd(xx_mean_LinCov(1:3,1),PP_LinCov(1:3,1:3, 1), 200);
samples_UT=mvnrnd(xx_mean_UT(1:3,1),PP_UT(1:3,1:3, 1), 200);
samples_MC=mvnrnd(xx_mean_MC(1:3,1),PP_MC(1:3,1:3, 1), 200);

scatter3(samples_LinCov(:,1), samples_LinCov(:,2), samples_LinCov(:,3), 2, 'o b');
scatter3(samples_UT(:,1), samples_UT(:,2), samples_UT(:,3), 2, 's r');
scatter3(samples_MC(:,1), samples_MC(:,2), samples_MC(:,3), 2, '^ k');
view(2)
title('@$\frac{T}{2}$', 'Interpreter', 'latex')
xlabel('x [km]', 'Interpreter', 'Latex', 'fontsize', 12)
ylabel('y [km]', 'Interpreter', 'Latex', 'fontsize', 12)
set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 10);

subplot(2,2,2)
plot3(xx(493:507, 1), xx(493:507,2), xx(493:507,3), 'Color',[0.8 0.8 0.8], 'LineWidth', 0.8)
hold on ; grid on ; 
samples_LinCov=mvnrnd(xx_mean_LinCov(1:3,3),(PP_LinCov(1:3,1:3, 3) + PP_LinCov(1:3,1:3, 3)')/2, 200);
samples_UT=mvnrnd(xx_mean_UT(1:3,3),PP_UT(1:3,1:3, 3), 200);
samples_MC=mvnrnd(xx_mean_MC(1:3,3),PP_MC(1:3,1:3, 3), 200);

scatter3(samples_LinCov(:,1), samples_LinCov(:,2), samples_LinCov(:,3), 3, 'o b');
scatter3(samples_UT(:,1), samples_UT(:,2), samples_UT(:,3), 3, 's r');
scatter3(samples_MC(:,1), samples_MC(:,2), samples_MC(:,3), 3, '^ k');
view(2)
title('@$\frac{3}{2}T$', 'Interpreter', 'latex')
xlabel('x [km]', 'Interpreter', 'Latex', 'fontsize', 12)
ylabel('y [km]', 'Interpreter', 'Latex', 'fontsize', 12)
set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 10);

subplot(2,2,3)
plot3(xx(493:507, 1), xx(493:507,2), xx(493:507,3), 'Color',[0.8 0.8 0.8], 'LineWidth', 0.8)
hold on ; grid on ; 
samples_LinCov=mvnrnd(xx_mean_LinCov(1:3,5),(PP_LinCov(1:3,1:3, 5) + PP_LinCov(1:3,1:3, 5)')/2, 200);
samples_UT=mvnrnd(xx_mean_UT(1:3,5),PP_UT(1:3,1:3, 5), 200);
samples_MC=mvnrnd(xx_mean_MC(1:3,5),PP_MC(1:3,1:3, 5), 200);

scatter3(samples_LinCov(:,1), samples_LinCov(:,2), samples_LinCov(:,3), 3, 'o b');
scatter3(samples_UT(:,1), samples_UT(:,2), samples_UT(:,3), 3, 's r');
scatter3(samples_MC(:,1), samples_MC(:,2), samples_MC(:,3), 3, '^ k');
view(2)
title('@$\frac{5}{2}T$', 'Interpreter', 'latex')
xlabel('x [km]', 'Interpreter', 'Latex', 'fontsize', 12)
ylabel('y [km]', 'Interpreter', 'Latex', 'fontsize', 12)
set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 10);

subplot(2,2,4)
plot3(xx(493:507, 1), xx(493:507,2), xx(493:507,3), 'Color',[0.8 0.8 0.8], 'LineWidth', 0.8)
hold on ; grid on ; 
samples_LinCov=mvnrnd(xx_mean_LinCov(1:3,7),(PP_LinCov(1:3,1:3, 7) + PP_LinCov(1:3,1:3, 7)')/2, 200);
samples_UT=mvnrnd(xx_mean_UT(1:3,7),PP_UT(1:3,1:3, 7), 200);
samples_MC=mvnrnd(xx_mean_MC(1:3,7),PP_MC(1:3,1:3, 7), 200);

scatter3(samples_LinCov(:,1), samples_LinCov(:,2), samples_LinCov(:,3), 3, 'o b');
scatter3(samples_UT(:,1), samples_UT(:,2), samples_UT(:,3), 3, 's r');
scatter3(samples_MC(:,1), samples_MC(:,2), samples_MC(:,3), 3, '^ k');
view(2)
title('@$\frac{7}{2}T$', 'Interpreter', 'latex')
xlabel('x [km]', 'Interpreter', 'Latex', 'fontsize', 12)
ylabel('y [km]', 'Interpreter', 'Latex', 'fontsize', 12)
set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 10);

sgtitle('Cloud of Points at Apogee (@ECI J2000)')
legend('Orbit', 'LinCov', 'UT','MC', 'Location', 'Best', 'Interpreter', 'Latex', 'fontsize', 15)


% PERIGEE - Cloud of Points Sampled with Propagated Covariance
figure(4)
subplot(2,2,1)
plot3(xx(1:8, 1), xx(1:8,2), xx(1:8,3), 'Color',[0.8 0.8 0.8], 'LineWidth', 0.8)
hold on ; grid on ; 
plot3(xx(end-8:end, 1), xx(end-8:end,2), xx(end-8:end,3), 'Color',[0.8 0.8 0.8], 'LineWidth', 0.8)
samples_LinCov=mvnrnd(xx_mean_LinCov(1:3,2),PP_LinCov(1:3,1:3, 2), 200);
samples_UT=mvnrnd(xx_mean_UT(1:3,2),PP_UT(1:3,1:3, 2), 200);
samples_MC=mvnrnd(xx_mean_MC(1:3,2),PP_MC(1:3,1:3, 2), 200);

scatter3(samples_LinCov(:,1), samples_LinCov(:,2), samples_LinCov(:,3), 2, 'o b');
scatter3(samples_UT(:,1), samples_UT(:,2), samples_UT(:,3), 2, 's r');
scatter3(samples_MC(:,1), samples_MC(:,2), samples_MC(:,3), 2, '^ k');
view(2)
title('@$T$', 'Interpreter', 'latex')
xlabel('x [km]', 'Interpreter', 'Latex', 'fontsize', 12)
ylabel('y [km]', 'Interpreter', 'Latex', 'fontsize', 12)
set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 10);

subplot(2,2,2)
plot3(xx(1:8, 1), xx(1:8,2), xx(1:8,3), 'Color',[0.8 0.8 0.8], 'LineWidth', 0.8)
hold on ; grid on ; 
plot3(xx(end-8:end, 1), xx(end-8:end,2), xx(end-8:end,3), 'Color',[0.8 0.8 0.8], 'LineWidth', 0.8)
samples_LinCov=mvnrnd(xx_mean_LinCov(1:3,4),(PP_LinCov(1:3,1:3, 4) + PP_LinCov(1:3,1:3, 4)')/2, 200);
samples_UT=mvnrnd(xx_mean_UT(1:3,4),PP_UT(1:3,1:3, 4), 200);
samples_MC=mvnrnd(xx_mean_MC(1:3,4),PP_MC(1:3,1:3, 4), 200);

scatter3(samples_LinCov(:,1), samples_LinCov(:,2), samples_LinCov(:,3), 3, 'o b');
scatter3(samples_UT(:,1), samples_UT(:,2), samples_UT(:,3), 3, 's r');
scatter3(samples_MC(:,1), samples_MC(:,2), samples_MC(:,3), 3, '^ k');
view(2)
title('@$2T$', 'Interpreter', 'latex')
xlabel('x [km]', 'Interpreter', 'Latex', 'fontsize', 12)
ylabel('y [km]', 'Interpreter', 'Latex', 'fontsize', 12)
set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 10);

subplot(2,2,3)
plot3(xx(1:8, 1), xx(1:8,2), xx(1:8,3), 'Color',[0.8 0.8 0.8], 'LineWidth', 0.8)
hold on ; grid on ; 
plot3(xx(end-8:end, 1), xx(end-8:end,2), xx(end-8:end,3), 'Color',[0.8 0.8 0.8], 'LineWidth', 0.8)
samples_LinCov=mvnrnd(xx_mean_LinCov(1:3,6),(PP_LinCov(1:3,1:3, 6) + PP_LinCov(1:3,1:3, 6)')/2, 200);
samples_UT=mvnrnd(xx_mean_UT(1:3,6),PP_UT(1:3,1:3, 6), 200);
samples_MC=mvnrnd(xx_mean_MC(1:3,6),PP_MC(1:3,1:3, 6), 200);

scatter3(samples_LinCov(:,1), samples_LinCov(:,2), samples_LinCov(:,3), 3, 'o b');
scatter3(samples_UT(:,1), samples_UT(:,2), samples_UT(:,3), 3, 's r');
scatter3(samples_MC(:,1), samples_MC(:,2), samples_MC(:,3), 3, '^ k');
view(2)
title('@$3T$', 'Interpreter', 'latex')
xlabel('x [km]', 'Interpreter', 'Latex', 'fontsize', 12)
ylabel('y [km]', 'Interpreter', 'Latex', 'fontsize', 12)
set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 10);

subplot(2,2,4)
plot3(xx(1:8, 1), xx(1:8,2), xx(1:8,3), 'Color',[0.8 0.8 0.8], 'LineWidth', 0.8)
hold on ; grid on ; 
plot3(xx(end-8:end, 1), xx(end-8:end,2), xx(end-8:end,3), 'Color',[0.8 0.8 0.8], 'LineWidth', 0.8)
samples_LinCov=mvnrnd(xx_mean_LinCov(1:3,8),(PP_LinCov(1:3,1:3, 8) + PP_LinCov(1:3,1:3, 8)')/2, 200);
samples_UT=mvnrnd(xx_mean_UT(1:3,8),PP_UT(1:3,1:3, 8), 200);
samples_MC=mvnrnd(xx_mean_MC(1:3,8),PP_MC(1:3,1:3, 8), 200);

scatter3(samples_LinCov(:,1), samples_LinCov(:,2), samples_LinCov(:,3), 3, 'o b');
scatter3(samples_UT(:,1), samples_UT(:,2), samples_UT(:,3), 3, 's r');
scatter3(samples_MC(:,1), samples_MC(:,2), samples_MC(:,3), 3, '^ k');
view(2)
title('@$4T$', 'Interpreter', 'latex')
xlabel('x [km]', 'Interpreter', 'Latex', 'fontsize', 12)
ylabel('y [km]', 'Interpreter', 'Latex', 'fontsize', 12)
set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 10);

sgtitle('Cloud of Points at Perigee (@ECI J2000)')
legend('Orbit','', 'LinCov', 'UT','MC', 'Location', 'Best', 'Interpreter', 'Latex', 'fontsize', 15)


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xxSTM_dot = TBP_STM_ode(~, xxSTM)
%
% xxSTM_dot = TBP_STM_ode(~, xxSTM)
%
% Two-Body Problem: Provides the ODE function for dynamic integration of
%                   the state transition matrix (STM) along with the state 
%                   vector for a Keplerian two-body problem.
%
% INPUT:
%   xxSTM               State + STM                 [42x1]
%                           - State                 [6x1]
%                           - Position              [3x1]           [km]
%                           - Velocity              [3x1]           [km/s]
%                           - STM                   [36x1] 
%
% OUTPUT:
%   xxSTM_dot          State and STM derivatives    [42x1] 
%
% AUTHOR:
% Davide Lanza
%

mu = 3.986004354360959e+05 ; % Gravitational parameter 

% Unpacking 
rr = xxSTM(1:3) ; % Position 
vv = xxSTM(4:6) ; % Velocity 

r = norm(rr) ; % Magnitude 

% State vector derivatives
xx_dot = [vv ; 
          -mu/(r^3) * rr] ; 

% Position extraction  
r1 = rr(1) ; 
r2 = rr(2) ;
r3 = rr(3) ;

% State transition matrix 
 A =    [                               0,                                0,                                0, 1, 0, 0 ; 
                                        0,                                0,                                0, 0, 1, 0 ;
                                        0,                                0,                                0, 0, 0, 1 ;
         (3*mu*r1^2)/(r)^(5) - mu/(r)^(3),             (3*mu*r1*r2)/(r)^(5),             (3*mu*r1*r3)/(r)^(5), 0, 0, 0 ;
                     (3*mu*r1*r2)/(r)^(5), (3*mu*r2^2)/(r)^(5) - mu/(r)^(3),             (3*mu*r2*r3)/(r)^(5), 0, 0, 0 ;
                     (3*mu*r1*r3)/(r)^(5),             (3*mu*r2*r3)/(r)^(5), (3*mu*r3^2)/(r)^(5) - mu/(r)^(3), 0, 0, 0 ] ;

% Computing STM derivative
STM = reshape(xxSTM(7:42), [6 6]) ; % Reshaping into a matrix
STM_dot = A * STM ; 

STM_dot = reshape(STM_dot, [36 1]); % Reshaped STM derivatives 

% Building the ouput 
xxSTM_dot = [xx_dot ; STM_dot] ; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xx_mean, P] = LinCov(xx0_mean, P0, et0, etf) 
%
% [xx_mean, P] = LinCov(xx0_mean, P0, et0, etf) 
%
% Uncertainty propagation with using a linear method
%
% INPUT:
%   xx0_mean               Initial mean state       [6x1]               
%                           - Position              [3x1]           [km]
%                           - Velocity              [3x1]           [km/s]
%   P0                     Initial covariance       [6x6]           [km, km/s]
%   et0, etf               Initial and final times  [1x1]           [s]
%
% OUTPUT:
%   xx_mean                Propagated mean vector   [6x1]           [km, km/s]  
%   P                      Propagated covariance    [6x6]           [km, km/s]
% 
% AUTHOR:
% Davide Lanza
%

xxSTM0 = [xx0_mean; reshape(eye(6), [36 1])] ; % Initial state & STM vector 

options = odeset('AbsTol', 1e-13, 'RelTol', 1e-13) ;
tspan = [et0 etf] ; 

[~, xxSTM] = ode113(@(t,xxSTM) TBP_STM_ode(t,xxSTM), tspan, xxSTM0, options) ;

xx_mean = xxSTM(end, 1:6)' ; 


STM = reshape(xxSTM(end,7:42), [6 6]);

P = STM * P0 * STM' ; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xx_mean, P, Sigma_pnt] = UT(Sigma_pnt0, et0, etf)
%
% [xx_mean, P, Sigma_pnt] = UT(Sigma_pnt0, et0, etf)
%
% Uncertainty propagation using unscented general transform 
%
% INPUT:
%   Sigma_pnt0             Non-Propagated points    [6x1]               
%                           - Position              [3x1]           [km]
%                           - Velocity              [3x1]           [km/s]
%   P0                     Initial covariance       [6x6]           [km, km/s]
%   et0, etf               Initial and final times  [1x1]           [s]
%
% OUTPUT:
%   xx_mean                Propagated mean vector   [6x1]           [km, km/s]  
%   P                      Propagated covariance    [6x6]           [km, km/s]
%   Sigma_pnt              Propagated sigma points             
% 
% AUTHOR:
% Davide Lanza
%

% Initial Sigma points composition 
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

options = odeset('AbsTol', 1e-13, 'RelTol', 1e-13) ;
% Sigma points propagation
Y = zeros(n,(2*n)+1) ; 
for i = 1 : (2*n)+1
    [~, Y_iesimo] = ode113(@(t, xx) TBP_ode(t, xx), [et0 etf], Sigma_pnt0(:,i), options) ; 
    Y(:,i) = Y_iesimo(end, :) ; 
end

% Weighted sample mean computation 
y_hat = zeros(n, 1) ; 
for i = 1 : (2*n + 1)
    y_hat = W_mean(i) * Y(:,i) + y_hat ; 
end

% Weighted sample covariance computation 
P_y = zeros(n,n) ; 
for i = 1 : (2*n + 1) 
    P_y = W_cov(i) * ((Y(:,i) - y_hat)*(Y(:,i) - y_hat)') + P_y ; 
end

xx_mean = y_hat ; 
P = P_y ; 

Sigma_pnt = Y ; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [xx_mean, P, R] = MC_sim(R0, et0, etf) 
%
% [xx_mean, P, R] = MC_sim(R0, et0, etf)
%
% Uncertainty propagation with Monte Carlo simulation
%
% INPUT:
%   R0            Initial multivariate normal random generator  [nx6]  
%   et0, etf      Initial and final times                       [1x1]     [s]
%
% OUTPUT:
%   xx_mean       Propagated mean vector                        [6x1] [km, km/s]
%   P             Propagated covariance                         [6x6] [km, km/s]
%   R             Propagated samples                            [nx6] [km, km/s]
%
% AUTHOR:
% Davide Lanza

n = length(R0) ; % Number of samples 

options = odeset('AbsTol', 1e-13, 'RelTol', 1e-13) ;

% Samples propagation 
Y = zeros(6,n) ; 
for i = 1 : n
    [~, Yiesimo] = ode113(@(t, xx) TBP_ode(t, xx), [et0 etf], R0(i,:), options) ; 
    Y(:,i) = Yiesimo(end, :) ; 
end

% Sample mean
xx_mean = zeros(6,1) ; 
for i = 1 : n
xx_mean = Y(:,i) + xx_mean ;
end
xx_mean = (1/n) * xx_mean ; 

P = zeros(6,6) ; 
% Sample covariance 
for i = 1 : n
P = (Y(:,i) - xx_mean) * (Y(:,i) - xx_mean)' + P ;
end
P = (1/(n-1)) * P ; 

R = Y' ; 

end




