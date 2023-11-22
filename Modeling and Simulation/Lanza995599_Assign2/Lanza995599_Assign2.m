% Modeling and Simulation of Aerospace Systems (2022/2023)
% Assignment # 2
% Author: Davide Lanza

%% EXERCISE 1
clearvars; close all; clc; 

%%% Exercise 1.2)
par = [1 1] ; % Guessed k and b parameters 
tspan = 0 : 0.01 : 10 ; % Time vector of integration
state0 = [0 0 0 0] ; % Initial state 

J1 = 0.2 ; % Inertia disk 1 [kg*m]
J2 = 0.1 ; % Inertia disk 2 [kg*m] 

% Linearization of the non linear equations
% Jacobian at [0 0 0 0] 

J = [0 1 0 0; -par(1)/J1 0 par(1)/J1 0; 0 0 0 1; par(1)/J2 0 -par(2)/J2 par(2)/J2*0*2] ; 

lambdaJ = eig(J) ; % Non-stiff problem --> ode45 is a proper integrator 

% Integration of the dynamic 
[tt, state] = ode45(@(t, state) ode_ex1(t, state, par), tspan, state0) ; 

% Angular accelerations computation
om_dot = acceleration(state, par) ; 

figure()
% Plot of rotations
subplot(3,1,1)
plot(tt, state(:,1), 'LineWidth', 1.5) % theta1
hold on; grid on; 
plot(tt, state(:,3) ,'--','LineWidth', 2.5) % theta2 

% Plot settings
set(gca, 'FontSize', 10) ;  
xlabel('time [s]', 'Fontsize', 15) ; 
ylabel('$\theta \; [rad]$', 'Interpreter','latex', 'Fontsize', 15) ; 
legend('$\theta_1$','$\theta_2$', 'Interpreter','latex', 'FontSize', 15); 

% Plot of angular velocities 
subplot(3,1,2)
plot(tt, state(:,2), 'LineWidth', 1.5) % om1 
hold on; grid on; 
plot(tt, state(:,4), 'LineWidth', 1.5) % om2 

% Plot settings
set(gca, 'FontSize', 10) ;  
xlabel('time [s]', 'Fontsize', 15) ; 
ylabel('$\dot \theta \; [\frac{rad}{s}]$', 'Interpreter','latex', 'Fontsize', 15) ;
legend('$\dot \theta_1$','$\dot \theta_2$', 'Interpreter','latex', 'FontSize', 15); 

% Plot of angular accelerations 
subplot(3,1,3)
plot(tt, om_dot(:,1), 'LineWidth', 1.5)
hold on; grid on; 
plot(tt, om_dot(:,2), 'LineWidth', 1.5)

% Plot settings 
set(gca, 'FontSize', 10) ;  
xlabel('time [s]', 'Fontsize', 15) ; 
ylabel('$\ddot \theta \; [\frac{rad}{s^2}]$', 'Interpreter','latex', 'Fontsize', 15) ;
legend('$\ddot \theta_1$','$\ddot \theta_2$', 'Interpreter','latex', 'FontSize', 15); 

%%% Exercise 1.3)

% Importing sampes from the txt file 
import = importdata("samples.txt") ;
samples = import.data ;

% Uncontrained minimization problem 
options = optimoptions('fminunc','Display','off') ; 
par0 = [3, 3] ; 
[par_corr, fval, flag] = fminunc(@retrace_data, par0, options) ;

fprintf(['--> EXERCISE 1 \n\n ' ...
    '   The fitting parameters are: \n']) ; 
fprintf('       k = %f \n', par_corr(1)) ; 
fprintf('       b = %f \n', par_corr(2)) ; 

% Integrating the problem with the corrected parameters 
[tt, state_corr] = ode45(@(t, state) ode_ex1(t, state, par_corr), tspan, state0) ; 

% Computin the acceleration 
om_dot = acceleration(state_corr, par_corr) ; 

% Plot of both the samples and the recomputed acceleration 1,2
figure()
plot(tt, samples(:,2),'k -.' ,'LineWidth', 2.5) ; 
hold on ; grid on ; 
plot(tt, om_dot(:,1),'b','LineWidth', 1) ; 
plot(tt(2:end), samples(2:end,3),'k -.' ,'LineWidth', 2.5) ; 
plot(tt(2:end), om_dot(2:end,2),'r', 'LineWidth', 1) ; 

% Plot settings
set(gca, 'FontSize', 15) ;  
xlabel('time [s]', 'Fontsize', 20) ; 
ylabel('$\ddot \theta \; [\frac{rad}{s^2}]$', 'Interpreter','latex', 'Fontsize', 20) ;
legend('$\ddot \theta_{samp,1}$','$\ddot \theta_1$', '$\ddot \theta_{samp,2}$','$\ddot \theta_2$','Interpreter','latex', 'FontSize', 20); 

% Computing the absolute error 
err_abs1 = abs(samples(:,2)-om_dot(:,1)) ; 
err_abs2 = abs(samples(:,3)-om_dot(:,2)) ; 

% Relative error computation 
% rel_err1 = (err_abs1(2:end)./abs(samples(2:end,2))) * 100 ; 
% rel_err2 = (err_abs2(2:end)./abs(samples(2:end,3))) * 100 ;

figure()
% Error acceleration disk 1
subplot(2,1,1)
plot(tt, err_abs1, 'LineWidth', 1.5) ; 
hold on ; grid on; 
set(gca, 'YScale', 'log')

set(gca, 'FontSize', 10) ;  
xlabel('time [s]', 'Fontsize', 15) ; 
ylabel('error $[\frac{rad}{s^2}]$', 'Interpreter','latex','Fontsize', 15) ;
legend('$|\ddot \theta_{samp,1}-\ddot \theta_1|$', 'Interpreter','latex', 'FontSize', 17);

% Error acceleration disk 2
subplot(2,1,2)
plot(tt, err_abs2, 'r','LineWidth', 1.5) ; 
grid on ; grid on ; 
set(gca, 'YScale', 'log')

% Plot settings
set(gca, 'FontSize', 10) ;  
xlabel('time [s]', 'Fontsize', 15) ; 
ylabel('error $[\frac{rad}{s^2}]$', 'Interpreter','latex','Fontsize', 15) ;
legend('$$|\ddot \theta_{samp,2}-\ddot \theta_2|$$', 'Interpreter','latex', 'FontSize', 17);  


%% EXERCISE 2
clearvars; close all; clc;

% Data
R1 = 100 ; % Resistance [Ohm]
R2_k = 10 ; % Resistence, corrent-dependent [Ohm/A]
L = 10 ; % Inductance [H]
C = 1e-3 ; % Capacitor [F]
f = 5 ; % Frequency [Hz]

% Packing the parameters into a vector  
par = [R1, R2_k, L, C, f] ; 

% Integration time and inizial condition 
tspan = 0 : 0.01 : 5 ; % Time vector of itegration 
state0 = [1, 0] ; % Initial state (V0, Vdot0)

%%% Exercise 2.1)  VOLTAGE SOURCE OFF
% Integration of the equations without the voltage source
[tt, state] = ode45(@(tt, state) ode_ex2(tt, state, par,'off'), tspan, state0) ; 

% Plot of the voltage V_c
figure(1)
plot(tt, state(:,1),'LineWidth', 1.5) ; 

% Plot settings
set(gca, 'FontSize', 15) ;  grid on  ; 
xlabel('time [s]', 'Fontsize', 20) ; 
ylabel('$V_c \;\; [V] $', 'Interpreter','latex', 'Fontsize', 20) ;
legend('Voltage $V_c$','Interpreter','latex', 'FontSize', 20); 


% Plot of the derivative of the voltage V_c_dot 
figure(2)
plot(tt, state(:,2),'LineWidth', 1.5) ; grid on ; 

% Plot settings
set(gca, 'FontSize', 15) ;  
xlabel('time [s]', 'Fontsize', 20) ; 
ylabel('$\dot V_c \;\; [\frac{V}{s}] $', 'Interpreter','latex', 'Fontsize', 20) ;
x = legend('Voltage derivative $\dot V_c$','Interpreter','latex', 'FontSize', 20, 'Location','southeast'); 


%%% Exercise 2.2) VOLTAGE SOURCE ON

% Integration of the equations with the voltage source ON
[tt, state2] = ode45(@(tt, state) ode_ex2(tt, state, par,'on'), tspan, state0) ; 

% Plot of the voltage V_c
figure(3)
plot(tt, state2(:,1),'LineWidth', 1.5) ; 

% Plot settings
set(gca, 'FontSize', 15) ;  grid on ; 
xlabel('time [s]', 'Fontsize', 20) ; 
ylabel('$V_c \;\; [V] $', 'Interpreter','latex', 'Fontsize', 20) ;
legend('Voltage $V_c$','Interpreter','latex', 'FontSize', 20); 
 

%% EXERCISE 3
clearvars; close all; clc;

% Material selection:
%   Inner lining: A286 Alloy
%   Conductor: 651 Silicon Bronze alloy
%   Interface: Zirconia
%   Insulator: Zirconia oxide
%   Outer coating: Inconel Alloy 718

% Flat plate hypotesis: A = 300x300 mm
A_e = 0.1 ; % m^2

% Layers properties (Units: [m], [W/mK], [J/kgK], [kg/m^3])
l1 = 0.003 ;    k1 = 13 ;
l2 = 0.03 ;     k2 = 57 ;   c2 = 380 ;  rho2 = 8900 ;  
l3 = 0.0002 ;   k3 = 2.7 ; 
l4 = 0.008 ;    k4 = 2 ;    c4 = 460 ;  rho4 = 5700 ;
l5 = 0.003 ;    k5 = 11.4 ; 

% Outer temperature 
T_outer = 293.15 ; % [K]
T_initial = 293.15 ; % Thermal equilibrium at t = 0 [K]

% Resistances 
R1 = l1/(k1*A_e) ;     R1_h = R1/2 ; % First layer resistance/half resistance 
R2 = l2/(k2*A_e) ;     R2_h = R2/2 ; % Second layer resistance/half resistance
R3 = l3/(k3*A_e) ;     R3_h = R3/2 ; % Third layer resistance/half resistance
R4 = l4/(k4*A_e) ;     R4_h = R4/2 ; % Fourth layer resistance/half resistance
R5 = l4/(k5*A_e) ;     R5_h = R5/2 ; % Fifth layer resistance/half resistance

%%% Exercise 3.3) 1 NODE MODEL

% Equivalent resistances (1N model)
Req1_1N = l1/(k1*A_e) + l2/(2*k2*A_e) ; % Before first capacitor 
Req2_1N = l2/(2*k2*A_e) + l3/(k3*A_e) + l4/(2*k4*A_e) ; % Within the two capacitors
Req3_1N = l4/(2*k4*A_e) + l5/(k5*A_e) ; % After second capacitor 

% Integration conditions 
tspan = 0 : 0.1 : 60 ; % Time integration vector 
T0_1N = T_initial*ones(2,1) ; % Initial temperature conditions 

% It describes the change in time of the temperature in the two nodes with the capacitors
[t, T] = ode45(@(t, T) ode_1N(t, T, T_outer, A_e, Req1_1N, Req2_1N, Req3_1N, rho2, rho4, c2, c4, l2, l4), tspan, T0_1N) ; 

% Computing inner wal temperatures vector 
T_inner = Inner_Temp(t) ; 

% Second layer - Conductor 
T2 = T(:,1) ; % From integration results

% Fourth layer - Insulator 
T4 = T(:,2) ; % From integration results

% First layer - Inner lining 
T1 = T_inner - R1_h * (T_inner - T2)/(R1+ R2_h) ; 

% Third layer - Interface
T3 = T4 + (R4_h+R3_h)*(T2-T4)/(R2_h + R3 + R4_h) ; 
 
% Fifth layer - Outer coating 
T5 = T_outer + (T4-T_outer)*R5_h / (R4_h + R5) ; 

% Computing outer wal temperatures vector 
T_out = 293.15*ones(length(t),1) ; % External temperature 

% Plot of temperature profile across sections (1 NODE)
figure()
plot(t, T_inner,'-- r', 'LineWidth',1.5) ; hold on ; grid on ; 
plot(t, T1, 'LineWidth',1.5) ; 
plot(t, T2, 'LineWidth',1.5) ; 
plot(t, T3, 'LineWidth',1.5) ;
plot(t, T4, 'LineWidth',1.5) ; 
plot(t, T5, 'LineWidth',1.5) ;
plot(t,T_out,'-- b ', 'LineWidth',1.5) ;

% Plot settings (1 NODE)
set(gca, 'FontSize', 15) ;  
legend('$T_{inner}(t)$','$T_1(t)$', '$T_2(t)$','$T_3(t)$','$T_4(t)$', '$T_5(t)$','$T_{outer}(t)$','Interpreter','latex', 'FontSize', 15) ;
xlabel ('Time  [s]','FontSize',15) ;
ylabel('Temperature  [K]','FontSize',15)

% Temperature computation at sections interface 
T12 = T_inner - R1*(T_inner - T2)/(R1+R2_h) ; 
T23 = T2 - R2_h*(T2-T4)/(R2_h+R3+R4_h) ; 
T34 = T4 + R4_h*(T2-T4)/(R2_h+R3+R4_h) ; 
T45 = T4 - R4_h*(T4-T_out)/(R4_h+R5) ; 

figure()
for i = [[1, 2, 5, 7, 9,15, 50, 300, 601]]
    plot([0,3, 3, 18, 18, 33, 33, 33.2, 33.2, 37.2, 37.2, 41.2, 41.2, 44.2 ], ...
         [T_inner(i),T12(i), T12(i), T2(i), T2(i), T23(i), T23(i), T34(i), T34(i), ...
         T4(i), T4(i), T45(i), T45(i), T_out(i)], '-o', 'LineWidth', 1.5);
    hold on ; grid on ;
end

set(gca, 'FontSize', 15) ;  
legend('$T(t = 0 s)$','$T(t = 0.1 s)$','$T(t = 0.4 s)$','$T(t = 0.6 s)$','$T(t = 0.8 s)$', ...
        '$T(t = 1.4 s)$','$T(t = 4.9 s)$','$T(t = 29.9 s)$','$T(t = 60 s)$','Interpreter','latex', 'FontSize', 15) ;
xlabel ('x [mm]','FontSize',15) ;
ylabel('Temperature  [K]','FontSize',15)


%%% Exercise 3.4) 2 NODES MODEL

% One-third resistance 
R2_th = R2/3 ; 
R4_th = R4/3 ; 

% Equivalent resistances (2N model) 
Req1_2N = l1/(k1*A_e) + l2/(3*k2*A_e) ; 
Req2_2N = l2/(3*k2*A_e) + l3/(k3*A_e) + l4/(3*k4*A_e) ; 
Req3_3N = l4/(3*k4*A_e) + l5/(k5*A_e) ; 

T0_2N =  T_initial*ones(4,1) ; % Initial temperature conditions

% Integrating the equations
[t, T_2N] = ode45(@(t, T) ode_2N(t, T, T_outer, A_e, Req1_2N, Req2_2N, Req3_3N, ...
                                 R2_th, R4_th, rho2, rho4, c2, c4, l2, l4) , tspan, T0_2N) ; 

% Unpacking the integration solutions

% Second layer - Conductor 
T2_a = T_2N(:,1) ; % From integration results
T2_b = T_2N(:,2) ;  

% Fourth layer - Insulator 
T4_a = T_2N(:,3) ; % From integration results
T4_b = T_2N(:,4) ; 

% First layer - Inner lining 
T1_2N = T_inner - R1_h * (T_inner - T2_a)/(R1 + R2_th) ; 

% Third layer - Interface
T3_2N = T4_a + (R3_h + R4_th) * (T2_b - T4_a)/(R2_th+R3+R4_th) ;
 
% Fifth layer - Outer coating 
T5_2N = T_out + R5_h * (T4_b - T_out)/(R4_th+R5) ; 

% Plot of temperature profile across sections (2 NODES)
figure()
plot(t, T_inner,'-- r', 'LineWidth',1.5) ; hold on ; grid on ; 
plot(t, T1_2N, 'LineWidth',1.5) ; 
plot(t, T2_a, 'LineWidth',1.5) ; 
plot(t, T2_b, 'LineWidth',1.5) ; 
plot(t, T3_2N, 'LineWidth',1.5) ;
plot(t, T4_a, 'LineWidth',1.5) ; 
plot(t, T4_b, 'LineWidth',1.5) ; 
plot(t, T5_2N, 'LineWidth',1.5) ;
plot(t,T_out,'-- b ', 'LineWidth',1.5) ;

% Plot settings (2 NODES)
set(gca, 'FontSize', 15) ;  
legend('$T_{inner}(t)$','$T_1(t)$', '$T_{2_1}(t)$','$T_{2_2}(t)$','$T_3(t)$', ...
       '$T_{4_1}(t)$','$T_{4_2}(t)$', '$T_5(t)$','$T_{outer}(t)$','Interpreter','latex', 'FontSize', 15) ;
xlabel ('Time  [s]','FontSize',15) ;
ylabel('Temperature  [K]','FontSize',15) ; 

% Temperature computation at sections interface 
T12_2N = T_inner - R1*(T_inner - T2_a)/(R1+R2_th) ; 
T23_2N = T2_b - R2_th*(T2_b-T4_a)/(R2_th+R3+R4_th) ; 
T34_2N = T4_b + R4_th*(T2_b-T4_a)/(R2_th+R3+R4_th) ; 
T45_2N = T4_b - R4_th*(T4_b-T_out)/(R4_th+R5) ; 

figure()
for i = [[1, 2, 5, 7, 9,15, 50, 300, 601]]
    plot([0,3,13,23,33,33.2, 35.86667, 38.5333, 41.2,44.2 ], ...
         [T_inner(i),T12_2N(i), T2_a(i), T2_b(i), T23_2N(i), T34_2N(i), ...
          T4_a(i), T4_b(i), T45_2N(i), T_out(i)], '-o', 'LineWidth', 1.5);
    hold on ; grid on ;
end

set(gca, 'FontSize', 15) ;  
legend('$T(t = 0 s)$','$T(t = 0.1 s)$','$T(t = 0.4 s)$','$T(t = 0.6 s)$','$T(t = 0.8 s)$', ...
        '$T(t = 1.4 s)$','$T(t = 4.9 s)$','$T(t = 29.9 s)$','$T(t = 60 s)$','Interpreter','latex', 'FontSize', 15) ;
xlabel ('x [mm]','FontSize',15) ;
ylabel('Temperature  [K]','FontSize',15)


%% EXERCISE 4
clearvars; close all; clc;

% Data
K_m = 20 ; % [N*m*rad/A]
R = 200 ; % Resistance 
L = 2e-3 ; % Inductance [H]
J1 = 0.5 ; % Inertia disk 1 [kg*m^2]
J2 = 0.3 ; % Inertia disk 2 [kg*m^2]
b = 0.1 ; % damper coefficient [kg*m^2/s]
k = 0.5 ; % spring param [N*m]

%%% Exericse 4.2)
% System matrix
A_e = [  -R/L      0  -K_m/L     0      0
          0      0       1     0      0
     K_m/J1  -k/J1   -b/J1  k/J1   b/J1
          0      0       0     0      1
          0   k/J2     b/J2 -k/J2 -b/J2 ] ; 


% Sysem eigenvalues
lambda = eig(A_e) ; % Stiff systeml

% Plotting eigenvalues
figure(1)
for i = 1 : 5
    plot(real(lambda(i)), imag(lambda(i)), 'b s', 'MarkerSize', 13, 'LineWidth',2) 
    hold on ; grid on ; 
end

% Plot settings
set(gca, 'XScale', 'log')
ha = gca;
ha.XAxisLocation = 'origin'; % Centering the reference axis at 0 
ha.YAxisLocation = 'origin';
set(gca, 'FontSize', 15) ; grid on ; 
legend('$\lambda_i$', 'Interpreter','latex', 'FontSize', 25); 
posx = xlabel('$Re(\lambda)$', 'Interpreter','latex', 'Fontsize', 20) ;
posx.Position = [-1e-9 -1.8 -1] ; 
posy = ylabel('$Im(\lambda)$', 'Interpreter','latex', 'Fontsize', 20) ; 
posy.Position = [-1e7 1.430511474942442e-06 -1] ; 
% title('System Eigenvalues', 'FontSize', 20)

%%% Exercise 4.4)
% Integration
tspan = 0 : 0.1 : 30 ; 
state0 = [0 0 0 0 0] ; 
[tt, xx] = ode15s(@(t, xx) ode_ex4(t, xx, A_e), tspan, state0) ; 

% Current plot
figure()
plot(tt, xx(:,1), LineWidth=1.5) ; 
set(gca, 'FontSize', 15) ; grid on ; 
legend('$i(t)$', 'Interpreter','latex', 'FontSize', 25); 
xlabel('t [s]', 'Interpreter','latex', 'Fontsize', 20) ; 
ylabel('i(t) [A]', 'Interpreter','latex', 'Fontsize', 20) ; 
% title('Current system response', 'FontSize', 20)

% Rotation plot
figure()
plot(tt, xx(:,2), LineWidth=1.5) ; 
hold on; grid on ; 
plot(tt, xx(:,4), LineWidth=1.5) ; 
set(gca, 'FontSize', 15) ; 
legend('$\theta_1$','$\theta_2$', 'Interpreter','latex', 'FontSize', 25); 
xlabel('$t \; [s]$', 'Interpreter','latex', 'Fontsize', 20) ; 
ylabel('$\theta \; [rad]$', 'Interpreter','latex', 'Fontsize', 20) ; 
% title('Angular position system response', 'FontSize', 20)

% Angular velocity plot
figure()
plot(tt, xx(:,3), LineWidth=1.5) ; 
hold on; grid on ; 
plot(tt, xx(:,5), LineWidth=1.5) ; 
set(gca, 'FontSize', 15) ; 
legend('$\dot {\theta_1}$','$\dot\theta_2$', 'Interpreter','latex', 'FontSize', 25); 
xlabel('$t \; [s]$', 'Interpreter','latex', 'Fontsize', 20) ; 
ylabel('$\dot \theta \; [rad/s]$', 'Interpreter','latex', 'Fontsize', 20) ;
% title('Angular velocity system response', 'FontSize', 20)

%%% Exercise 4.5)
options = optimoptions('fminunc','Display','off') ; 
par0 = [20, 200] ;  % Initial guessed parameters (K_m and R)

% Minimization problem given the objective function 
[par_corr, fval, exitflag] = fminunc(@matching_fcn, par0, options) ;

fprintf(['--> EXERCISE 4 \n\n ' ...
    '   The fitting parameters are: \n']) ; 
fprintf('       K_m = %f \n', par_corr(1)) ; 
fprintf('       R = %f \n', par_corr(2)) ;

% Extracting the resulting parameters
K_m = par_corr(1) ; 
R = par_corr(2) ; 

% Recomputing the system matrix 
A_corr = [   -R/L     0  -K_m/L     0      0
                0     0       1     0      0
           K_m/J1 -k/J1   -b/J1  k/J1   b/J1
                0     0       0     0      1
                0  k/J2    b/J2 -k/J2  -b/J2 ] ; 

% Final integration of the problem with the fitting parameters
[tt, xx_corr] = ode15s(@(t, xx) ode_ex4(t, xx, A_corr), tspan, state0) ; 

profile = importdata("Profile.txt") ;

% Absolute error plot
figure()
plot(profile(:,1), abs(profile(:,2)-xx_corr(:,5)), 'LineWidth',1.5) ; 
grid on ; 
set(gca, 'FontSize', 15) ; set(gca, 'YScale', 'log')
legend('$|\dot \theta_2 - \dot \theta_{prof}|$', 'Interpreter','latex', 'FontSize', 25); 
xlabel('$t \; [s]$', 'Interpreter','latex', 'Fontsize', 20) ; 
ylabel('$ error \; [rad/s]$', 'Interpreter','latex', 'Fontsize', 20) ;
% title('Integration results vs Profile ', 'FontSize', 20)

% Plot to compare the two curves 
figure()
plot(tt, profile(:,2), 'k -.' ,'LineWidth', 2.5) ; 
hold on; 
plot(tt, xx_corr(:,5),'r','LineWidth', 1) ; 

% Plot settings
set(gca, 'FontSize', 15) ;  
xlabel('time [s]', 'Fontsize', 20) ; 
ylabel('$\dot \theta \; [\frac{rad}{s}]$', 'Interpreter','latex', 'Fontsize', 20) ;
legend('$\dot \theta_{prof}$','$\dot \theta_2$', 'Interpreter','latex', 'FontSize', 20); 


%% EXERCISE 5
clearvars; close all; clc;

% Heat exchanger data to compute the matrix 
L_e = 0.5 ; % Length of the pipe inside the exchanger [m]
A_e = 0.1 ; % Exchange area [m^2]
rho = 1000 ; % Density of the incompressible fluid [kg/m^3]
c_w = 4186 ;  % Specifc heat [J/kgK]
D = 0.02 ; % Tube section [m^2]
V = pi/4 * D^2 * L_e ; % Volume of the pipe inside the heat exchanger [m^3]
l1 = 0.01 ; k1 = 395; % First layer thermal properties 
l2 = 0.025 ; k2 = 310 ; c2 = 100 ; rho2 = 8620 ; % Second layer thermal properties
l3 = 0.01 ; k3 = 125 ; % Third layer thermal properties 
h = 20 ; % Heat transfer coefficient [W/m^2 K]
k = 941 ; % Spring coefficient 
m_k = 2 ; % mass [kg]
r_k = 1 ; % friction coeff [Ns/m
% Equivalent resistance: conductive res. of layer 1 + counductive resistance of half of layer 2
Req1 = l1/(k1*A_e) + l2/(2*k2*A_e) ; 
% Equivalent resistance: Convective resistance + conductive res. of layer 2 + counductive resistance of half of layer 2
Req2 = 1/(h*A_e) + l3/(k3*A_e) + l2/(2*k2*A_e) ; 


% A matrix of the associayed homogenous system
A = zeros(4,4) ; 
A11 = -1/(Req1*rho2*A_e*l2*c2) -1/(Req2*rho2*A_e*l2*c2) ; 
A12 = 1/(Req2*rho2*A_e*l2*c2) ; 
A21 = 1/(Req2*c_w*rho*V) ; 
A22 = - 1/(Req2*c_w*rho*V) ; 

A = [A11 A22      0        0;
     A21 A22      0        0 ;
     0     0      0        1 ;
     0     0 -k/m_k -r_k/m_k ] ; 

lambda = eig(A) ; % --> non stiff homogenous system
% Integration conditions 
tspan = [0 25] ; 
state0 = [320, 293.15, 0, 0] ; % Settings xc and vc at zero 

% Integration of the problem 
[tt, state] = ode45(@(t, state) ode_ex05(t, state), tspan, state0) ; 

for kk = 1:length(tt)
        [~,parout(kk,:)] = ode_ex05(tt(kk),state(kk,:));
end

Q1 = parout(:, end-1) ; % Outgoing flow rate(from tank)
Q9 = parout(:, end) ; % Incoming flow rate 

Vdot_T = -Q1 + Q9 ; % Volume derivative of the tank [m^3 / s]

% Temperatures plot (Just the two state variables)
figure()
plot(tt, state(:,2),'LineWidth', 1.5) ; 
hold on; grid on ; 
plot(tt, state(:,1),'LineWidth', 1.5) ; 

% Plot settings
set(gca, 'FontSize', 15) ;  
legend('$T_2$','$T_{fluid}$','Interpreter','latex', 'FontSize', 20) ;
xlabel ('Time  [s]','FontSize',15) ;
ylabel('Temperature  [K]','FontSize',15)

% Piston stroke plot
figure()
plot(tt, state(:,3),'LineWidth', 1.5) ;  
grid on ; 

% Plot settings
set(gca, 'FontSize', 15) ;  
xlabel('time [s]', 'Fontsize', 15) ; 
ylabel('$x_k \; [m]$', 'Interpreter','latex', 'Fontsize', 22) ; 
legend('Piston stroke', 'FontSize', 15);

% Piston speed plot
figure()
plot(tt, state(:,4),'LineWidth', 1.5) ; 
grid on ; 

% Plot settings
set(gca, 'FontSize', 15) ;  
xlabel('time [s]', 'Fontsize', 15) ; 
ylabel('$v_k \; [\frac{m}{s}]$', 'Interpreter','latex', 'Fontsize', 22) ; 
legend('Piston speed', 'FontSize', 15);

% Plotting the pressure cascade
figure()
plot(tt, parout(:,1:9)/10^6,'LineWidth', 1.5) ;
% Plot settings
set(gca, 'FontSize', 15) ;  
xlabel('time [s]', 'Fontsize', 15) ; 
ylabel('pressure [MPa]','FontSize', 20) ; 
legend('$P_1$','$P_2$','$P_3$','$P_4$','$P_5$','$P_6$','$P_7$','$P_8$','$P_9$','Interpreter','latex', 'FontSize', 18) ;

% Sublotting the pressures based on the order of magnitude (to zoom in)
figure()
subplot(3,1,1)
plot(tt, parout(:,[2 3 4])/10^6,'LineWidth', 1.5) ;
hold on; grid on ; 
set(gca, 'FontSize', 15) ;  
xlabel('time [s]', 'Fontsize', 15) ; 
ylabel('P [MPa]','FontSize', 15) ;  
legend('$P_2$','$P_3$','$P_4$','Interpreter','latex', 'FontSize', 15) ;

subplot(3,1,2)
plot(tt, parout(:,[1 9])/10^6,'LineWidth', 1.5) ;
grid on ; 
set(gca, 'FontSize', 15) ;  
xlabel('time [s]', 'Fontsize', 15) ; 
ylabel('P [MPa]','FontSize', 15) ;  
legend('$P_1$','$P_9$','Interpreter','latex', 'FontSize', 15) ;

subplot(3,1,3)
plot(tt, parout(:,[5 6 7 8])/10^6,'LineWidth', 1.5) ;
grid on ; 
set(gca, 'FontSize', 15) ;  
xlabel('time [s]', 'Fontsize', 15) ; 
ylabel('P [MPa]','FontSize', 15) ; 
legend('$P_5$','$P_6$','$P_7$','$P_8$','Interpreter','latex', 'FontSize', 15) ;

% Volume flow rate plot
figure()
plot(tt, parout(:, [12 13]),'LineWidth', 1.5) ;
grid on ; 
 
% Plot settings
set(gca, 'FontSize', 15) ;  
xlabel('time [s]', 'Fontsize', 15) ; 
ylabel('Q  $[\frac{m^3}{s}]$','Interpreter','latex', 'FontSize', 20) ; 
legend('$Q_1$','$Q_7$','Interpreter','latex', 'FontSize', 20) ;

% Volume flow rate from the tank 
figure()
plot(tt, Vdot_T, 'LineWidth',1.5) ; 
grid on ;

% Plot settings
set(gca, 'FontSize', 15) ;  
xlabel('time [s]', 'Fontsize', 15) ; 
ylabel('$\dot V_T [\frac{m^3}{s}]$','Interpreter','latex', 'FontSize', 20) ; 
legend('Fluid volume variation','FontSize', 20) ;

%% --------------------------------------------------------------------------------------------%
%------------------------------------------ FUNCTIONS------------------------------------------%
%----------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXERCISE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function state_dot = ode_ex1(t, state, par)
%
% ode_ex1.m - Ordinary differential equations of the system
%             described in exercise 1 fir the NON LINEAR damper model
%
% PROTOTYPE:
%  state_dot = ode_ex1(t, state, par)
% 
% INPUT:
%  t        [1x1]       Time of integration                         [s]               
%  state    [4x1]       State vector:
%                           theta1     = angular position disk 1    [rad]     
%                           theta1_dot = angular velocity disk 1    [rad/s]
%                           theta2     = angular position disk 2    [rad]
%                           theta2_dot = angular velocity disk 2    [rad/s]
%  par      [2x1]       Parameters of the system: 
%                           k = stiffness                           [N*m/rad] 
%                           b = viscous friction                    [N*m/(rad/s)]
% 
% OUTPUT:
%  state_dot [4x1]      Derivatives of the state:
%                           om1     = angular velocity disk 1       [rad/s] 
%                           om1_dot = angular acceleration disk 1   [rad/s^2]
%                           om2     = angular velocity disk 2       [rad/s]
%                           om2_dot = angular acceleration disk 2   [rad/s^2] 
% 
% AUTHOR:
%  Lanza Davide
% 
% VERSION:
%  25-11-2022
% 

% Inertia disks [kg*m]
J1 = 0.2 ; 
J2 = 0.1 ;  

T = 0.1 ;  %  Exteral toruge [N*m] 

k = par(1) ; % Torisional spring constant
b = par(2) ; % Torsional damping coefficient

% Output vector initialization
state_dot = zeros(4,1) ; % 

% Angular velocities 
om1 = state(2) ; % first disk 
om2 = state(4) ; % second disk 

% Getting the four derivatives of the state 
state_dot(1) = om1 ; % theta1_dot

state_dot(2) = k/J1 * (state(3) - state(1))  ; % om1_dot

state_dot(3) = om2 ;  % theta2_dot

state_dot(4) = 1/J2 * (-k*(state(3)-state(1)) - sign(om2)*b*om2^2 + T) ; % om2_dot 


end

%----------------------------------------------------------------------------------------------%


function [state_dot] = ode_ex1_Lin(t, state, par)
%
% ode_ex1.m - Ordinary differential equations of the system
%             described in exercise 1 for the LINEAR damper model
%
% PROTOTYPE:
%  state_dot = ode_ex1(t, state, par)
% 
% INPUT:
%  t        [1x1]       Time of integration                         [s]               
%  state    [4x1]       State vector:
%                           theta1     = angular position disk 1    [rad]     
%                           theta1_dot = angular velocity disk 1    [rad/s]
%                           theta2     = angular position disk 2    [rad]
%                           theta2_dot = angular velocity disk 2    [rad/s]
%  par      [2x1]       Parameters of the system: 
%                           k = stiffness                           [N*m/rad] 
%                           b = viscous friction                    [N*m/(rad/s)]
% 
% OUTPUT:
%  state_dot [4x1]      Derivatives of the state:
%                           om1     = angular velocity disk 1       [rad/s] 
%                           om1_dot = angular acceleration disk 1   [rad/s^2]
%                           om2     = angular velocity disk 2       [rad/s]
%                           om2_dot = angular acceleration disk 2   [rad/s^2] 
% 
% AUTHOR:
%  Lanza Davide
% 
% VERSION:
%  25-11-2022
% 

% Inertia disks [kg*m]
J1 = 0.2 ; 
J2 = 0.1 ; 

k = par(1) ; 
b = par(2) ; 

T = 0.1 ;  %  Exteral toruge [N*m] 

state_dot = zeros(4,1) ; % Initialization

% Vector 
bb = [0 0 0 T]' ; 

% Matrix that describes the linear system 
A = [    0 1     0     0
     -k/J1 0  k/J2     0
         0 0     0     1
      k/J2 0 -k/J2 -b/J2] ; 

% Packing the output: xdot = Ax + b
state_dot = A*state + bb ;

end

%----------------------------------------------------------------------------------------------%


function om_dot = acceleration(state, par)
%
% acceleration.m - computes angular accelerations of both disks one and two
%
% PROTOTYPE:
%  om_dot = acceleration(state, par)
% 
% INPUT:               
%  state    [nx4]       State vector:
%                           theta1     = angular position disk 1    [rad]     
%                           theta1_dot = angular velocity disk 1    [rad/s]
%                           theta2     = angular position disk 2    [rad]
%                           theta2_dot = angular velocity disk 2    [rad/s]
%  par      [2x1]       Parameters of the system: 
%                           k = stiffness                           [N*m/rad] 
%                           b = viscous friction                    [N*m/(rad/s)]
% 
% OUTPUT:
%  om_dot  [nx2]        Angular accelerations
%                           om1_dot = angular acceleration disk 1   [rad/s^2]
%                           om2_dot = angular acceleration disk 2   [rad/s^2] 
% 
% AUTHOR:
%  Lanza Davide
% 
% VERSION:
%  25-11-2022
%

% Inertia disks [kg*m]
J1 = 0.2 ; 
J2 = 0.1 ;  

T = 0.1 ;  %  Exteral toruge [N*m] 

k = par(1) ; % Torisional spring constant
b = par(2) ; % Torsional damping coefficient

om_dot = zeros(length(state), 2) ; % Inizialization 

om_dot(:, 1) = k/J1 * (state(:,3) - state(:,1))  ; % om1_dot

om_dot(:, 2) = 1/J2 * (-k*(state(:,3)-state(:,1)) - sign(state(:,4))*b.*state(:,4).^2 + T) ; % om2_dot 

end

%----------------------------------------------------------------------------------------------%

function obj = retrace_data(par)
%
% retrace_data.m - Objective function to be given to the fminunc.m function
%                  in order to minimize the stated problem
%
% PROTOTYPE:
%  obj = retrace_data(par)
% 
% INPUT:
%  par      [2x1]       Parameters of the system: 
%                           k = stiffness                           [N*m/rad] 
%                           b = viscous friction                    [N*m/(rad/s)]
% 
% OUTPUT:
%  obj      [1x1]       Sum of the errors                            [rad/s^2]
% 
% AUTHOR:
%  Lanza Davide
% 
% VERSION:
%  25-11-2022
%

% Inertia disks [kg*m]
J1 = 0.2 ; 
J2 = 0.1 ;  

T = 0.1 ;  %  Exteral toruge [N*m] 

tspan = 0 : 0.01 : 10 ; % Vector of times
state0 = [0 0 0 0] ;  % Initial guessed state

% Importing samples from txt file
import = importdata("samples.txt") ;
samples = import.data ;

% Integrating the ODEs
[tt, state] = ode45(@(t, state) ode_ex1(t, state, par), tspan, state0) ; 

% Computing the accelerations
om_dot = acceleration(state, par) ;

% Computing the two error vectors
err_abs = [abs(om_dot(:,1) - samples(:,2)), abs(om_dot(:,2) - samples(:,3))] ; 

% Summing the sum of the two error vectors
obj = sum(sum(err_abs)) ; % Function to be minimezed 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXERCISE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function state_dot = ode_ex2(t, state, par, volt_source)
%
% ode_ex2.m - Ordinary differential equations of the system
%             described in exercise 2. 
%
% PROTOTYPE:
%  state_dot = ode_ex2(t, state, par, volt_source)
% 
% INPUT:
%  t        [1x1]       Time of integration                         [s]               
%  state    [2x1]       State vector:
%                           V_c     = Voltage across the capacitor  [V]     
%                           V_c_dot = Voltage derivative "  "       [V/s]
%  par      [5x1]       Parameters of the system: 
%                           R1   = Resistance                       [Ohm]
%                           R2_k = Resistence, i-dep                [Ohm/A]
%                           L    = Inductance                       [H]
%                           C    = Capacitor                        [F]
%                           f    = Frequency                        [Hz]
% 
% OUTPUT:
%  state_dot [2x1]      Derivatives of the state:
%                           V_c_dot = Voltage derivative            [V/s] 
%                           V_c_dotdot = double volt derivative     [V/s^2] 
% 
% AUTHOR:
%  Lanza Davide
% 
% VERSION:
%  25-11-2022
% 

% Unpacking the parameters vector 
R1 = par(1) ; % Resistance [Ohm]
R2_k = par(2) ; % Resistence, corrent-dependent [Ohm/A]
L = par(3) ; % Inductance [H]
C = par(4) ; % Capacitor [F]
f = par(5) ; % Frequency [Hz]

% Voltage source on-off cases
switch volt_source
    case 'on'
        v_dot = @(t) cos(2*pi*f*t)*2*pi*f*atan(t) + (1/(1+t^2))*sin(2*pi*f*t) ; 
    
    case 'off' 
        v_dot = @(t) 0 ; 
end

% Output inizialization 
state_dot = zeros(2,1) ; 

state_dot(1) = state(2) ; % Vc_dot 

% numerator and denominator to visually simplify the expression 
denom = -R1*C - 2*R2_k*C^2*state(2) ;
numer = state(2) + R1*state(1)/L + R1*R2_k*C^2*state(2)^2/L - v_dot(t) ; 

state_dot(2) = numer/denom ; % Vc_dotdot 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXERCISE 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T_dot = ode_1N(t, T, T_out, A, Req1, Req2, Req3, rho2, rho4, c2, c4, l2, l4)
%
% ode_1N.m - Ordinary differential equations of the system
%            described in exercise 3 with one node per layer positioned 
%            in the center. 
%
% PROTOTYPE:
%  T_dot = ode_1N(t, T, T_out, A, Req1, Req2, Req3, rho2, rho4, c2, c4, l2, l4)
% 
% INPUT:
%  t        [1x1]       Time of integration                     [s]               
%  T        [2x1]       State vector:
%                           T2 = Temperature at node 2          [K]   
%                           T4 = Temperature at node 4          [K]
%  T_out    [1x1]       Outer temperature                       [K]
%  A        [1x1]       Cross section                           [m^2]
%  Req_i    [1x1]       Equivalent resistences                  [K/W] 
%  rho_i    [1x1]       Material density                        [kg/m^3]
%  c_i      [1x1]       Heat capacity                           [J/kgK]
%  l_i      [1x1]       length                                  [m]
%
% OUTPUT:
%  T_dot [2x1]          Derivatives of the state:
%                           T2_dot = Temp. deriv. node 2        [K/s] 
%                           T4_dot = Temp. deriv. node 4        [K/s]
% 
% AUTHOR:
%  Lanza Davide
% 
% VERSION:
%  25-11-2022
% 

% Output initialization
T_dot = zeros(2,1) ; 

% Building the ramp inner temperature
if t >= 0 && t <= 1 
    T_inner = 293.15 + 980*t ; 
else
    T_inner = 1273.15 ; 
end

% Computing output derivatives 
T_dot(1) = (Req2*(T_inner - T(1)) - Req1*(T(1)-T(2))) / (Req1*Req2*rho2*A*l2*c2) ; % T2_dot

T_dot(2) = (Req3*(T(1) - T(2)) - Req2*(T(2)-T_out)) / (Req2*Req3*rho4*A*l4*c4) ; % T4_dot 

end

%----------------------------------------------------------------------------------------------%


function T_dot = ode_2N(t, T, T_out, A, Req1, Req2, Req3, R2_th, R4_th, rho2, rho4, c2, c4, l2, l4)
%
% ode_2N.m - Ordinary differential equations of the system
%            described in exercise 3 with one node per layer positioned 
%            in the center and two nodes for layer 2 & 4 positioned at
%            one third and two third of the length.
%
% PROTOTYPE:
%  T_dot = ode_2N(t, T, T_out, A, Req1, Req2, Req3, R2_th, R4_th, rho2, rho4, c2, c4, l2, l4)
% 
% INPUT:
%  t        [1x1]       Time of integration                     [s]               
%  T        [4x1]       State vector:
%                           T2_a = Temp. at node 2a             [K]  
%                           T2_b = Temp. at node 2b             [K]  
%                           T4_a = Temp. at node 4a             [K]  
%                           T4_b = Temp. at node 2a             [K]  
%  T_out    [1x1]       Outer temperature                       [K]
%  A        [1x1]       Cross section                           [m^2]
%  Req_i    [1x1]       Equivalent resistences                  [K/W] 
%  R_th_i   [1x1]       One-third of resistence                 [K/W]
%  rho_i    [1x1]       Material density                        [kg/m^3]
%  c_i      [1x1]       Heat capacity                           [J/kgK]
%  l_i      [1x1]       length                                  [m]
%
% OUTPUT:
%  T_dot [4x1]          Derivatives of the state:
%                           T2a_dot = Temp. deriv. node 2a      [K/s] 
%                           T2b_dot = Temp. deriv. node 2b      [K/s] 
%                           T4a_dot = Temp. deriv. node 4a      [K/s] 
%                           T4b_dot = Temp. deriv. node 4b      [K/s] 
% 
% AUTHOR:
%  Lanza Davide
% 
% VERSION:
%  25-11-2022
%

% Output initialization
T_dot = zeros(4,1) ; 

% Building the ramp inner temperature
if t >= 0 && t <= 1 
    T_in = 293.15 + 980*t ; 
else
    T_in = 1273.15 ; 
end

% Computing output derivatives 
T_dot(1) = (R2_th*(T_in-T(1)) - Req1 * (T(1)-T(2)))/(Req1*R2_th*rho2*A*l2*c2/2) ; % T2a_dot
T_dot(2) = (Req2*(T(1)-T(2)) - R2_th * (T(2)-T(3)))/(Req2*R2_th*rho2*A*l2*c2/2) ; % T2b_dot
T_dot(3) = (R4_th*(T(2)-T(3)) - Req2 * (T(3)-T(4)))/(Req2*R4_th*rho4*A*l4*c4/2) ; % T4a_dot
T_dot(4) = (Req3*(T(3)-T(4)) - R4_th * (T(4)-T_out))/(Req3*R4_th*rho4*A*l4*c4/2) ; % T4b_dot

end

%----------------------------------------------------------------------------------------------%

function T_in = Inner_Temp(tspan)
%
% Inner_Temp.m - Computes the inner temperature vector knowing that the
%                in the first second it ramps and then it's constant.
%
% PROTOTYPE:
%  T_in = Inner_Temp(tspan)
% 
% INPUT:
%  tspan    [nx1]       Vector of times             [s]
%
% OUTPUT:
%  T_in     [nx1]       Inner temperature vector    [K]
% 
% AUTHOR:
%  Lanza Davide
% 
% VERSION:
%  25-11-2022
%

% Initialization of the output
T_in = zeros(length(tspan), 1) ; 

% Ramp inner temperature
for i = 1 : length(tspan)
    if (tspan(i) >= 0) &&  (tspan(i) <= 1) % Ramp temperature in the first second
        T_in(i) = 293.15 + 980*tspan(i) ;
    else
        T_in(i) = 1273.15 ; % Costante temperature after one second 
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXERCISE 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xx_dot] = ode_ex4(t, xx, A)
%
% ode_ex4.m - Ordinary differential equations of the system
%             described in exercise 4 
%
% PROTOTYPE:
%  [xx_dot] = ode_ex4(t, xx, A)
% 
% INPUT:
%  t        [1x1]       Time of integration                         [s]               
%  x        [5x1]       State vector:
%                           i     = current                         [A]     
%                           theta1 = rotation disk 1                [rad]
%                           om1     = angular velocity disk 1       [rad/s]
%                           theta2 = rotation disk 2                [rad]
%                           om1     = angular velocity disk 1       [rad/s]
%  A        [5x5]        System matrix
% 
% OUTPUT:
%  xx_dot  [5x1]       Derivatives of the state:
%                           i_dot     = current derivative                 [A/s] 
%                           theta1_dot = angular velocity disk 1           [rad/s]
%                           om1_dot     = angular acceleration disk 1      [rad/s^2]
%                           theta2_dot = angular veocity disk 2            [rad/s] 
%                           om2_dot = angular acceleration disk 2          [rad/s^2]
% 
% AUTHOR:
%  Lanza Davide
% 
% VERSION:
%  28-11-2022
%

v0 = 2 ; % Voltage [V]
om = 5 ; % Frequency [Hz]
beta = 0.2 ; % Parameter [Hz]
L = 2e-3 ; % Inductance [H]

xx_dot = zeros(5,1) ; % Initialization

% Matricial form for linear systems: xdot = Ax + b
b = [1/L * (v0*cos(om*t)*exp(-beta*t)); 0; 0; 0; 0] ; 

xx_dot = A*xx + b ; 

% extended form 

% xx_dot(1) = 1/L * (v0*cos(om*t)*exp(-beta*t) - R*xx(1) + K_m*xx(3)) ; % di/dt
% 
% xx_dot(2) = xx(3) ; % dtheta1/dt
% 
% xx_dot(3) = 1/J1* (K_m*xx(1) + k*(xx(4) - xx(2)) + b*(xx(5) - xx(3))) ; %dOM1/dt
% 
% xx_dot(4) = xx(5) ; % dtheta2/dt 
% 
% xx_dot(5) = 1/J2 * (-k*(xx(4) - xx(2)) - b*(xx(5) - xx(3))) ; %dOM2/dt

end

%----------------------------------------------------------------------------------------------%

function [obj] = matching_fcn(par)
%
% matching_fcn.m - Objective function to be given to the fminunc.m function
%                  in order to minimize the stated problem
%
% PROTOTYPE:
%  [obj] = matching_fcn(par)
% 
% INPUT:
%  par      [2x1]       Parameters of the system: 
%                           K_m = stiffness                           [N*m/rad] 
%                           R  = viscous friction                    [N*m/(rad/s)]
% 
% OUTPUT:
%  obj      [1x1]       Sum of the errors                            [rad/s^2]
% 
% AUTHOR:
%  Lanza Davide
% 
% VERSION:
%  28-11-2022
%
profile = importdata("Profile.txt") ;

% Given data

K_m = par(1) ; % Parameter 
R = par(2) ; % Resistance 
v0 = 2 ; % Voltage [V]
om = 5 ; % Frequency [Hz]
beta = 0.2 ; % Parameter [Hz]
L = 2e-3 ; % Inductance [H]

J1 = 0.5 ; % Inertia disk 1 [kgm^2]
J2 = 0.3 ; % Inertia disk 2 [kgm^2]
b = 0.1 ; % Dumper parameter [kgm^2/s]
k = 0.5 ; % Spring coefficient [Nm]

% Matrix describing the system 
A = [   -R/L      0  -K_m/L     0     0
           0      0      1     0     0
      K_m/J1  -k/J1  -b/J1  k/J1  b/J1
           0      0     0     0      1
           0   k/J2  b/J2 -k/J2  -b/J2 ] ;

% Integrating the ODEs
[tt, xx] = ode15s(@(t, xx) ode_ex4(t, xx, A), [0 : 0.1 : 30], [0 0 0 0 0]) ; 

% Computing the absolute error vector
err = abs(xx(:,5) - profile(:,2)) ; 

% Summing the errors to minimize a scalar 
obj = sum(err) ; 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXERCISE 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [state_dot, parout] = ode_ex05(t, state)
%
% ode_ex05.m - Ordinary differential equations of the system
%             described in exercise 5 containing other parameters
%
% PROTOTYPE:
%  [state_dot, parout] = ode_ex05(t, state)
% 
% INPUT:
%  t        [1x1]       Time of integration                         [s]               
%  state    [4x1]       State vector:
%                           T_2     = Temperature secondo layer     [K]     
%                           V_fluid = Temperature of the fluid      [K]
%                           x_k = piston stroke                     [m]
%                           v_k = piston speed                      [m/s]
%
% OUTPUT
%  par      [13x1]       Parameters of the system: 
%                           P_i = Pressure at i                     [Pa]
%                           Q_dot = Heat flow                       [W] 
%                           Q_i = volumetric flow rate              [m^3/s]
%                           
%
% AUTHOR:
%  Lanza Davide
% 
% VERSION:
%  3-12-2022
%

%%%%%%%%%%%%
%%% DATA %%%
%%%%%%%%%%%%

% Fluid 
rho = 1000 ; % Density of the incompressible fluid [kg/m^3]
c_w = 4186 ;  % Specifc heat [J/kgK]

% Lines
k_cv = 2 ; % Pressure drop coeff across the check valve [-]
D = 0.02 ; % Diameter of the lines [m]

% Length ([m]) and friction factor ([-]) of the branches
L_T1 = 0.5 ; f_T1 = 0.032 ; 
L_34 = 1.5 ; f_34 = 0.032 ; 
L_56 = 2.7 ; f_56 = 0.040 ; 
L_78 = 2.5 ; f_78 = 0.028 ; 
L_9T = 1 ;   f_9T = 0.032 ; 

% Tank
P_T = 0.1e6 ; % Constat tank pressure [Pa]

% Pitstons
N = 9 ; % Number of pistons
D_p = 0.007 ; % Diameter [m]
A_p = pi*D_p^2/4 ; % Area of a piston [m^2]
d_p = 0.015 ; % Distance shaft [m]
P_nom = 506625 ; % Nominal pressure [Pa]


% Pilot piston:
D_k = 0.01 ; % Diameter [m]
A_k = pi * D_k^2/4 ; % Area [m^2]
l_c = 0.1 ; % Control lever length [m]
theta_max = 20*pi/180 ; % Maximum angle of the control plate [rad]
n = 4000/60 ; % Rotation speed [rps]
m_k = 2 ; % Equivalent mass [kg]
F0 = 5 ; % Pre-loaded force [N]
r_k = 1 ; % Friction coefficient [Ns/m]
d_k = 0.001 ; % Diameter of the pipe [m]
A_k_pipe = pi*d_k^2/4 ; % Section of the pipe [m^2]
k_p = 2.5 ; % Pilot pipe head loss [-]

% Distributor 
k_d = 15 ; % Pressure drop coeff across distributor [-]
d0 = 0.010 ; % Diameter [m]

% Filter 
k_f = 35 ; % Pressure drop coeff across filter [-]
k_l = 2.5/100 ; % Leaking coefficient 

% Heat exchanger 
L_e = 0.5 ; % Length of the pipe inside the exchanger [m]
A = 0.1 ; % Exchange area [m^2]

l1 = 0.01 ; k1 = 395; % First layer thermal properties 
l2 = 0.025 ; k2 = 310 ; c2 = 100 ; rho2 = 8620 ; % Second layer thermal properties
l3 = 0.01 ; k3 = 125 ; % Third layer thermal properties 
 
h = 20 ; % Heat transfer coefficient [W/m^2 K]
T0 = 350; % Constant [K]
k_T = 20; % Coefficient [-]
om = 5 ; % Frequency [Hz]
k = 1000 ; % [W]


%%%%%%%%%%%%%%%%%%%%%%%%%%% HEAT EXCHANGER ODES %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Equivalent resistance: conductive res. of layer 1 + counductive resistance of half of layer 2
Req1 = l1/(k1*A) + l2/(2*k2*A) ; 

% Equivalent resistance: Convective resistance + conductive res. of layer 2 + counductive resistance of half of layer 2
Req2 = 1/(h*A) + l3/(k3*A) + l2/(2*k2*A) ; 

T = T0 + k_T*cos(om*t) ; 
state_dot = zeros(4,1) ; % Output initialization 

% Derivative of the state vector (Temperature at node 2 and at node of the fluid) 
state_dot(1) = (T - state(1))/(Req1*rho2*A*l2*c2) - (state(1)-state(2))/(Req2*rho2*A*l2*c2) ; % T2_dot 
state_dot(2) = (state(1) - state(2))/(Req2*c_w*rho*pi/4 * D^2 * L_e) ; % T_fluid_dot

Qdot = (state(1) - state(2))/Req2 ;   % Heat exchanged [W]


%%%%%%%%%%%%%%%%%%%%%%%%% PISTON STROKE AND SPEED %%%%%%%%%%%%%%%%%%%%%%%%%


A = pi*D^2/4 ; % Cross section of the tube lines [m^2]

x = state(3); % Piston stroke
v = state(4); % Piston speed 

c_max = l_c*tan(theta_max);  

% Physical constraint related to the piston
if x < 0
    x = 0;
end
if x <= 0 && v < 0
    x = 0;
    v = 0;
end
if x > c_max
    x = c_max;
end
if x >= c_max && v > 0
    x = c_max;
    v = 0;
end

s = d_p/l_c * (c_max - x) ; % displacement
K = (P_nom*A_k - F0)/c_max ; % Spring coefficient 

Q1 = n * N * A_p * s ; % Flow rate from 1 to 6
Q7 = (1-k_l)*Q1 ; % Flow rate from 7 to T

% Pressure cascade computation 
P9 = P_T + (1/2 * f_9T * L_9T / D * rho * Q7*abs(Q7)/(A^2)) ;
P8 = P9/(exp(Qdot/k)) ;
P7 = P8 + (1/2 * f_78 * L_78 / D * rho * Q7*abs(Q7)/(A^2)) ;
P6 = P7 + (1/2 * k_f * rho * Q1^2/(A^2)) ; 
P5 = P6 + (1/2 * f_56 * L_56 / D * rho * Q1*abs(Q1)/(A^2)) ;

r0 = d0/2 ; % Radius 

% A_d(z) variation over time which is opening radially in 2 second
if t <2
    z = r0*t/2 ; 
    alpha = 2*acos(z/r0) ; 
    a = r0^2/2*(alpha-sin(alpha)) ; 
    A_d = pi*r0^2 - a ; 
else
    A_d = pi*r0^2 ; 
end

% Pressure cascade computation
P4 = P5 + (1/2 * k_d * rho * Q1*abs(Q1)/(A_d)^2) ; 
P3 = P4 + (1/2 * f_34 * L_34 / D * rho * Q1*abs(Q1)/(A^2)) ; 
P2 = P3 + (1/2 * k_cv * rho * Q1*abs(Q1)/(A^2)) ; 
P_k = P2 - 1/2 * k_p * rho * (v*A_k/A_k_pipe)*abs(v*A_k/A_k_pipe) ;
P1 = P_T - (1/2 * f_T1 * L_T1 / D * rho * Q7*abs(Q1)/(A^2)) ;

% Packagin the output 
state_dot(3) = v;
state_dot(4) = (P_k*A_k - K*x - F0 - r_k*v)/m_k;

% Other parameters that will be used out of the ode function 
parout = [P1 P2 P3 P4 P5 P6 P7 P8 P9 P_k Qdot Q1 Q7];

end





