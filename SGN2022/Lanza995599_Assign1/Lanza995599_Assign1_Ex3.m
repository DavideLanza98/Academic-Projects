% Spacecraft Guidance and Navigation (2022/2023)
% Assignment # 1
% Exercise # 3
% Author: Davide Lanza

%% EXERCISE 3.2: FOUR THRUSTERS CONFIGURARTION

% Clear workspace, close all figures, and clear command window
clearvars; close all; clc;

% Clear kernels
cspice_kclear(); 

% loading SPICE kernels:
cspice_furnsh('assignment01.tm') ; 

% Launch date and initial state 
t0 = '2022-08-03-12:45:20.000 UTC' ; 
et0 = cspice_str2et(t0) ; % Transforming in ephemeris time [s] 
xx0 = cspice_spkezr('Earth',et0,'ECLIPJ2000','NONE','SSB') ; % Initial state vector [km], [km/s]
rr0 = xx0(1:3) ; % Extracting initial position  [km]                                                   
vv0 = xx0(4:6); % Extracting initial velocity  [km/s]                                                      

% Other data
m0 = 1500 ; % Spacecraft mass [kg]
Tmax = 4*0.15*10^-3; % 4 Thruster * maximum thrust  [kN]
Isp = 3000 ; % Specific impulse [s]
muS = cspice_bodvrd('Sun','GM',1) ; % Sun gravitational parameter [km^3/s^2]
g0 = 9.8067 * 1e-3 ; % Gravitational constant [km/s^2]  

% Adimensionalization (scaling)                                                           
AU = 149597870.7 ; % Astronomical unit [km]                                                   
DU = 1*AU ; % Distance unit [km]                                                            
VU = sqrt(muS/DU) ; % Velocity unit [km/s]                                                   
TU  = DU/VU ; % Time unit [s]                                                            

par.muS = muS*(TU^2/DU^3) ; % Nondimensional Sun parameter[-]
par.rr0 = rr0/DU ; % Nondimensional initial position [-]                                                        
par.vv0 = vv0/DU*TU;  % Nondimensional initial velocity [-]    
par.Tmax = Tmax/(DU/TU^2);  % Nondimensional thrust [-]   
par.Isp = Isp/TU ; % Nondimensional specific impulse [-]
par.g0 = g0*(TU^2/DU) ; % Nondimensional gravit. constant [-]  
par.m0 = m0 ; % Nondimensional mass [-]

par.et0 = et0/TU ; % Nondimensional initial time [-]

% TPBVP Solver: Initial Guess Generation and Shooting Function Root Search

% Loop until error is below tolerance
% error = 1;
% error_tol = 1e-9; % Set error tolerance for exit criteria
% 
% while error > error_tol
% 
%     % Initial guess for the fsolve problem 
%     etf_rand = par.et0 + ((250 + 100*rand(1)) * cspice_spd) / TU;
%     lambda0_rand = zeros(7,1); % resetting 
%     lambda0_rand(1:6) = 10.^(2*(2*rand(6,1)-1)).*(2*rand(6,1)-1);
%     lambda0_rand(7) = 10^(2*(2*rand(1)-1));
% 
%     % Finding the roots
%     options = optimset('Algorithm', 'levenberg-marquardt', 'Display', ...
%                         'iter', 'MaxIter', 150, 'MaxFunEvals', 800); 
%     Y_opt = fsolve(@(Y) shootingFun(Y, par), [lambda0_rand; etf_rand], options);
%     optimal_lambda0 = Y_opt(1:7);
%     optimal_etf = Y_opt(8);
% 
%     % Check of the loop 
%     error = norm(shootingFun(Y_opt, par));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: Comment out the lines below and uncomment the loop above if
%          you want to check the code execution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optimal initial guesses from commented loop (from a random run)
optimal_etf = 1.471852999097268e2 ;

optimal_lambda0 = [ -54.986852462677860 ; 
                    39.040322506207210  ; 
                    1.764377616436101   ; 
                    -67.903937561702340 ; 
                    -27.376777801571176 ; 
                    -0.072505475878433  ; 
                    0.005333942883151 ] ;

% Time of Flight [days]
ToF = (optimal_etf - par.et0)*TU/cspice_spd;

% Final state error with respect to Mars 
residual = shootingFun([optimal_lambda0; optimal_etf], par);

err_rr = norm(residual(1:3)) * DU ; % Final potision error[km]
err_vv = norm(residual(4:6)) * DU/TU*1e3 ;  % Final velocity error [m/s]

% Final time string
Final_date = datestr(seconds(optimal_etf*TU),'dd-mmm-yy HH:MM:SS') ; 

% Display results
fprintf('Solutions of the four-thruster TPBVP are:\n\n')
fprintf('   λ = [');
fprintf('%f ',optimal_lambda0(1:7))
fprintf(']''\n')
fprintf('   Arrival date: %s UTC \n',Final_date)
fprintf('   Time of flight: %f days\n', ToF)
fprintf('   Position error: %f km\n', err_rr)
fprintf('   Velocity error: %e m/s\n', err_vv)
fprintf('------------------------------------------------------------------------------------\n')

% Integrating state and the costate, given the initial condition found 
options = odeset('AbsTol',1e-13 , 'RelTol', 1e-13); 
tspan = [par.et0 optimal_etf] ; 
Y0 = [par.rr0; par.vv0 ; par.m0; optimal_lambda0] ; 
[time_vec , Y] = ode113(@(t,Y) TPBVP_ode(t,Y,par),tspan , Y0, options);

% Propagation in the usefull tspan
tspan_E = TU*(par.et0 : 0.1 : optimal_etf) ; 
tspan_M = TU*(par.et0 : 0.1 : optimal_etf) ; 

rr_E_prop = cspice_spkpos('Earth', tspan_E,'ECLIPJ2000', 'NONE', 'SSB'); 
rr_M_prop = cspice_spkpos('Mars', tspan_M,'ECLIPJ2000', 'NONE', 'SSB'); 

% Plotting orbits in 2D
figure(1)
plot3(Y(:,1), Y(:,2), Y(:,3), 'color', '#A3A3A3', 'LineWidth',2)
hold on ; grid on ; 
plot3(rr_E_prop(1,:)/DU, rr_E_prop(2,:)/DU, rr_E_prop(3,:)/DU, '--', 'LineWidth', 1.5, 'Color', '#ADD8E6')
plot3(rr_M_prop(1,:)/DU, rr_M_prop(2,:)/DU, rr_M_prop(3,:)/DU, '--', 'LineWidth', 1.5, 'Color', '#FFA07A') 

scatter3(0, 0, 0, 'filled', 'SizeData', 140,'MarkerFaceColor','#FFD700') % Sun
scatter3(rr_E_prop(1,1)/DU, rr_E_prop(2,1)/DU, rr_E_prop(3,1)/DU, 'filled', 'SizeData', 140,'MarkerFaceColor','#0072BD') % Earth 
scatter3(rr_M_prop(1, end)/DU, rr_M_prop(2, end)/DU, rr_M_prop(3, end)/DU, 'filled', 'SizeData', 140,'MarkerFaceColor','#D62728') % Mars
view(2)

% Plot settings
set(gca,'FontSize',12)
legend('Transfer orbit','Earth orbit','Mars orbit', 'Sun','Departure from Earth','Arrival at Mars', 'FontSize',14) 
xlabel('$x$ [-]','Interpreter','latex','FontSize', 20)
ylabel('$y$ [-]','Interpreter','latex','FontSize', 20)
zlabel('$z$ [-]','Interpreter','latex','FontSize', 20)
title('Time optimal Earth-Mars transfer orbit (@Sun ECLIPTIC J2000)', 'Interpreter', 'Latex', 'fontsize', 18)

% Computing the thrust pointing unit vector  
alpha_cap = zeros(length(time_vec), 3) ; % Allocation 
for i=1:size(Y,1)
    alpha_cap(i,:) = -Y(i,11:13)./norm(Y(i,11:13)) ; 
end

% Plotting orbits in 3D with the pointing thrust unity vector 
figure(2)
plot3(Y(:,1), Y(:,2), Y(:,3),'color', '#A3A3A3', 'LineWidth',2)
hold on ; grid on ; 
plot3(rr_E_prop(1,:)/DU, rr_E_prop(2,:)/DU, rr_E_prop(3,:)/DU, '--', 'LineWidth',1.5, 'Color', '#ADD8E6')
plot3(rr_M_prop(1,:)/DU, rr_M_prop(2,:)/DU, rr_M_prop(3,:)/DU, '--', 'LineWidth',1.5, 'Color', '#FFA07A') 

scatter3(0, 0, 0, 'filled', 'SizeData', 140,'MarkerFaceColor','#FFD700') % Sun
scatter3(rr_E_prop(1,1)/DU, rr_E_prop(2,1)/DU, rr_E_prop(3,1)/DU, 'filled', ...
         'SizeData', 140,'MarkerFaceColor','#0072BD') % Earth 
scatter3(rr_M_prop(1, end)/DU, rr_M_prop(2, end)/DU, rr_M_prop(3, end)/DU, 'filled', ...
         'SizeData', 140,'MarkerFaceColor','#D62728') % Mars

quiver3(Y(1:2:end,1), Y(1:2:end,2), Y(1:2:end,3), alpha_cap(1:2:end,1), ...
    alpha_cap(1:2:end,2), alpha_cap(1:2:end,3), 'Color', '#1E90FF', 'LineWidth', 0.7);

daspect([1,1,0.15]); % Scaling the z axis differently
view(3)

% Plot settings
set(gca,'FontSize',12)
legend('Transfer orbit','Earth orbit','Mars orbit', 'Sun','Departure from Earth', ...
        'Arrival at Mars','Thrust pointing vector','Location','northeast', 'FontSize',14) 
xlabel('$x$ [-]', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('$y$ [-]', 'Interpreter', 'latex', 'FontSize', 20)
zlabel('$z$ [-]', 'Interpreter', 'latex', 'FontSize', 20)
% title('Thrust pointing unit vector (@Sun ECLIPTIC J2000)', 'Interpreter', 'Latex', 'fontsize', 18)

% Hamiltonian computation and plot
H_time = zeros(size(Y,1),1) ; 
for i = 1 : size(Y,1)

H_time(i) = 1 + dot(Y(i,8:10),Y(i,4:6)) -  (par.muS/norm(Y(i,1:3))^3) * dot(Y(i,11:13), ...
            Y(i,1:3)) + (par.Tmax/(par.Isp*par.g0)) * (-norm(Y(i,11:13)) * ...
            par.Isp*par.g0/Y(i,7) - Y(i,14)) ;
end

figure(3)
plot((time_vec-par.et0)*TU/cspice_spd, H_time, 'LineWidth', 2) ; 
grid on ; axis equal 

set(gca,'FontSize',12)
xlabel('time [days]','FontSize',20, 'Interpreter', 'Latex')
ylabel('H [-]','FontSize', 20, 'Interpreter', 'Latex')
legend('Hamiltonian','FontSize',14)
title('Time Evolution of the Hamiltonian along the Optimal Trajectory', 'Interpreter', 'Latex', 'fontsize', 18)

% Computation and plot of thrust angles (in-plane and out-of-plane angles)
phi = zeros(length(time_vec),1); % In-plane angle 
theta = zeros(length(time_vec),1); % Out-of-plane angle

for i = 1 : length(time_vec)
    phi(i) = rad2deg(acos(alpha_cap(i,1)/sqrt(alpha_cap(i,1)^2 + alpha_cap(i,2)^2)));   % In-plane  
    theta(i) = rad2deg(atan(alpha_cap(i,3)/sqrt(alpha_cap(i,1)^2 + alpha_cap(i,2)^2))); % Out-of-plane 

end

figure(4)
plot((time_vec-par.et0)*TU/cspice_spd, phi(:), 'LineWidth', 2)
hold on ; grid on
plot((time_vec-par.et0)*TU/cspice_spd, theta(:), 'LineWidth', 2)

set(gca,'FontSize',12)
xlim([0,ToF])
xlabel('time [days]', 'Interpreter', 'Latex', 'fontsize', 20)
ylabel('$\phi$, $\theta$ [deg]', 'Interpreter', 'latex', 'FontSize', 20)
legend('In-Plane Angle $\phi$', 'Out-of-Plane Angle $\theta$', 'Interpreter', 'Latex', 'fontsize', 14)
title('Thrust Angles Profile', 'Interpreter', 'Latex', 'fontsize', 18)

% Switching function computation and plot
S_t = zeros(length(time_vec),1) ; 
for i = 1 : length(time_vec)
    S_t(i) = -norm(Y(i, 11:13))*par.Isp*par.g0/Y(i,7) - Y(i,14) ; 
end

figure(5)
plot((time_vec-par.et0)*TU/cspice_spd, S_t, LineWidth=2)
grid on ;  

set(gca,'FontSize',12)
xlabel('time [days]','FontSize',20, 'Interpreter', 'Latex')
ylabel('$S_t$ [-]','FontSize', 20, 'Interpreter', 'Latex')
legend('Switching function','FontSize',14)
title('Time Evolution of the switching function', 'Interpreter', 'Latex', 'fontsize', 18)
xlim([0,ToF]) 


%% EXERCISE 3.3: THREE THRUSTERS CONFIGURARTION

% Clear workspace, close all figures, and clear command window
clearvars; close all; clc;

% Clear kernels
cspice_kclear(); 

% loading SPICE kernels:
cspice_furnsh('assignment01.tm') ; 

% Launch date and initial state 
t0 = '2022-08-03-12:45:20.000 UTC' ; 
et0 = cspice_str2et(t0) ; % Transforming in ephemeris time [s] 
xx0 = cspice_spkezr('Earth',et0,'ECLIPJ2000','NONE','SSB') ; % Initial state vector [km], [km/s]
rr0 = xx0(1:3) ; % Extracting initial position  [km]                                                   
vv0 = xx0(4:6); % Extracting initial velocity  [km/s]                                                      

% Other data
m0 = 1500 ; % Spacecraft mass [kg]
Tmax = 3*0.15*10^-3; % 4 Thruster * maximum thrust  [kN]
Isp = 3000 ; % Specific impulse [s]
muS = cspice_bodvrd('Sun','GM',1) ; % Sun gravitational parameter [km^3/s^2]
g0 = 9.8067 * 1e-3 ; % Gravitational constant [km/s^2]  

% Adimensionalization (scaling)                                                           
AU = 149597870.7 ; % Astronomical unit [km]                                                   
DU = 1*AU ; % Distance unit [km]                                                            
VU = sqrt(muS/DU) ; % Velocity unit [km/s]                                                   
TU  = DU/VU ; % Time unit [s]                                                            

par.muS = muS*(TU^2/DU^3) ; % Nondimensional Sun parameter[-]
par.rr0 = rr0/DU ; % Nondimensional initial position [-]                                                        
par.vv0 = vv0/DU*TU;  % Nondimensional initial velocity [-]    
par.Tmax = Tmax/(DU/TU^2);  % Nondimensional thrust [-]   
par.Isp = Isp/TU ; % Nondimensional specific impulse [-]
par.g0 = g0*(TU^2/DU) ; % Nondimensional gravit. constant [-]  
par.m0 = m0 ; % Nondimensional mass [-]

par.et0 = et0/TU ; % Nondimensional initial time [-]

% TPBVP Solver: Initial Guess Generation and Shooting Function Root Search
% Loop until error is below tolerance
% error = 1;
% error_tol = 1e-9; % Set error tolerance for exit criteria
% while error > error_tol
% 
%     % Initial guess for the fsolve problem 
%     etf_rand = par.et0 + ((350 + 300*rand(1)) * cspice_spd) / TU;
%     lambda0_rand = zeros(7,1); % resetting 
%     lambda0_rand(1:6) = 10.^(2*(2*rand(6,1)-1)).*(2*rand(6,1)-1);
%     lambda0_rand(7) = 10^(2*(2*rand(1)-1));
% 
%     % Finding the roots
%     options = optimset('Algorithm', 'levenberg-marquardt', 'Display', ...
%                         'iter', 'MaxIter', 150, 'MaxFunEvals', 1500); 
%     Y_opt = fsolve(@(Y) shootingFun(Y, par), [lambda0_rand; etf_rand], options);
%     optimal_lambda0 = Y_opt(1:7);
%     optimal_etf = Y_opt(8);
% 
%     % Check of the loop 
%     error = norm(shootingFun(Y_opt, par));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: Comment out the lines below and uncomment the loop above if
%          you want to check the code execution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optimal initial guesses from commented loop (from a random run)
optimal_etf = 1.500856472088991e2 ;

optimal_lambda0 = [ -61.1708076171535     ; 
                    54.5531170745652      ; 
                    0.557456086980271     ; 
                    -80.52850261005790    ; 
                    -45.8226867181350     ; 
                    0.602307742093763     ; 
                    0.00724183831305695 ] ;


% Time of Flight [days]
ToF = (optimal_etf - par.et0)*TU/cspice_spd;

% Final state error with respect to Mars 
residual = shootingFun([optimal_lambda0; optimal_etf], par);

err_rr = norm(residual(1:3)) * DU ; % Final potision error[km]
err_vv = norm(residual(4:6)) * DU/TU*1e3 ;  % Final velocity error [m/s]

% Final time string
Final_date = datestr(seconds(optimal_etf*TU),'dd-mmm-yy HH:MM:SS') ; 

% Display results
fprintf('Solutions of the three-thruster TPBVP are:\n\n')
fprintf('   λ = [');
fprintf('%f ',optimal_lambda0(1:7))
fprintf(']''\n')
fprintf('   Arrival date: %s UTC \n',Final_date)
fprintf('   Time of flight: %f days\n', ToF)
fprintf('   Position error: %f km\n', err_rr)
fprintf('   Velocity error: %e m/s\n', err_vv)
fprintf('------------------------------------------------------------------------------------\n')

% Integrating state and the costate, given the initial condition found 
options = odeset('AbsTol',1e-13 , 'RelTol', 1e-13); 
tspan = [par.et0 optimal_etf] ; 
Y0 = [par.rr0; par.vv0 ; par.m0;optimal_lambda0] ; 
[time_vec , Y] = ode113(@(t,Y) TPBVP_ode(t,Y,par),tspan , Y0, options);

% Propagation in the usefull tspan
tspan_E = TU*(par.et0 : 0.1 : optimal_etf) ; 
tspan_M = TU*(par.et0 : 0.1 : optimal_etf) ; 

rr_E_prop = cspice_spkpos('Earth', tspan_E,'ECLIPJ2000', 'NONE', 'SSB'); 
rr_M_prop = cspice_spkpos('Mars', tspan_M,'ECLIPJ2000', 'NONE', 'SSB'); 

% Plotting orbits in 2D
figure(1)
plot3(Y(:,1), Y(:,2), Y(:,3), 'color', '#A3A3A3', 'LineWidth',2)
hold on ; grid on ; 
plot3(rr_E_prop(1,:)/DU, rr_E_prop(2,:)/DU, rr_E_prop(3,:)/DU, '--', 'LineWidth', 1.5, 'Color', '#ADD8E6')
plot3(rr_M_prop(1,:)/DU, rr_M_prop(2,:)/DU, rr_M_prop(3,:)/DU, '--', 'LineWidth', 1.5, 'Color', '#FFA07A') 

scatter3(0, 0, 0, 'filled', 'SizeData', 140,'MarkerFaceColor','#FFD700') % Sun
scatter3(rr_E_prop(1,1)/DU, rr_E_prop(2,1)/DU, rr_E_prop(3,1)/DU, 'filled', 'SizeData', 140,'MarkerFaceColor','#0072BD') % Earth 
scatter3(rr_M_prop(1, end)/DU, rr_M_prop(2, end)/DU, rr_M_prop(3, end)/DU, 'filled', 'SizeData', 140,'MarkerFaceColor','#D62728') % Mars
view(2)

% Plot settings
set(gca,'FontSize',12)
legend('Transfer orbit','Earth orbit','Mars orbit', 'Sun','Departure from Earth','Arrival at Mars', 'FontSize',14) 
xlabel('$x$ [-]','Interpreter','latex','FontSize', 20)
ylabel('$y$ [-]','Interpreter','latex','FontSize', 20)
zlabel('$z$ [-]','Interpreter','latex','FontSize', 20)
title('Time optimal Earth-Mars transfer orbit (@Sun ECLIPTIC J2000)', 'Interpreter', 'Latex', 'fontsize', 18)

% Computing the thrust pointing unit vector  
alpha_cap = zeros(length(time_vec), 3) ; % Allocation 
for i=1:size(Y,1)
    alpha_cap(i,:) = -Y(i,11:13)./norm(Y(i,11:13)) ; 
end

% Plotting orbits in 3D with the pointing thrust unity vector 
figure(2)
plot3(Y(:,1), Y(:,2), Y(:,3),'color', '#A3A3A3', 'LineWidth',2)
hold on ; grid on ; 
plot3(rr_E_prop(1,:)/DU, rr_E_prop(2,:)/DU, rr_E_prop(3,:)/DU, '--', 'LineWidth',1.5, 'Color', '#ADD8E6')
plot3(rr_M_prop(1,:)/DU, rr_M_prop(2,:)/DU, rr_M_prop(3,:)/DU, '--', 'LineWidth',1.5, 'Color', '#FFA07A') 

scatter3(0, 0, 0, 'filled', 'SizeData', 140,'MarkerFaceColor','#FFD700') % Sun
scatter3(rr_E_prop(1,1)/DU, rr_E_prop(2,1)/DU, rr_E_prop(3,1)/DU, 'filled', ...
         'SizeData', 140,'MarkerFaceColor','#0072BD') % Earth 
scatter3(rr_M_prop(1, end)/DU, rr_M_prop(2, end)/DU, rr_M_prop(3, end)/DU, 'filled', ...
         'SizeData', 140,'MarkerFaceColor','#D62728') % Mars

quiver3(Y(1:2:end,1), Y(1:2:end,2), Y(1:2:end,3), alpha_cap(1:2:end,1), ...
    alpha_cap(1:2:end,2), alpha_cap(1:2:end,3), 'Color', '#1E90FF', 'LineWidth', 0.7);

daspect([1,1,0.15]); % Scaling the z axis differently
view(3)

% Plot settings
set(gca,'FontSize',12)
legend('Transfer orbit','Earth orbit','Mars orbit', 'Sun','Departure from Earth', ...
        'Arrival at Mars','Thrust pointing vector','Location','northeast', 'FontSize',14) 
xlabel('$x$ [-]', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('$y$ [-]', 'Interpreter', 'latex', 'FontSize', 20)
zlabel('$z$ [-]', 'Interpreter', 'latex', 'FontSize', 20)
% title('Thrust pointing unit vector (@Sun ECLIPTIC J2000)', 'Interpreter', 'Latex', 'fontsize', 18)

% Hamiltonian computation and plot
H_time = zeros(size(Y,1),1) ; 
for i = 1 : size(Y,1)

H_time(i) = 1 + dot(Y(i,8:10),Y(i,4:6)) -  (par.muS/norm(Y(i,1:3))^3) * dot(Y(i,11:13), ...
            Y(i,1:3)) + (par.Tmax/(par.Isp*par.g0)) * (-norm(Y(i,11:13)) * ...
            par.Isp*par.g0/Y(i,7) - Y(i,14)) ;
end

figure(3)
plot((time_vec-par.et0)*TU/cspice_spd, H_time, 'LineWidth', 2) ; 
grid on ; axis equal 

set(gca,'FontSize',12)
xlabel('time [days]','FontSize',20, 'Interpreter', 'Latex')
ylabel('H [-]','FontSize', 20, 'Interpreter', 'Latex')
legend('Hamiltonian','FontSize',14)
title('Time Evolution of the Hamiltonian along the Optimal Trajectory', 'Interpreter', 'Latex', 'fontsize', 18)

% Computation and plot of thrust angles (in-plane and out-of-plane angles)
phi = zeros(length(time_vec),1); % In-plane angle 
theta = zeros(length(time_vec),1); % Out-of-plane angle

for i = 1 : length(time_vec)
    phi(i) = rad2deg(acos(alpha_cap(i,1)/sqrt(alpha_cap(i,1)^2 + alpha_cap(i,2)^2)));   % In-plane  
    theta(i) = rad2deg(atan(alpha_cap(i,3)/sqrt(alpha_cap(i,1)^2 + alpha_cap(i,2)^2))); % Out-of-plane 

end

figure(4)
plot((time_vec-par.et0)*TU/cspice_spd, phi(:), 'LineWidth', 2)
hold on ; grid on
plot((time_vec-par.et0)*TU/cspice_spd, theta(:), 'LineWidth', 2)

set(gca,'FontSize',12)
xlim([0,ToF])
xlabel('time [days]', 'Interpreter', 'Latex', 'fontsize', 20)
ylabel('$\phi$, $\theta$ [deg]', 'Interpreter', 'latex', 'FontSize', 20)
legend('In-Plane Angle $\phi$', 'Out-of-Plane Angle $\theta$', 'Interpreter', 'Latex', 'fontsize', 14)
title('Thrust Angles Profile', 'Interpreter', 'Latex', 'fontsize', 18)

% Switching function computation and plot
S_t = zeros(length(time_vec),1) ; 
for i = 1 : length(time_vec)
    S_t(i) = -norm(Y(i, 11:13))*par.Isp*par.g0/Y(i,7) - Y(i,14) ; 
end

figure(5)
plot((time_vec-par.et0)*TU/cspice_spd, S_t, LineWidth=2)
grid on ;  

set(gca,'FontSize',12)
xlabel('time [days]','FontSize',20, 'Interpreter', 'Latex')
ylabel('$S_t$ [-]','FontSize', 20, 'Interpreter', 'Latex')
legend('Switching function','FontSize',14)
title('Time Evolution of the switching function', 'Interpreter', 'Latex', 'fontsize', 18)
xlim([0,ToF]) 


%% FUNCTIONS

function [Y_dot] = TPBVP_ode(~,Y,par)
%
% [Y_dot] = TPBVP_ode(~,Y,par)
% 
% It provides the odefun the optimal TPBVP
%
% INPUT
%  t                Time                                                                                         
%  Y                State and costate vector:                   [14x1]             
%                     -State:
%                             r          Position               [3x1]    
%                             v          Velocity               [3x1]    
%                             m          Mass                   [1x1]     
%                     -Costate:
%                             lambda_r   Position costate       [3x1]    
%                             lambda_v   Velocity costate       [3x1]    
%                             lambda_m   Mass costate           [1x1]   
%  par              Parameters of the problem                   [struct]
%
% OUTPUT
%  Y_dot            Derivative of state and costate                [14x1]
%
% NOTE: all the quantities are nondimensional 
%
% AUTHOR:
%  Davide Lanza
%


% Parameters
muS = par.muS ;
Tmax = par.Tmax ; 
Isp= par.Isp ; 
g0 = par.g0 ;

% Extract the state vectors
rr = Y(1:3) ;   
vv = Y(4:6) ;

m = Y(7) ; % Extract Mass

% Extract costate vectors
lambda_rr = Y(8:10) ;
lambda_vv = Y(11:13) ;

r = sqrt(rr(1)^2+rr(2)^2+rr(3)^2) ; 
lambda_vv_norm = sqrt(lambda_vv(1)^2+lambda_vv(2)^2+lambda_vv(3)^2) ; 
u_star = 1 ;  % Time-optimal throattling factor 

% State vector derivatives
dxxdt = [                                                 vv(1);
                                                          vv(2);
                                                          vv(3);
        -muS*rr(1)/r^3-(u_star*Tmax/m)*(lambda_vv(1)/lambda_vv_norm);
        -muS*rr(2)/r^3-(u_star*Tmax/m)*(lambda_vv(2)/lambda_vv_norm);
        -muS*rr(3)/r^3-(u_star*Tmax/m)*(lambda_vv(3)/lambda_vv_norm);
                                               -u_star*Tmax/(Isp*g0) ] ;
 
dlambda_rdt = -3*muS*dot(rr,lambda_vv)*rr/r^5 + muS*lambda_vv/r^3 ;
dlambda_vdt = - lambda_rr ;
dlambda_mdt = - u_star*lambda_vv_norm*Tmax/m^2 ;

% Packing out the final vector 


Y_dot =  [                                         dxxdt;
                                             dlambda_rdt ; 
                                            dlambda_vdt; 
                                            dlambda_mdt ] ;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = shootingFun(yy, par) 
%     
% F = shootingFun(yy, par)
% 
% It provides the shooting function for the time optimal 
% TPBVP for continous low-thrust guidance in the 2BP.
%
% INPUT: 
%  yy              Unknowns array:
%                     - lambda_r                [3x1]   
%                     - lambda_v                [3x1]   
%                     - lambda_m                [1x1]     
%                     - tf                      [1x1]        
% par             Problem parameters   
%
% OUTPUT:
%  F              Errors of the TPBVP           [8x1]
%
% NOTE: the problem is nondimensional
%
% Author:
%  Davide Lanza
%

% Parameters 
muS =par.muS ;
Tmax = par.Tmax ; 
Isp= par.Isp ; 
et0 = par.et0 ; 
rr0 =par.rr0 ;
vv0 = par.vv0 ; 
g0 = par.g0 ;
m0 = par.m0 ;
u_star = 1 ; % Time optimal 

% Adimenstional unit 
DU=149597870700*1e-3;
TU=5.0229*1e6; 

%Propagation
Y0=[rr0; vv0; m0; yy(1:7)] ;
options = odeset('RelTol', 2.5e-13, 'AbsTol', 2.5e-13);
[~, YY] = ode113(@(t,Y) TPBVP_ode(t, Y, par), [et0 yy(8)], Y0, options);

% Unpacking the final conditions 
YY_f = YY(end,:)' ; 
rr_f = YY_f(1:3)  ; % Position
vv_f = YY_f(4:6)  ; % Velocity 
m_f= YY_f(7) ;

lambda_rrf = YY_f(8:10);                    % [-]
lambda_vvf = YY_f(11:13);                   % [-]
lambda_mf = YY_f(14);  % [-]

% Mars State vector 
xx_M=cspice_spkezr('Mars', yy(8)*TU, 'ECLIPJ2000', 'NONE', 'SSB');

% Allocation 
rr_M = zeros(3,1) ; % Position 
vv_M = zeros(3,1) ; % Velocity 

rr_M(1:3) = xx_M(1:3)/DU ; % Position vector of Mars
vv_M(1:3) = xx_M(4:6)/DU*TU ; % Velocity vector of Mars
aa_M=-(muS/norm(rr_M)^3)*rr_M; % Acceleration vector of Mars

% Hamiltonian computation at final time 
H = 1+dot(lambda_rrf,vv_f)-muS/norm(rr_f)^3*dot(rr_f, lambda_vvf)...
    -u_star*Tmax*norm(lambda_vvf)/m_f-lambda_mf*u_star*Tmax/(Isp*g0);

% Building shooting function vector
F = [ rr_f-rr_M ;
      vv_f-vv_M  ; 
      lambda_mf ;
      H - lambda_rrf'*vv_M- lambda_vvf'*aa_M ] ; 


end


