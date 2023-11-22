% Spacecraft Guidance and Navigation (2022/2023)
% Assignment # 1
% Exercise # 2
% Author: Davide Lanza

%% EXERCISE 2.1

% Clear workspace, close all figures, and clear command window
clearvars; close all; clc;

% Clear kernels
cspice_kclear(); 

% Physical constants (from Topputo 2013)
mu = 1.21506683e-2 ;    % Earth-Moon mass parameter             [-] 
m_s = 3.28900541e5 ;    % Scaled mass of the Sun                [-]
rho = 3.88811143e2 ;    % Scaled Sun–(Earth + Moon) distance    [-]
w_s = -9.25195985e-1 ;  % Scaled angular velocity of the Sun    [-]
Re = 6378e3 ;           % Mean Earth's radius                   [m]
hi = 167e3 ;            % Altitude of departure orbit           [m]
DU = 3.84405000e8 ;     % Distance unit                         [m]

% Transformed parameters
alpha = 1.5*pi ; % Angle defined on Earth circular parking orbit 
beta = 1.41 ;    % initial-to-circular velocity ratio 
t_i = 0 ;        % Initial time 
delta = 7 ;      % Transfer duration

% Initial guess solutions 
r0 = (Re + hi)/DU ;          % [-]
v0 = beta*sqrt((1-mu)/r0) ;  % [-]

x0 = r0*cos(alpha) - mu ;    % [-]
y0 = r0*sin(alpha) ;         % [-]
v_x0 = -(v0-r0)*sin(alpha) ; % [-]
v_y0 = (v0-r0)*cos(alpha) ;  % [-] 

% Integration settings
xx0 = [x0; y0; v_x0; v_y0] ; % Initial state [4x1]
tspan = [t_i, t_i + delta] ; % Time span for integration 
par = [mu, m_s, rho, w_s] ; % Parameters

% Dynamic integration in barycentre rotating frame
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[tt, xx] = ode113(@PBRFBP_ode, tspan, xx0, options);

% Plot in the barycentre rotating frame
figure(1)
plot(xx(:,1), xx(:,2), 'k','LineWidth', 1.5)
hold on ; grid on ;

% Scattering Planets
scatter(-mu, 0, 'filled', 'SizeData', 140,'MarkerFaceColor','#1E90FF') 
scatter(1-mu, 0, 'filled', 'SizeData', 130,'MarkerFaceColor','#D3D3D3') 

set(gca,'FontSize',12)
xlabel('x [-]', 'FontSize', 20)
ylabel('y [-]', 'FontSize', 20)
legend('Transfer Orbit','Earth','Moon', 'FontSize', 18)
title('Transfer Orbit (@Earth-Moon Barycentre Rotating Frame)' , 'Interpreter', 'Latex', 'fontsize', 18)

% Conversion to the Earth-centred intertial frame (from appendix 1) 
X1 = zeros(length(xx), 1) ; 
X2 = zeros(length(xx), 1) ; 
for i = 1 : length(xx)
    X1(i) = (xx(i,1) + mu)*cos(tt(i)) - xx(i,2)*sin(tt(i)) ;
    X2(i) = (xx(i,1) + mu)*sin(tt(i)) + xx(i,2)*cos(tt(i)) ;
end

% Plot in the Earth-Centred intertial frame
figure(2)
plot(X1(:), X2(:), 'k','LineWidth', 1.5)
hold on ; grid on ; 

scatter(0, 0, 'filled', 'SizeData', 150, 'MarkerFaceColor','#1E90FF') 
scatter(1, 0, 'filled', 'SizeData', 150, 'MarkerFaceColor','#D3D3D3') 

set(gca,'FontSize',12)
xlabel('x [-]', 'FontSize', 20)
ylabel('y [-]', 'FontSize', 20)
legend('Transfer Orbit','Earth','Moon', 'FontSize', 18)
title('Transfer Orbit (@Earth-centred intertial frame)' , 'Interpreter', 'Latex', 'fontsize', 18)

%% EXERCISE 2.2

% Clear workspace, close all figures, and clear command window
clearvars; close all; clc;

% Clear kernels     
cspice_kclear(); 

% Physical constants (from Topputo 2013)
mu = 1.21506683e-2 ;    % Earth-Moon mass parameter             [-] 
w_s = -9.25195985e-1 ;  % Scaled angular velocity of the Sun    [-]
Re = 6378e3 ;           % Mean Earth's radius                   [m]
hi = 167e3 ;            % Altitude of departure orbit           [m]
DU = 3.84405000e8 ;     % Distance unit                         [m]
TU = 4.34811305 ;       % Time unit                             [days]

% Transformed parameters
alpha = 1.5*pi ;    % Angle defined on Earth circular parking orbit 
beta = 1.41 ;       % initial-to-circular velocity ratio 
delta = 7 ;         % Transfer duration
t_i = 0 ;           % Initial time 
t_f = t_i + delta ; % Final time 

% Initial guess solutions 
r0 = (Re + hi)/DU ;          % [-]
v0 = beta*sqrt((1-mu)/r0) ;  % [-]

x0 = r0*cos(alpha) - mu ;    % [-]
y0 = r0*sin(alpha) ;         % [-]
v_x0 = -(v0-r0)*sin(alpha) ; % [-]
v_y0 = (v0-r0)*cos(alpha) ;  % [-] 

% Initial NLP variable vector 
xx_i = [x0; y0; v_x0; v_y0; t_i; t_f] ; 

% Lower and Upper boundaries
v_c = sqrt((1-mu)/r0 ); % Circular velocity around Earth [-]
LB = [ (-r0-mu) -r0 (-sqrt(2)*v_c + r0) (-sqrt(2)*v_c + r0) 0 0 ];
UB = [ (r0-mu) r0 (sqrt(2)*v_c - r0) (sqrt(2)*v_c - r0) abs(2*pi/w_s) 23];

% Settings
options_NoGrad = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'active-set', ...
                           'MaxFunctionEvaluations', 1500, 'ConstraintTolerance', 1.5e-7) ;  
warning('off');

% Minimization of the Objective function without gradient provided
[x_NoGrad, fval_NoGrad, ~, out_NoGrad] = fmincon(@(xx_i) myObjFnc(xx_i), xx_i, [], [], [], [], LB, UB, @(xx_i) myNlConstr(xx_i), options_NoGrad) ; 

% Settings
options_Grad = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'active-set', ...
                           'MaxFunctionEvaluations', 1500, 'SpecifyConstraintGradient',true, ...
                           'SpecifyObjectiveGradient', true, 'ConstraintTolerance', 1.5e-7) ;  
warning('off');

% Minimization of the Objective function with gradient provided 
[x_Grad, fval_Grad, ~, out_Grad] = fmincon(@(xx_i) myObjFnc(xx_i), xx_i, [], [], [], [], LB, UB, @(xx_i) myNlConstr(xx_i), options_Grad) ; 

% Print results without gradient
fprintf('<strong>Results for No-Gradient Method:</strong>\n');
fprintf('       Total impulse: <strong>%f m/s</strong>\n', fval_NoGrad *((DU)/(TU*(3600*24))));
fprintf('       Time of flight: <strong>%f days</strong>\n',  (x_NoGrad(end) - x_NoGrad(end-1))*TU);
fprintf('       Number of iterations: <strong>%d</strong>\n', out_NoGrad.iterations);

fprintf('--------------------------------------------------\n');

% Print results with gradient
fprintf('\n<strong>Results for Gradient Method:</strong>\n');
fprintf('       Total impulse: <strong>%f m/s</strong>\n', fval_Grad *((DU)/(TU*(3600*24))));
fprintf('       Time of flight: <strong>%f days</strong>\n',  (x_Grad(end) - x_Grad(end-1))*TU);
fprintf('       Number of iterations: <strong>%d</strong>\n', out_Grad.iterations);

% Plots
% No Gradient solution 
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[~, xx_NoGrad] = ode113(@PBRFBP_ode, [x_NoGrad(5), x_NoGrad(6)], x_NoGrad(1:4), options);

figure(1) 
plot(xx_NoGrad(:,1), xx_NoGrad(:,2), 'LineWidth', 1.5)
grid on; hold on ; 

% Scattering Planets
scatter(-mu, 0, 'filled', 'SizeData', 140,'MarkerFaceColor','#1E90FF') 
scatter(1-mu, 0, 'filled', 'SizeData', 130,'MarkerFaceColor','#D3D3D3') 

% Settings
set(gca,'FontSize',12)
xlabel('x [-]', 'FontSize', 20)
ylabel('y [-]', 'FontSize', 20)
legend('Transfer Orbit' ,'Earth', 'Moon' ,'FontSize', 14)
title('Simple Shooting (a) (@Earth-Moon Barycenter Rotating frame' , 'Interpreter', 'Latex', 'fontsize', 18)

% Gradient solution
[~ , xx_Grad] = ode113(@PBRFBP_ode, [x_Grad(5), x_Grad(6)], x_Grad(1:4), options);
figure(2) 
plot(xx_Grad(:,1), xx_Grad(:,2), 'LineWidth', 1.5)
grid on; hold on ; 

% Scattering Planets
scatter(-mu, 0, 'filled', 'SizeData', 140,'MarkerFaceColor','#1E90FF') 
scatter(1-mu, 0, 'filled', 'SizeData', 130,'MarkerFaceColor','#D3D3D3') 

set(gca,'FontSize',12)
xlabel('x [-]', 'FontSize', 20)
ylabel('y [-]', 'FontSize', 20)
legend('Transfer Orbit' ,'Earth', 'Moon' ,'FontSize', 14)
title('Simple Shooting (b) (@Earth-Moon Barycenter Rotating frame' , 'Interpreter', 'Latex', 'fontsize', 18)


%% EXERCISE 2.3

% Clear workspace, close all figures, and clear command window
clearvars; close all; clc;

% Clear kernels     
cspice_kclear(); 

% Physical constants (from Topputo 2013)
mu = 1.21506683e-2 ;    % Earth-Moon mass parameter             [-] 
m_s = 3.28900541e5 ;    % Scaled mass of the Sun                [-]
rho = 3.88811143e2 ;    % Scaled Sun–(Earth + Moon) distance    [-]
w_s = -9.25195985e-1 ;  % Scaled angular velocity of the Sun    [-]
l_em = 3.84405000e8 ;   % Earth-Moon distance                   [m]
om_em = 2.66186135e-6 ; % Earth-Moon angular velocity           [s^-1]
Re = 6378e3 ;           % Mean Earth's radius                   [m]
Rm = 1738e3 ;           % Mean Moon's radius                    [m]
hi = 167e3 ;            % Altitude of departure orbit           [m]
hf = 100e3 ;            % Altitude of arrival orbit             [m]
DU = 3.84405000e8 ;     % Distance unit                         [m]
TU = 4.34811305 ;       % Time unit                             [days]
VU = 1.02323281e3 ;     % Speed unit                            [m/s]

% Transformed parameters
alpha = 1.5*pi ;    % Angle defined on Earth circular parking orbit 
beta = 1.41 ;       % initial-to-circular velocity ratio 
delta = 7 ;         % Transfer duration
t_i = 0 ;           % Initial time 
t_f = t_i + delta ; % Final time 

% Initial guess solutions 
r0 = (Re + hi)/DU ;          % [-]
v0 = beta*sqrt((1-mu)/r0) ;  % [-]

x0 = r0*cos(alpha) - mu ;    % [-]
y0 = r0*sin(alpha) ;         % [-]
v_x0 = -(v0-r0)*sin(alpha) ; % [-]
v_y0 = (v0-r0)*cos(alpha) ;  % [-] 

% Time grid definition
N = 4 ;                      % Number of points
tt_grid = zeros(1,N) ; % Allocation 
for j = 1 : N
    tt_grid(j) = t_i + ((j - 1)/(N - 1)) * (t_f - t_i) ; % Time grid of N points
end

% Initial state of each segment 
xx_i = [x0; y0; v_x0; v_y0] ;  % Initial state (@ node 1)
xx0(:,1) = xx_i ; 
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); % Settings

for j = 1 : (N - 1) % Initial states at node 2, 3 and 4

    [~, xx] = ode113(@PBRFBP_ode, [tt_grid(j), tt_grid(j+1)], xx0(:,j), options) ; % Dynamic integration
    
    xx0(:,j+1) = xx(end,:) ; % Building a matrix which contains all the initial conditions
end

% Initial NLP variable vector 
yy0 = [reshape(xx0, [16 1]); t_i; t_f] ; % Inital vector [18x1]

% Lower and Upper boundaries
v_c = sqrt((1-mu)/r0 ); % Circular velocity around Earth [-]
LB = [(-r0-mu) -r0 (-sqrt(2)*v_c + r0) (-sqrt(2)*v_c + r0) -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf 0 0];
UB = [ (r0-mu) r0 (sqrt(2)*v_c - r0) (sqrt(2)*v_c - r0) inf inf inf inf inf inf inf inf inf inf inf inf abs(2*pi/w_s) 23];

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'active-set', ...
                        'MaxFunctionEvaluations', 1500, 'ConstraintTolerance', 1.5e-7, ...
                        'SpecifyConstraintGradient', true, 'SpecifyObjectiveGradient',true) ;  

% Minimization of the Objective function
[yy, fval, flag, out] = fmincon(@(yy0) myObjFnc_MS(yy0), yy0, [], [], [], [], LB, UB, @(yy0) myNlConstr_MS(yy0), options) ;

% Print results
fprintf('\n<strong>Results for Multiple Shooting Method:</strong>\n') ;
fprintf('       Total impulse: <strong>%f m/s</strong>\n', fval*((DU)/(TU*(3600*24)))) ;
fprintf('       Time of flight: <strong>%f days</strong>\n',  (yy(18) - yy(17))*TU) ;
fprintf('       Number of iterations: <strong>%d</strong>\n', out.iterations) ;

% Plot
xx1 = yy(1:4) ; % Initial state @ node 1
xx2 = yy(5:8) ; % Initial state @ node 2
xx3 = yy(9:12) ; % Initial state @ node 3
xx4 = yy(13:16) ; % Initial state @ node 4

xx_opt = [xx1 xx2 xx3 xx4] ; 

% New time grid definition
tt_opt_grid = zeros(1, N) ; 
for j = 1 : N
    tt_opt_grid(j) = yy(end-1) + ((j - 1)/(N - 1)) * (yy(end) - yy(end-1)) ; % Time grid of N points
end

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12) ;

figure(1)
for j = 1 : N-1

    % Propagating each segment 
    [~, prop] = ode113( @PBRFBP_ode, [tt_opt_grid(j) tt_opt_grid(j+1)], xx_opt(:,j), options) ;
    
    plot(prop(:,1), prop(:,2), 'LineWidth', 2)
    hold on ; grid on; axis equal
    scatter(prop(1,1), prop(1,2),'filled', 'SizeData', 40, 'MarkerFaceColor','#8B0000') 
end
scatter(prop(end,1), prop(end,2),'filled', 'SizeData', 30, 'MarkerFaceColor','#8B0000')

% Scattering Planets
scatter(-mu, 0, 'filled', 'SizeData', 140, 'MarkerFaceAlpha',0.7 ,'MarkerFaceColor','#1E90FF') 
scatter(1-mu, 0, 'filled', 'SizeData', 130, 'MarkerFaceAlpha',0.8,'MarkerFaceColor','#D3D3D3') 

% Settings
set(gca,'FontSize',12)
xlabel('x [-]', 'FontSize', 20)
ylabel('y [-]', 'FontSize', 20)
legend('First segment', '', 'Second segment', '', 'Third segment','','Nodes', ...
       'Earth', 'Moon' ,'FontSize', 14)
title('Multiple Shooting Orbit (@Earth-Moon Rotating Frame)' , 'Interpreter', 'Latex', 'fontsize', 18)

%% FUNCTIONS

function xx_dot = PBRFBP_ode(tt, xx)
%
% xx_dot = PBRFBP_ode(tt, xx)
%
% Planar Bicircular Restricted Four-Body Problem : Provide the odefun for
%                                                  the dynamic integration
%
% INPUT:
%   tt         Time
%   xx         State:                                     [4x1]      
%                  - Position                             [2x1]          
%                  - Velocity                             [2x1]      
%
% OUTPUT:
%   xx_dot     State derivative                           [4x1]      
%
% AUTHOR:
%   Davide Lanza
%

% Physical parameters (from Topputo 2013)
mu = 1.21506683e-2 ;    % Earth-Moon mass parameter             [-] 
m_s = 3.28900541e5 ;    % Scaled mass of the Sun                [-]
rho = 3.88811143e2 ;    % Scaled Sun–(Earth + Moon) distance    [-]
w_s = -9.25195985e-1 ;  % Scaled angular velocity of the Sun    [-]

% State extraction [-]
x = xx(1) ;     % Position 
y = xx(2) ; 
v_x = xx(3) ;   % Velocity
v_y = xx(4) ; 

% Magnitude distances [-]
r1 = sqrt((x + mu)^2 + y^2) ;
r2 = sqrt((x + mu - 1)^2 + y^2) ; 
r3 = sqrt((x - rho*cos(w_s*tt))^2 + (y - rho*sin(w_s*tt))^2) ; 

% PBRFBP dynamics as a 1st order ODE: xx_dot = f(xx)
xx_dot = [ v_x ;
           v_y ;
           2*v_y + x - (1 - mu)*(x + mu)/r1^3 - mu*(x + mu - 1)/r2^3 - m_s*(x - rho*cos(w_s*tt))/r3^3 - m_s*cos(w_s*tt)/rho^2 ;
           -2*v_x + y - (1 - mu)*y/r1^3 - mu*y/r2^3 - m_s*(y - rho*sin(w_s*tt))/r3^3 - m_s*sin(w_s*tt)/rho^2                ] ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xxSTM_dot = PBRFBP_STM_ode(tt, xxSTM)
%
% xxSTM_dot = PBRFBP_STM_ode(tt, xxSTM, par))
%
% Planar Bicircular Restricted Four-Body Problem : Provide the odefun for
%                                                  the dynamic integration
%
% INPUT:
%   t          Time
%   xxSTM      State and State transition matrix        [20x1]      
%                - Position                             [2x1]          
%                - Velocity                             [2x1]  
%                - Reshaped STM                         [16x1]
%
% OUTPUT:
%   xx_dot     Derivatives of the state and the STM     [20x1]      
%
% AUTHOR:
%   Davide Lanza
%

% Physical parameters (from Topputo 2013)
mu = 1.21506683e-2 ;    % Earth-Moon mass parameter             [-] 
m_s = 3.28900541e5 ;    % Scaled mass of the Sun                [-]
rho = 3.88811143e2 ;    % Scaled Sun–(Earth + Moon) distance    [-]
w_s = -9.25195985e-1 ;  % Scaled angular velocity of the Sun    [-]

% State extraction [-]
x = xxSTM(1) ;     % Position 
y = xxSTM(2) ; 
v_x = xxSTM(3) ;   % Velocity
v_y = xxSTM(4) ; 

% Magnitude distances [-]
r1 = sqrt((x + mu)^2 + y^2) ;
r2 = sqrt((x + mu - 1)^2 + y^2) ; 
r3 = sqrt((x - rho*cos(w_s*tt))^2 + (y - rho*sin(w_s*tt))^2) ; 

% PBRFBP dynamics as a 1st order ODE: xx_dot = f(xx)
xx_dot = [ v_x ;
           v_y ;
           2*v_y + x - (1 - mu)*(x + mu)/r1^3 - mu*(x + mu - 1)/r2^3 - m_s*(x - rho*cos(w_s*tt))/r3^3 - m_s*cos(w_s*tt)/rho^2 ;
           -2*v_x + y - (1 - mu)*y/r1^3 - mu*y/r2^3 - m_s*(y - rho*sin(w_s*tt))/r3^3 - m_s*sin(w_s*tt)/rho^2                ] ;


% Variational approach: A = df/dxx (f = xx_dot) 
A = [                                                                        0                                                               0    1     0     
                                                                             0                                                               0    0     1     
     1 - (1-mu)/r1^3 + 3*(1-mu)*(x+mu)^2/r1^5 - mu/r2^3 + 3*mu*(x-1+mu)^2/r2^5                   3*(1-mu)*(x+mu)*y/r1^5 + 3*mu*(x+mu-1)*y/r2^5    0     2     
                                 3*(1-mu)*(x+mu)*y/r1^5 + 3*mu*(x+mu-1)*y/r2^5   1 - (1-mu)/r1^3 + 3*(1-mu)*y^2/r1^5 - mu/r2^3 + 3*mu*y^2/r2^5   -2     0 ];

% Reshaping into a matrix
STM = reshape(xxSTM(5:20), [4 4]) ;

% Computing STM derivative
STM_dot = A * STM;

% Building the output 
xxSTM_dot(1:4,1) = xx_dot ; % State derivatives
xxSTM_dot(5:20,1) = reshape(STM_dot, [16 1]) ; % Reshaped STM derivatives

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, GF] = myObjFnc(xx)
%
% [f, GF] = myObjFnc(xx)
% 
% SIMPLE SHOOTING: Objective function to minimize the total impulse of the mission 
%
% INPUT:   
%   xx         NLP variables                  [6x1]
%                - state vector               [4x1] 
%                - Initial and final time     [2x1]
%
% OUTPUT:   
%   f         Total impulse                   [1x1]
%   GF        Gradient of the obj function    [6x1]
%
% AUTHOR:
%   Davide Lanza
%

% Physical parameters (from Topputo 2013)
mu = 1.21506683e-2 ;    % Earth-Moon mass parameter             [-] 
Re = 6378e3 ;           % Mean Earth's radius                   [m]
Rm = 1738e3 ;           % Mean Moon's radius                    [m]
hi = 167e3 ;            % Altitude of departure orbit           [m]
hf = 100e3 ;            % Altitude of arrival orbit             [m]
DU = 3.84405000e8 ;     % Distance unit                         [m]

t_i = xx(5) ; 
t_f = xx(6) ; 

% Initial state + STM
STM0 = reshape(eye(4), [16,1]) ; % [16x1]

xxSTM0 = [xx(1:4); STM0] ; % [20x1]
tspan = [t_i t_f] ; 
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Dynamic integration 
[~, xxSTM] = ode113(@PBRFBP_STM_ode, tspan, xxSTM0, options) ;

% Initial state
x_i = xx(1) ;  % Position 
y_i = xx(2) ; 
v_xi = xx(3) ; % Velocity 
v_yi = xx(4) ; 

% Final state
x_f = xxSTM(end,1) ;
y_f = xxSTM(end,2) ;
v_xf = xxSTM(end,3) ;
v_yf = xxSTM(end,4) ;

% Physical Distances
r_i = (Re + hi)/DU ;    
r_f = (Rm + hf)/DU ; 

% Objective function (Total cost computation)
DV_i = sqrt((v_xi - y_i)^2 + (v_yi + x_i + mu)^2) - sqrt((1-mu)/r_i)  ;   % Initial impulse   
DV_f = sqrt((v_xf - y_f)^2 + (v_yf + x_f + mu - 1)^2) - sqrt((mu)/r_f) ;  % Final impulse 

% Cost to be minimized (total impulse)
f = DV_i + DV_f ;  

% Gradient computation
if nargout > 1
    % Extracting the state transition matrix
    STM = reshape(xxSTM(end, 5:20), [4 4]) ;

    % Computing the derivative of the state transition matrix at the initial and final times
    xxSTM_dot_i = PBRFBP_STM_ode(t_i, xxSTM0) ;
    xxSTM_dot_f = PBRFBP_STM_ode(t_f, [xxSTM(end,1:4)'; STM0]) ;

    % Extracting the partial derivatives of the state vector at the initial and final times
    f_xi = xxSTM_dot_i(1:4) ;
    f_xf = xxSTM_dot_f(1:4) ;

    % Computing the gradient of DeltaV_i
    z_i1 = v_xi - y_i ; 
    z_i2 = v_yi + x_i + mu ;

    z_f1 = v_xf - y_f ;
    z_f2 = v_yf + x_f + mu - 1 ;

    dDV_i_dx_i = (1/sqrt((z_i1)^2 + (z_i2)^2)) * [z_i2 ;-z_i1 ;z_i1 ;z_i2] ;
    dDV_i_dt_i = 0 ;
    dDV_i_dt_f = 0 ;

    % Computing the gradient of DeltaV_f
    dDV_f_dx_i = STM' * ((1/sqrt( z_f1^2 + z_f2^2)) * [z_f2; -z_f1; z_f1; z_f2]) ;
    dDV_f_dt_i = ((1/sqrt(z_f1^2 + z_f2^2)) * [z_f2 ;-z_f1 ; z_f1; z_f2])' * (-STM * f_xi) ;
    dDV_f_dt_f = ((1/sqrt(z_f1^2 + z_f2^2)) * [z_f2; -z_f1; z_f1; z_f2])' * f_xf ;

    % Collect the partial derivatives into a matrix
    GF = [dDV_i_dx_i + dDV_f_dx_i ; 
          dDV_i_dt_i + dDV_f_dt_i ; 
          dDV_i_dt_f + dDV_f_dt_f];  

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c, Ceq, J, Jeq] = myNlConstr(xx)
%
% [c, Ceq, J, Jeq] = myNlConstr(xx)
%
% SIMPLE SHOOTING: non linear constraints
%
% INPUT:   
%   xx         NLP variables                    [6x1]
%                - state vector                 [4x1] 
%                - Initial and final time       [2x1]
%
% OUTPUT:   
%   C         Non linear inequality constraints
%   Ceq       Non linear equality constraints
%   J         Gradient of c
%   Jeq       Gradient of Ceq
%
% AUTHOR:
%   Davide Lanza
%

% Physical parameters (from Topputo 2013)
mu = 1.21506683e-2 ;    % Earth-Moon mass parameter             [-] 
Re = 6378e3 ;           % Mean Earth's radius                   [m]
Rm = 1738e3 ;           % Mean Moon's radius                    [m]
hi = 167e3 ;            % Altitude of departure orbit           [m]
hf = 100e3 ;            % Altitude of arrival orbit             [m]
DU = 3.84405000e8 ;     % Distance unit                         [m]

t_i = xx(5) ; 
t_f = xx(6) ; 
% Initial state + STM
STM0 = reshape(eye(4), [16,1]) ; % [16x1]

xxSTM0 = [xx(1:4); STM0] ; % [20x1]
tspan = [t_i t_f] ; 
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Dynamic integration 
[~, xxSTM] = ode113(@PBRFBP_STM_ode, tspan, xxSTM0, options) ;

% Initial state
x_i = xx(1) ;  % Position 
y_i = xx(2) ; 
v_xi = xx(3) ; % Velocity 
v_yi = xx(4) ; 

% Final state
x_f = xxSTM(end,1) ;
y_f = xxSTM(end,2) ;
v_xf = xxSTM(end,3) ;
v_yf = xxSTM(end,4) ;

% Physical Distances
r_i = (Re + hi)/DU ;    
r_f = (Rm + hf)/DU ; 

% Equality constraints (Ceq using fmincon notation)
Ceq = [(x_i + mu)^2 + y_i^2 - r_i^2 ;                              % PSI(x_i) 
       (x_i + mu)*(v_xi - y_i) + y_i*(v_yi + x_i + mu) ;
       (x_f + mu - 1)^2 + y_f^2 - r_f^2 ;                          % PSI(x_f) 
       (x_f + mu - 1)*(v_xf - y_f) + y_f*(v_yf + x_f + mu - 1) ] ; 

% c = [t_i - t_f] ; 
c = [] ; 

% Gradient implementation
if nargout > 2
    % Extracting the state transition matrix from the final state vector
    STM = reshape(xxSTM(end, 5:20), [4 4]) ;
    
    % Computing the derivative of the state transition matrix at the initial and final times
    xxSTM_dot_i = PBRFBP_STM_ode(t_i, xxSTM0) ; % Derivative at initial time
    xxSTM_dot_f = PBRFBP_STM_ode(t_f, [xxSTM(end,1:4)'; STM0]) ; % Derivative at final time
    
    % Extracting the partial derivatives of the state vector at the initial and final times from the state transition matrix derivatives
    f_xi = xxSTM_dot_i(1:4) ; % Derivative of state vector at initial time
    f_xf = xxSTM_dot_f(1:4) ; % Derivative of state vector at final time
        
    % Partial derivatives of constraint function wrt xx_i (state)
    dc1dx1 = [(2*(x_i + mu)) ; (2*y_i) ; 0 ; 0 ]; 
    dc2dx1 = [v_xi ; v_yi ; (x_i + mu) ; y_i ];
    dc3dx1 = [(2*(x_f+mu-1)) ; 2*y_f ; 0 ; 0 ]' * STM; % [1x4]
    dc4dx1 = [v_xf ; v_yf ; (x_f+mu-1) ; y_f ]' * STM;
    
    % Partial derivatives ofconstraint function wrt t_i 
    dc1dt1 = 0;
    dc2dt1 = 0;
    dc3dt1 = [(2*(x_f + mu - 1)) ;2*y_f ; 0 ; 0 ]' * (- STM * f_xi);
    dc4dt1 = [v_xf ; v_yf ;(x_f + mu - 1) ; y_f ]' * (- STM * f_xi);

    % Partial derivative of constraint function wrt t_f
    dc1dt2 = 0;
    dc2dt2 = 0;
    dc3dt2 = [(2*(x_f + mu - 1)) ; 2*y_f ; 0 ; 0]' * f_xf;
    dc4dt2 = [v_xf ; v_yf ; (x_f + mu - 1) ; y_f]' * f_xf;
    
    % Combining the partial derivatives into a matrix
    Jeq = [ dc1dx1 dc2dx1 dc3dx1' dc4dx1';
             dc1dt1 dc2dt1 dc3dt1 dc4dt1;
             dc1dt2 dc2dt2 dc3dt2 dc4dt2 ];

    % J = [0; 0; 0; 0; 1; -1] ; % Deleting it since it just slow down the
    % code but with the same results
    J = [] ; 
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, GF] = myObjFnc_MS(yy)
%
% [f, GF] = myObjFnc_MS(yy)
% 
% MULTIPLE SHOOTING: Objective function to minimize the total impulse of the mission 
%
% INPUT:   
%   yy         NLP variables                  [18x1]
%                - state vector               [16x1] 
%                - Initial and final time     [2x1]
%
% OUTPUT:   
%   f         Total impulse                   [1x1]
%   GF        Gradient of the obj function    [18x1]
%
% AUTHOR:
%   Davide Lanza
%

% Physical parameters (from Topputo 2013)
mu = 1.21506683e-2 ;    % Earth-Moon mass parameter             [-] 
Re = 6378e3 ;           % Mean Earth's radius                   [m]
Rm = 1738e3 ;           % Mean Moon's radius                    [m]
hi = 167e3 ;            % Altitude of departure orbit           [m]
hf = 100e3 ;            % Altitude of arrival orbit             [m]
DU = 3.84405000e8 ;     % Distance unit                         [m]

% Unpack the NLP variable:
% Initial state
x_i = yy(1) ;  % Position 
y_i = yy(2) ; 
v_xi = yy(3) ; % Velocity 
v_yi = yy(4) ; 

% Final state
x_f = yy(13) ;  % Position 
y_f = yy(14) ;  
v_xf = yy(15) ; % Velocity 
v_yf = yy(16) ; 

% Physical Distances
r_i = (Re + hi)/DU ;    
r_f = (Rm + hf)/DU ; 

% Substitutions
z_i1 = v_xi - y_i ;
z_i2 = v_yi + x_i + mu ;

z_f1 = v_xf - y_f ;
z_f2 = v_yf + x_f + mu - 1 ;

% Objective function (Total cost computation)
DV_i = sqrt((z_i1)^2 + (z_i2)^2) - sqrt((1-mu)/r_i)  ;   % Initial impulse   
DV_f = sqrt((z_f1)^2 + (z_f2)^2) - sqrt((mu)/r_f) ;  % Final impulse 

% Cost to be minimized (total impulse)
f = DV_i + DV_f ;  

% Gradient Computation
df_dx_i = (1/sqrt((z_i1)^2 + (z_i2)^2)) * [z_i2 ;-z_i1 ;z_i1 ;z_i2] ;
    
df_dx_j = zeros(8,1) ; % Null terms (j = 2, 3)
    
df_dx_f = (1/sqrt((z_f1)^2 + (z_f2)^2)) * [z_f2 ;-z_f1 ;z_f1 ;z_f2] ;
    
df_dt_i = 0 ; 
df_dt_f = 0 ; 
    
% Assembling into the output vecotor
GF = [df_dx_i; df_dx_j; df_dx_f; df_dt_i; df_dt_f] ; % [18x1]

    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C, Ceq, J, Jeq] = myNlConstr_MS(yy)
%
% [C, Ceq, J, Jeq] = myNlConstr_MS(yy)
%
% MULTIPLE SHOOTING: non linear constraints
%
% INPUT:   
%   yy         NLP variables                    [18x1]
%                - state vector                 [16x1] 
%                - Initial and final time       [2x1]
%
% OUTPUT:   
%   C         Non linear inequality constraints
%   Ceq       Non linear equality constraints
%   J         Gradient of c
%   Jeq       Gradient of Ceq
%
% AUTHOR:
%   Davide Lanza
%

% Physical parameters (from Topputo 2013)
mu = 1.21506683e-2 ;    % Earth-Moon mass parameter             [-] 
Re = 6378e3 ;           % Mean Earth's radius                   [m]
Rm = 1738e3 ;           % Mean Moon's radius                    [m]
hi = 167e3 ;            % Altitude of departure orbit           [m]
hf = 100e3 ;            % Altitude of arrival orbit             [m]
DU = 3.84405000e8 ;     % Distance unit                         [m]


% Unpacking the NLP vector 
xx = reshape(yy(1:16), [4,4]) ; % Matrix containing state vectors @ each node
t_i = yy(end-1) ; 
t_f = yy(end) ; 

% Time grid definition
N = 4 ;                      % Number of points

tt_grid = zeros(1, N) ; 
for j = 1 : N
    tt_grid(j) = t_i + ((j - 1)/(N - 1)) * (t_f - t_i) ; % Time grid of N points
end

% INEQUALITY CONSTRAINTS: to avoid solutions which impact either the Earth or the Moon.
C = zeros(9, 1);

C(1:8) = [ (Re/DU)^2 - (xx(1,1) + mu)^2 - xx(2,1)^2 ;
      (Rm/DU)^2 - (xx(1,1) + mu - 1)^2 - xx(2,1)^2 ;
      (Re/DU)^2 - (xx(1,2) + mu)^2 - xx(2,2)^2 ;
      (Rm/DU)^2 - (xx(1,2) + mu - 1)^2 - xx(2,2)^2 ;
      (Re/DU)^2 - (xx(1,3) + mu)^2 - xx(2,3)^2 ;
      (Rm/DU)^2 - (xx(1,3) + mu - 1)^2 - xx(2,3)^2 ;
      (Re/DU)^2 - (xx(1,4) + mu)^2 - xx(2,4)^2 ;
      (Rm/DU)^2 - (xx(1,4) + mu - 1)^2 - xx(2,4)^2] ;

C(end) = t_i - t_f ; %  transfer duration has to be strictly positive

% EQUALITY CONSTRAINTS: - defects constraints
%                       - continuity of the position @ t_i and t_f 

% Initial matrix containing states and STM @ node 1, 2, 3 [20x3]
xxSTM_0 = [               xx(:,1),                xx(:,2),              xx(:,3)   ;
           reshape(eye(4), [16,1]), reshape(eye(4), [16,1]), reshape(eye(4), [16,1])] ; 

phi = zeros(4, N-1) ; % Flux matrix allocation 
STM = zeros(16, N-1) ; % State transition matrix allocation

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

for i = 1 : (N - 1)

    % Dynamic integration
    [~, xxSTM] = ode113(@PBRFBP_STM_ode, [tt_grid(i) tt_grid(i+1)], xxSTM_0(:,i), options) ; 
    
    phi(:,i) = xxSTM(end, 1:4) ; % Flux @ final time of the segment
    STM(:,i) = xxSTM(end, 5:20) ; % State transition matrix @ final time of the segment
    
end

Defect =[phi(:,1) - xx(:,2) ;
        phi(:,2) - xx(:,3) ;
        phi(:,3) - xx(:,4)] ;

% Physical Distances
r_i = (Re + hi)/DU ;    
r_f = (Rm + hf)/DU ; 

% Renaming for semplicity 
x_i = xx(1,1) ;         x_f = xx(1,4) ; 
y_i = xx(2,1) ;         y_f = xx(2,4) ;
v_xi = xx(3,1) ;        v_xf = xx(3,4) ;
v_yi = xx(4,1) ;        v_yf = xx(4,4) ; 

% continuity of the position @ t_i and t_f and velocity condition
cont = [(x_i + mu)^2 + y_i^2 - r_i^2 ;
        (x_i + mu)*(v_xi - y_i) + y_i*(v_yi + x_i + mu) ;
        (x_f + mu - 1)^2 + y_f^2 - r_f^2 ; 
        (x_f + mu - 1)*(v_xf - y_f) + y_f*(v_yf + x_f + mu - 1)] ;

% Assemblying 
Ceq = [Defect ;
       cont] ; 

     % Gradient implementation for inequality constraints
J = zeros(9,18);
J(1:2,1:4) = [ -2*(xx(1,1)+mu)   -2*xx(2,1)  0   0 ;
             -2*(xx(1,1)+mu-1)   -2*xx(2,1)  0   0] ;
    
J(3:4,5:8) = [ -2*(xx(1,2)+mu)   -2*xx(2,2)   0   0;
             -2*(xx(1,2)+mu-1)   -2*xx(2,2)   0   0];
    
J(5:6,9:12) = [ -2*(xx(1,3)+mu)   -2*xx(2,3)   0   0;
              -2*(xx(1,3)+mu-1)   -2*xx(2,3)   0   0];
    
J(7:8,13:16) = [-2*(xx(1,4)+mu)  -2*xx(2,4)  0  0;
               -2*(xx(1,4)+mu -1)  -2*xx(2,4)  0  0];
    
J(9,17:18) = [1 -1]; % time elements 
    
% Transposing the matrix
J = J' ; 
    
% Gradient implementation for equality constraints  
phi = zeros(4, N-1) ; % Flux matrix allocation 
fun_phi = zeros(4,3); 
fun = zeros(4,3) ; 
STM = zeros(16, N-1) ; % State transition matrix allocation
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    
for i = 1 : (N - 1)
        % Dynamic integration
        [~, xxSTM] = ode113(@PBRFBP_STM_ode, [tt_grid(i), tt_grid(i+1)], xxSTM_0(:,i), options); 
        
        phi(:,i) = xxSTM(end, 1:4) ; % Flux @ final time of the segment
        STM(:,i) = xxSTM(end, 5:20) ; % State transition matrix @ final time of the segment
        xxSTM_dot = PBRFBP_STM_ode(tt_grid(i), xxSTM_0(:,i)) ; 
        fun(:,i) = xxSTM_dot(1:4) ;
        xxSTM_dot_phi = PBRFBP_STM_ode(tt_grid(i+1), [phi(:,i); reshape(eye(4), [16,1])]) ; 
        fun_phi(:,i) = xxSTM_dot_phi(1:4) ;   
end
    
% Allocation matrix 
Jeq = zeros(4*N, ((4*N) + 2));
    
% Populating the matrix with values from STM matrix
Jeq(1:4,1:4) = reshape(STM(:,1), [4 4]);
Jeq(5:8,5:8) = reshape(STM(:,2), [4 4]);
Jeq(9:12,9:12) = reshape(STM(:,3), [4 4]);
    
% Populating the matrix with values from identity matrix
Jeq(1:4,5:8) = -eye(4);
Jeq(5:8,9:12) = -eye(4);
Jeq(9:12,13:16) = -eye(4);
    
% Populating the matrix with values calculated using fun and fun_phi vectors
Jeq(1:4,17) = reshape( STM(:,1), [4 4] ) * fun(:,1)* (-(N-1)/(N-1)) + fun_phi(:,1) * (N-2)/(N-1);
Jeq(1:4,18) = fun_phi(:,1) * (1/(N-1));
    
Jeq(5:8,17) = reshape( STM(:,2), [4 4] ) * fun(:,2)*(-(N-2)/(N-1)) + fun_phi(:,2) * (N-3)/(N-1);
Jeq(5:8,18) = reshape( STM(:,2), [4 4] ) * fun(:,2) * -(1/(N-1)) + fun_phi(:,2) * (2/(N-1));
    
Jeq(9:12,17) = reshape( STM(:,3), [4 4] ) * fun(:,3)* (-(N-3)/(N-1)) + fun_phi(:,3) * (N-4)/(N-1);
Jeq(9:12,18) = reshape( STM(:,3), [4 4] ) * fun(:,3) * -(2/(N-1)) + fun_phi(:,3) * (3/(N-1));
    
% Populating the matrix with values for the initial and final positions and velocities
Jeq(13:14,1:4) = [2*(x_i + mu) , 2*y_i , 0 , 0;
                       v_xi , v_yi , x_i + mu , y_i ];
Jeq(15:16,13:16) = [2*(x_f + mu - 1) , 2*y_f , 0, 0;
                         v_xf , v_yf , (x_f + mu - 1) , y_f];
    
% Transposing the matrix
Jeq = Jeq';

end

