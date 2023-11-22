% Spacecraft Guidance and Navigation (2022/2023)
% Assignment # 1
% Exercise # 1
% Author: Davide Lanza

%% EXERCISE 1.1

% Clear workspace, close all figures, and clear command window
clearvars; close all; clc;

% Define physical parameter of the Sun-Earth system
mu = 3.0359e-6 ;

% Define function to compute the second libration point (adimensional rotating referance frame)

dUdx = @(x) x - (1 - mu)*(x + mu)./abs(x + mu).^ 3 - mu*(x + mu - 1)./abs(x + mu - 1).^ 3 ;

% Set options for fsolve function to stop when the root is found with a tolerance of 10^-10
options = optimset('TolX', 1e-10) ;

% Compute the second libration point
L2 = fzero(dUdx, [0.9, 1.1], options) ;

% Display the result with 10-digit accuracy
fprintf('x-coordinate of the Lagrange point L2: %.10f\n', L2) ;

% Defining plotting intervals
x1 = linspace(-2, 0-0.001, 10000) ;
x2 = linspace(0+0.00001, 1-0.00001, 10000) ; 
x3 = linspace(1+0.00001, 2, 10000) ;

% Plot function and the second libration point
figure(1)
plot(x1, dUdx(x1), 'k',  x2, dUdx(x2), 'k', x3, dUdx(x3), 'k', 'LineWidth', 1.5) ;
xlim([-1.8 1.8]) ; % Plot settings
ylim([-50 50]) ;
hold on; grid on;
scatter(L2, 0, 'filled', 'SizeData', 200, 'MarkerFaceColor','#8B0000') 
scatter(-mu, 0, 'filled', 'SizeData', 200, 'MarkerFaceColor','#FFE600') 
scatter(1-mu, 0, 'filled', 'SizeData', 200, 'MarkerFaceColor','#0080FF') 

% Plot settings
set(gca,'FontSize',12)
xlabel('$x \space [-]$','Interpreter','latex','FontSize', 25)
ylabel('$\frac{\partial U}{\partial x} \space [-] $','Interpreter','latex','FontSize', 25)
legend('','','', '$L_2$', 'Sun', 'Earth','Interpreter','latex','FontSize', 20)
title('Second Libration point (@Sun-Earth Barycenter rotating frame)' , 'Interpreter', 'Latex', 'fontsize', 18)

% Zoom in around (1,0)
figure(2)
plot(x1, dUdx(x1), 'k',  x2, dUdx(x2), 'k', x3, dUdx(x3), 'k', 'LineWidth', 1.5) ;
xlim([0.96 1.04]) ; % Plot settings
ylim([-1.5 1.5]) ;
hold on; grid on; 
scatter(L2, 0, 'filled', 'SizeData', 200, 'MarkerFaceColor','#8B0000') 
scatter(-mu, 0, 'filled', 'SizeData', 200, 'MarkerFaceColor','#FFE600') 
scatter(1-mu, 0, 'filled', 'SizeData', 200, 'MarkerFaceColor','#0080FF') 

% Plot settings
set(gca,'FontSize',12)
xlabel('$x \space [-]$','Interpreter','latex','FontSize', 25)
ylabel('$\frac{\partial U}{\partial x} \space [-]$','Interpreter','latex','FontSize', 25)
legend('','','', '$L_2$', '', 'Earth', 'Interpreter','latex','FontSize', 20)
title('zoom in around L2' , 'Interpreter', 'Latex', 'fontsize', 18)

%% EXERCISE 1.2

% Clear workspace, close all figures, and clear command window
clearvars; close all; clc;

% Define physical parameter of the Sun-Earth system
mu = 3.0359e-6 ;

% Inital guess state
xx0 = [1.008296144180133; 0; 0.001214294450297; 0; 0.010020975499502; 0] ; 

% Second libration point (from previous point)
L2 = [1.0100701875 0 0] ; 

% Differential correction scheme
Tol = 1e-8 ; % Tolerance
[N_it, xx0STM0] = diff_corr(xx0, mu, Tol) ;  

% CRTBP integration 
options = odeset('RelTol', 1e-13,'AbsTol', 1e-13,'Events', @myEvent) ; % Setting the event options
[tt, xxSTM, tt_e, xxSTM_e] = ode113( @(tt, xxSTM) CRTBP_ode(xxSTM, mu), [0 2*pi], xx0STM0, options) ; % ODE integration

% Plotting the solution 
figure(1)
plot3(xxSTM(:,1), xxSTM(:,2), xxSTM(:,3),'k','LineWidth',2) ; 
hold on; 
plot3(xxSTM(:,1), -xxSTM(:,2), xxSTM(:,3),'--k','LineWidth',2) ; 


scatter3(xx0(1), xx0(2), xx0(3), 'filled', 'SizeData', 100, 'MarkerFaceColor','#1a237e' ) 
scatter3(L2(1), L2(2), L2(3), 'filled', 'SizeData', 100, 'MarkerFaceColor','#8B0000' ) 

% Plot settings 
grid on; axis equal;
set(gca,'FontSize',12)
xlabel('x [-]','FontSize', 20) ; ylabel('y [-]','FontSize', 20) ; zlabel('z [-]','FontSize', 20) ;
lgd = legend('Half-halo orbit', 'Symmetric halo orbit', 'Initial condition', ...
            '$L_2$','Interpreter','latex','FontSize', 16) ; 
lgd.Location = 'northwest' ; 
title('Periodic Halo Orbit (@Sun-Earth Barycenter rotating frame)' , 'Interpreter', 'Latex', 'fontsize', 18)



%% EXERCISE 1.3

% Clear workspace, close all figures, and clear command window
clearvars; close all; clc;

% Define physical parameter of the Sun-Earth system
mu = 3.0359e-6 ;

% Inital guess state
xx0 = [1.008296144180133; 0; 0.001214294450297; 0; 0.010020975499502; 0] ; 

% Second libration point (from previous point)
L2 = [1.0100701875 0 0] ; 


% Tolerance
Tol = 1e-8 ; 

z0_f = 0.0046; % Final initial z0 guess
zz = linspace(xx0(3), z0_f, 20) ; % Vector containing all the initial z0 

% Retrieving and plotting the family of halo orbits
figure(1) 
scatter3(L2(1), L2(2), L2(3), 'filled', 'SizeData', 100, 'MarkerFaceColor','#8B0000' )
hold on ; grid on; axis equal;

for i = 1:length(zz)

    xx0(3) = zz(i) ; % Gradually increasing z0  
    [N_it, xx0STM0] = diff_corr( xx0, mu, Tol); % Differential correction scheme
    options = odeset('RelTol', 1e-13,'AbsTol', 1e-13,'Events', @myEvent) ; % Setting the event options
    [tt, xxSTM, tt_e, xxSTM_e] = ode113( @(tt, xxSTM) CRTBP_ode(xxSTM, mu), [0 2*pi], xx0STM0, options) ; % ODE integration

    plot3(xxSTM(:,1), xxSTM(:,2), xxSTM(:,3), 'k','LineWidth', 1.5) ; 
    plot3(xxSTM(:,1), -xxSTM(:,2), xxSTM(:,3), '--k','LineWidth', 1.5) ; 

    
    xx0 = xx0STM0(1:6); % Updating the initial conditions
end
 
% Plot settings
grid on; axis equal;
set(gca,'FontSize',12)
xlabel('x [-]','FontSize', 17) ; ylabel('y [-]','FontSize', 17) ; zlabel('z [-]','FontSize', 17) ;
lgd = legend('$L_2$', 'Family of halo orbits', 'Interpreter', 'latex','FontSize', 18) ; 
lgd.Location = 'northeast' ; 

%% FUNCTIONS

function xxSTM_dot = CRTBP_ode(xxSTM, mu)
%
% dxxSTMdt = CRTBP_ode(xxSTM, mu)
%
% 3D Circular Restricted Three-Body Problem : Provide the odefun for integration 
%                                             of the 3D CRTBP
%
% INPUT:
%   xxSTM      Vector containing:                  [42x1]      
%                  - State of the system           [6x1]       [-]
%                  - State transition matrix       [36x1]      [-]
%   mu         Sun-Earth physical parameter        [1x1]       [-]           
%
% OUTPUT:
%   xxSTM_dot  Derivatives                         [42x1]      [-]
%
% AUTHOR:
%   Davide Lanza
%

% Position extraction  
x = xxSTM(1) ;
y = xxSTM(2) ;
z = xxSTM(3) ;

% Velocity extraction
vx = xxSTM(4) ;
vy = xxSTM(5) ;
vz = xxSTM(6) ;

r1 = sqrt((x+mu)^2 + y^2 + z^2) ; % Sun-S/c magnitude distance
r2 = sqrt((x+mu-1)^2 + y^2 + z^2) ; % Earth-S/c magnitude distance 

% CRTBP dynamics as a 1st order ODE: xx_dot = f(xx)
xx_dot = [                                              vx
                                                        vy
                                                        vz
          2*vy + x - (1-mu)*(x+mu)/r1^3 - mu*(x+mu-1)/r2^3
                     -2*vx + y - (1-mu)*y/r1^3 - mu*y/r2^3 
                                -(1-mu)*z/r1^3 - mu*z/r2^3 ] ;

% Variational approach: A = df/dxx (f = xx_dot) 
A = [                                                                        0                                                                 0                                                           0     1     0     0
                                                                             0                                                                 0                                                           0     0     1     0
                                                                             0                                                                 0                                                           0     0     0     1
     1 - (1-mu)/r1^3 + 3*(1-mu)*(x+mu)^2/r1^5 - mu/r2^3 + 3*mu*(x-1+mu)^2/r2^5                     3*(1-mu)*(x+mu)*y/r1^5 + 3*mu*(x+mu-1)*y/r2^5               3*mu*z*(mu+x-1)/r2^5 + 3*(1-mu)*z*(mu+x)/r1^5     0     2     0
                                 3*(1-mu)*(x+mu)*y/r1^5 + 3*mu*(x+mu-1)*y/r2^5     1 - (1-mu)/r1^3 + 3*(1-mu)*y^2/r1^5 - mu/r2^3 + 3*mu*y^2/r2^5                           3*(1-mu)*y*z/r1^5 + 3*mu*y*z/r2^5    -2     0     0 
                                 3*mu*z*(mu+x-1)/r2^5 + 3*(1-mu)*z*(mu+x)/r1^5                                 3*(1-mu)*y*z/r1^5 + 3*mu*y*z/r2^5   3*mu*z^2/r2^5 + 3*(1-mu)*z^2/r1^5 - mu/r2^3 - (1-mu)/r1^3     0     0     0 ];

% Reshaping into a matrix
STM = reshape(xxSTM(7:42), [6 6]) ; 

% Computing STM derivative
STM_dot = A * STM ; 

% Building the output 
xxSTM_dot(1:6,1) = xx_dot ; % State derivatives 
xxSTM_dot(7:42,1) = reshape(STM_dot, [36 1]) ; % Reshaped STM derivatives 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [position, isterminal, direction] = myEvent(~, xxSTM)
%
% [position, isterminal, direction] = myEvent(t, xxSTM, mu)
%
% 'Events' option of the odeset function to specify an event function.
%
% INPUT:
%   t           time                              [1x1]       [s]
%   xxSTM     Vector containing:                  [42x1]      
%                  - State of the system          [6x1]       [-]
%                  - State transition matrix      [36x1]      [-]
%   mu        Sun-Earth physical parameter        [1x1]       [-]   
%
% OUTPUT:
%   position       
%   isterminal
%   direction 
%
% AUTHOR:
% Davide Lanza
%

position = xxSTM(2) ; % We want second component to to zero 
isterminal = 1 ; % Halt integration 
direction = -1 ; % Locates only zeros where the event function is decreasing
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N_it, xxSTM_0] = diff_corr(xx_0, mu, Tol)
%
% xx0STM0 = diff_corr(xx0, mu, Tol)
%
% Differential correction scheme : Find the correction to Delta xx  
%
% INPUT:
%   xx0       State vector                        [6x1]       [-]
%   mu        Sun-Earth physical parameter        [1x1]       [-]    
%   Tol       Tolerance                           [1x1]
%
% OUTPUT:
%   xx0STM0   Vector containing the corrected     [42x1]
%             state vector + STM 
%   N_it      Number of iteration                 [1x1]       [-]
%
% AUTHOR:
%   Davide Lanza
%

STM_0 = reshape(eye(6), [36 1]) ; % Initializing and reshaping STM into a vector 
xxSTM_0 = [ xx_0; STM_0 ] ; % Initial conditions
err = Tol + 1 ; % Entering the cycle 
N_it = 0 ; 

while err > Tol
    % Propagating xx and STM until reaching the event condition  
    options = odeset('RelTol', 1e-13,'AbsTol', 1e-13,'Events', @myEvent) ;
    [~, ~, ~, xxSTM_e] = ode113(@(tt, xxSTM) CRTBP_ode(xxSTM, mu), [0 2*pi], xxSTM_0, options) ; 
    
    xx_f = xxSTM_e(1:6) ; % Final state vector 
    STM_f = reshape( xxSTM_e(7:42), [6 6] ) ; % Reashping the final STM into a matrix 

    % Correction scheme
     A = [STM_f(4,1) STM_f(4,5) 
         STM_f(6,1) STM_f(6,5)] ; % Just the relevant components 

    c = -[xx_f(4); xx_f(6)] ; % vx, vz

    Deltaxx = A\c ; % Correction vector 

    % Updating the initial conditions
    xx_0(1) = xx_0(1) + Deltaxx(1) ;
    xx_0(5) = xx_0(5) + Deltaxx(2) ;
    xxSTM_0 = [xx_0; STM_0] ;

    % Error computation 
    err = max(abs(Deltaxx(1)/xx_f(1)), abs(Deltaxx(2)/xx_f(5))) ;
    
    % Iteration counter
    N_it = N_it + 1 ; 

end

end