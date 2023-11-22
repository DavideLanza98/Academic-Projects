% Modeling and Simulation of Aerospace Systems (2022/2023)
% Assignment # 1
% Author: Davide Lanza 

%% EXERCISE 1
clearvars; close all; clc; 

% Defining function 
f = @(x) sin(x) + x - 1;
x = linspace(-10,10,100);

% Plot of the function 
plot(x,f(x),'b','LineWidth',2); % Plotting the function
grid on; hold on;
plot(x, zeros(1,length(x)),'-- k','LineWidth',2) ; 
set(gca,'FontSize',12)

% a and b guesses such that f(a)f(b) < 0 from the function plot
a = -5 ; b = 5;
plot(a,f(a),'o r','LineWidth',7) ; % Plot of the initial condition
plot(b,f(b),'o r','LineWidth',7) ; % Plot of the final condition 
plot(0.510973, 0, '* c','LineWidth',7) % Solution we are looking for 

% Plot settings 
lgd = legend('$f(x) = sin(x) + x - 1$','$y = 0$','$f(a), f(b)$','','$f(x^*)$','Interpreter','latex');
xlab = xlabel('$x$', 'Interpreter', 'latex');
ylab = ylabel('$f(x)$', 'Interpreter', 'latex');
% tit = title('Function representation');
lgd.FontSize = 15; lgd.Location = 'southeast';
xlab.FontSize = 20; ylab.FontSize = 20; 
% tit.FontSize = 20;
ha = gca;
ha.XAxisLocation = 'origin'; % Centering the reference axis at 0 
ha.YAxisLocation = 'origin';
ha.YLim = [-7 7]; ha.XLim = [-8 8];

%%% 1.1) Bisection method (use of Bisection.m)
Tol = 10^-10; n = 8 ; % Tollerance and accuracy requested 

tic
[~, N_it1, acc_bis] = Bisection(f, a, b, Tol, n) ;
Time_bis = toc;

figure() ; 
plot(acc_bis,'o','LineWidth',5) ; 
set(gca, 'YScale', 'log') ; hold on; grid on; 
set(gca,'FontSize',12)

%%% 1.2) Secant method (use of Secant.m)
nmax = 100;

tic
[~, N_it2, acc_sec] = Secant(f, -5, 5, nmax, Tol, n);
Time_sec = toc ;
plot(acc_sec,'o','LineWidth',5) ;

%%% 1.3) Regula Falsi method (use of Regula_Falsi.m)
x0 = a; x1 = b; 
tic
[~, N_it3, acc_falsi] = Regula_Falsi(f, x0, x1, nmax, Tol, n);
Time_falsi = toc ; 
plot(acc_falsi,'o','LineWidth',5) ;

plot(linspace(0,30,100), ones(1,100)*10^-8,'--k','LineWidth',2)

legend('Bisection','Seacant',' Regula Falsi','$Accuracy = 10^-8$','Interpreter','latex','FontSize',15)

xlabel('$N_{it}$', 'Interpreter', 'Latex', 'fontsize', 20);
ylabel('$|x_k - x_{k-1}|$', 'Interpreter', 'Latex', 'fontsize', 20);
% title('Accuracy convergence','FontSize', 20)

% Print of computational times 
fprintf('COMPUTATIONAL TIME [s]\n')
fprintf('   Bisection : %-12.8f\n', Time_bis); % Printing the final root
fprintf('   Secant : %-12.8f\n', Time_sec); % Printing the final root
fprintf('   Regula Falsi : %-12.8f\n', Time_falsi); % Printing the final root
fprintf('-------------------------\n')

%% EXERCISE 2
clearvars; close all; clc; 

% Function definition 
f = @(x) [x(1).^2 + x(2) - 5 ; x(2).^2 - x(1)] ; 
f1 = @(x) x(1).^2 + x(2) - 5 ;
f2 = @(x) x(1).^2 - x(2) ;

% Plot of the function 
z1 = @(x1,x2) x1.^2 + x2 - 5 ;
zhandle = fcontour(z1,'b','LineWidth',1.5) ; 
zhandle.LevelList = 0 ; 
grid on ;  hold on ; 

z2 = @(x1,x2) x2.^2 - x1 ;
zhandle2 = fcontour(z2,'k','LineWidth',1.5) ; 
zhandle2.LevelList = 0 ; 

% Initial guesses
xx0_1 = [2; 2] ;  
xx0_2 = [2; -2] ;
plot(xx0_1(1),xx0_1(2),'k s','LineWidth',7) ; plot(xx0_2(1),xx0_2(2),'k s','LineWidth',7) ;

%%% Exercise 2.1)
% Analytic computation of partial derivatives 

J = @(x) [2*x(1) 1; -1 2*x(2)] ; % Jacobian: matrix of partial derivatives
  
Tol = 10^-10; kmax = 100 ;

% First solution 
[xx1,N_iter1,err1] = newton(f, J, xx0_1, Tol, kmax) ; 
plot(xx1(1,end), xx1(2,end),'r o','LineWidth', 7) ;
% Second solution
[xx2,N_iter2, err2] = newton(f, J, xx0_2, Tol, kmax) ; 
plot(xx2(1,end), xx2(2,end),'r o','LineWidth', 7) ;

% Plot
legend('$f_1(x_1,x_2) = 0$','$f_2(x_1,x_2) = 0$','','Initial guess','Analytic derivatives', 'Interpreter', 'latex','FontSize', 12) ;
set(gca,'FontSize',12)
xlabel('$x_1$','Interpreter', 'latex','FontSize',20)
ylabel('$x_2$','Interpreter', 'latex','FontSize',20)

%%% Exercise 2.2)
% Estimation of derivatives through finite differences both forward and centered 

type1 = 'forward' ; 
type2 = 'centred' ; 

% First solution
[xx3, N_iter3, err3] = newton_diff_fin(f, xx0_1, Tol, kmax, type1); % forward
[xx4, N_iter4, err4] = newton_diff_fin(f, xx0_1, Tol, kmax, type2); % Centred


% Second solution 
[xx5, N_iter5, err5] = newton_diff_fin(f, xx0_2, Tol, kmax, type1); % Forward
[xx6, N_iter6, err6] = newton_diff_fin(f, xx0_2, Tol, kmax, type2); % Centred

% Error plots 
figure; 
plot(1:N_iter1, abs(err1)','b','LineWidth',2)
set(gca, 'YScale', 'log') ; hold on; grid on; 

plot(1:N_iter1, abs(err3)','r o','MarkerSize',15,'LineWidth',2)

%plot(1:N_iter1, abs(err4)','g o','MarkerSize',15,'LineWidth',2) % Not
%needed in the discussion 

legend('Analytic','Forward','FontSize', 15) ;
xlabel('No. of Iterations','FontSize', 15)
ylabel('Error','FontSize',15)

% Print of results 
fprintf('FIRST SOLUTION\n')
fprintf('   Analytic : [%-12.8f; %-12.8f]\n', xx1(1,end),xx1(2,end)); % Analytic
fprintf('   Forward : [%-12.8f; %-12.8f]\n', xx3(1,end),xx3(2,end)); % Forward
fprintf('   Centred : [%-12.8f; %-12.8f]\n', xx4(1,end),xx4(2,end)); % Centred
fprintf('-------------------------\n') 
fprintf('SECOND SOLUTION\n')
fprintf('   Analytic : [%-12.8f; %-12.8f]\n', xx2(1,end),xx2(2,end)); % Analytic
fprintf('   Forward : [%-12.8f; %-12.8f]\n', xx5(1,end),xx5(2,end)); % Forward
fprintf('   Centred : [%-12.8f; %-12.8f]\n', xx6(1,end),xx6(2,end)); % Centred
fprintf('-------------------------\n')


%% EXERCISE 3
clearvars; close all; clc;

% 3.1) Look at Heun.m  

% 3.2) Heun: solution and comparison
f = @(x,t) -((t^2-1) / (t^2+1)) * x; % Function definition x_dot = f*x

x_an = @(t) exp(2 * atan(t) - t) ; % Analytic solution
x0 = 1 ; t0 = 0; tf = 2; % Initial condition, initial and final time

h = [0.5 0.2 0.05 0.01] ;  % [h1 h2 h3 h4]

% Heun's integration with each h (Heun.m)
TIME_Heun = [] ; % Inizialization 
ERR_Heun = [] ; 

% Solving with Heun method 
for i = 1 : length(h) 
    tic ;
    [t_h, xx_Heun] = Heun(f, x0, h(i), t0, tf); 
    Time = toc ;
    TIME_Heun = [TIME_Heun; Time] ; % Computational time 

    err = abs(xx_Heun - x_an(t_h)) ; % Comparing Heun and analytic solution
    ERR_Heun = [ERR_Heun; norm(err)] ;   % Saving (tt will be used in point 4)
    semilogy(t_h, err,'LineWidth',1.5) ; % Plot in semilog scale --> much more interesting 
    hold on; grid on;
    semilogy(t_h, err,'x r','MarkerSize', 10) ; 
end

% Plot settings
lgd = legend('$h_1 = 0.5$','', '$h_2 = 0.2$','', '$h_3 = 0.05$','', ...
             '$h_4 = 0.01$', 'Interpreter', 'latex') ;
lgd.Location = 'southeast';  lgd.FontSize = 15 ; 
set(gca,'FontSize',12)
xlabel('$t_h$', 'Interpreter', 'latex', 'FontSize', 20) ;
ylabel('$|x_{an} - x_{Heun}|$', 'Interpreter', 'latex', 'FontSize', 20) ; 
title('RK2 integration error', 'FontSize',15)


% title('Heun integration error', 'FontSize',15)

% 3.3) RK4: solution and comparison (RK4.m)
figure() ; % New plot
TIME_RK4 = [] ; % Inizialization 
ERR_RK4 = [] ;

for i = 1 : length(h) 
    tic ; 
    [t_h, xx_RK4] = RK4(f, x0, h(i), t0, tf);  
    Time = toc ;
    TIME_RK4 = [TIME_RK4; Time] ; % Computational time 

    err = abs(xx_RK4 - x_an(t_h)) ; % Comparing Heun and analytic solution
    ERR_RK4 = [ERR_RK4 norm(err)] ;

    semilogy(t_h, err,'LineWidth',1.5) ; % Plot in semilog scale --> much more interesting 
    hold on; grid on;
    semilogy(t_h, err,'x r','MarkerSize',10) ;  
end

% Plot settings
lgd = legend('$h_1 = 0.5$','', '$h_2 = 0.2$','', '$h_3 = 0.05$','', ...
             '$h_4 = 0.01$', 'Interpreter', 'latex') ;
lgd.Location = 'southeast';  lgd.FontSize = 15 ; 
set(gca,'FontSize',12)
xlabel('$t_h$', 'Interpreter', 'latex', 'FontSize', 20) ;
ylabel('$|x_{an} - x_{RK4}|$', 'Interpreter', 'latex', 'FontSize', 20)
title('RK4 integration error', 'FontSize',15)

% 3.4) Trade off between CPU time & integraiton error
% !!!!! Important note: several simulations were carried out to get a
% meaningful graph of computational times. ONLY one has been reported here
% in order to not weigh down the code !!!!!

figure() ; 
loglog(ERR_Heun, TIME_Heun, '. r', 'MarkerSize',20) ;
hold on; axis equal; 
loglog(ERR_RK4, TIME_RK4, '. b', 'MarkerSize',20) ;

legend('Heun', 'RK4','FontSize',15) ; 
set(gca,'FontSize',12)
xlabel('$||error||$', 'Interpreter', 'latex', 'FontSize', 20) ;
ylabel('CPU time [s]','FontSize', 20) ;
title('Trade off: CPU time & Integration time','FontSize', 15) ;

%% EXERCISE 4
clearvars; close all; clc;

% Data 
A = @(alpha) [0 1; -1 2*cos(alpha)] ; 

%%% Exercise 4.1) Look at the report
% F_RK2 = @(alpha, h) (eye(2) + h*A(alpha) + (h^2/2)*A(alpha)^2);

%%% Exercise 4.2) 
alpha = pi ; % Given by the problem 
A = A(alpha) ; % Substituing the known alpha

F_RK2_pi = @(h) eye(2) + h * A + (h^2/2) * A^2; % F_RK2 Operator for alpha = pi;
g = @(h) max(abs(eig(eye(2) + h * A + (h^2/2) * A^2))) - 1 ; % Functional definition 

% Plotting the functional to guess the initial value in the fzero solver
hh = linspace(-2, 4, 1000) ; 
gg = [] ; 

for h = hh
    g = @(h) max(abs(eig(eye(2) + h * A + (h^2/2) * A^2))) - 1 ; % Functional(h)
    gg = [gg g(h)] ; 
end

plot(hh, gg, 'LineWidth',1.5) ; % Functional plot
grid on; hold on;
plot(2,0,'s','LineWidth',10) ; plot(0,0,'s', 'LineWidth',10)  
plot(linspace(-2,4,100), zeros(100,1),'--k','LineWidth',1.2)

% Plot settings 
txt1 = legend('Functional $g(h)$','good initial guess','trivial solution', ...
           '$g = 0$','Interpreter','latex','FontSize',15);     
xlabel('$h$','Interpreter','latex','FontSize', 20);
ylabel('$g(h)$','Interpreter','latex','FontSize', 20);
title('Functional plot','FontSize',15) ;

h_guess = 2; % The one estimated previously (graphically)
options = optimset('Display','off') ;
h_pi = fsolve(g, h_guess, options) ; % h > 0 s.t. g = 0 <--> max(abs(eig(F_RK2))) = 1;

%%% Exercise 4.3)

% Plotting the functional to choose the starting h values 
figure() ;
gg = [] ; GG = [] ; 
hh = linspace(-5, 5, 1000) ;
for alpha = [0 : pi/4 : pi]
    A = @(alpha) [0 1; -1 2*cos(alpha)] ; 
    A = A(alpha) ;
    for h = hh 
    g = @(h) max(abs(eig(eye(2) + h * A + (h^2/2) * A^2))) - 1 ; % Functional 
    gg = [gg g(h)] ; 
    end
    GG = [GG ; gg] ; % Updating 
    gg = [] ; % Cleaning 
end

for i = 1 : size(GG,1)
    plot(hh, GG(i,:), 'LineWidth',1.5) ;
    grid on; hold on;
end
plot(linspace(-5,5,100), zeros(100,1),'--k','LineWidth',1.2) % g = 0 line
plot(-2,0,'s','LineWidth',5) ; plot(0,0,'s', 'LineWidth',5) ; plot(2,0,'s','LineWidth',5) ;
lgd = legend('$\alpha = 0$', '$\alpha = \pi/4$', '$\alpha = \pi/2$', ...
             '$\alpha = 3/4\pi$','$\alpha = \pi$', 'Interpreter','latex','FontSize',15); 
a = annotation('textarrow',[.22, .3], [.25,.6]);
a.Interpreter = 'latex' ; a.String = {'$\alpha$ increases'} ; 
a.Color = 'black' ; a.FontSize = 15 ; a.LineWidth = 1.5 ;
xlabel('$h$','Interpreter','latex','FontSize', 20);
ylabel('$g_{\alpha}(h)$','Interpreter','latex','FontSize', 20);
title('Functional($\alpha$)','Interpreter','latex','FontSize',20,'FontWeight','bold') ;
rectangle('Position',[1.5 -1 1 2],'LineStyle','--') ;
rectangle('Position',[-0.5 -1 1 2],'LineStyle','--') ;

% From previous consideration we know a correct initial guess
h_guess = [0 2] ;
hh_alpha = [] ;
LAMBDA_alpha = [] ;

for alpha = 0 : pi/128 : pi
    A = @(alpha) [0 1; -1 2*cos(alpha)] ;
    A = A(alpha) ; 
    lambda_alpha = eig(A) ;
    LAMBDA_alpha = [LAMBDA_alpha; lambda_alpha] ; 

    g = @(h) max(abs(eig(eye(2) + h * A + (h^2/2) * A^2))) - 1 ;

    if (alpha > (pi/2 - pi/64)) && (alpha < (pi/2 + pi/64))
        h_alpha = fzero(g,h_guess(1), options) ;
    else
        h_alpha = fzero(g,h_guess(2), options) ;
    end
    hh_alpha = [hh_alpha; h_alpha] ;

end

figure() ;
j = 0 ; LAMBDA_h = [] ; 
for i = 1 : 2 : length(LAMBDA_alpha)
    j = j + 1;
    lambda_h = hh_alpha(j) * LAMBDA_alpha(i:i+1) ; 
    LAMBDA_h = [LAMBDA_h; lambda_h] ;    
end

% Plotting stability margin for RK2 method 
plot(real(LAMBDA_h(1 : 2 : length(LAMBDA_h))), imag(LAMBDA_h(1 : 2 : length(LAMBDA_h))), ...
     ' b ', 'LineWidth',1.5 ) ;  
hold on ; axis equal ; grid on ; 
plot(real(LAMBDA_h(2 : 2 : length(LAMBDA_h))), imag(LAMBDA_h(2 : 2 : length(LAMBDA_h))), ...
     ' b ', 'LineWidth',1.5 ) ;  
legend('Stability domain of RK2')

ha = gca;
ha.XAxisLocation = 'origin';
ha.YAxisLocation = 'origin';

% Filling up the stable region for RK2
p3 = fill(real(LAMBDA_h(1:2:length(LAMBDA_h))), imag(LAMBDA_h(1:2:length(LAMBDA_h))),'g', 'FaceAlpha',0.2);
p3 = fill(real(LAMBDA_h(2:2:length(LAMBDA_h))), imag(LAMBDA_h(2:2:length(LAMBDA_h))),'g', 'FaceAlpha',0.2);


% Exercise 4.4) 
% We repeat the same procedure as for RK2 
A = @(alpha) [0 1; -1 2*cos(alpha)] ; 
F_RK4 = @(h) (eye(2) + h*A + (h^2/2)*A^2 + (h^3/6)*A^3 + (h^4/24)*A^4);

figure(4) ;
gg = [] ; GG = [] ; 
hh = linspace(-5, 5, 1000) ;

for alpha = [0 : pi/6 : pi]
    A = @(alpha) [0 1; -1 2*cos(alpha)] ; 
    A = A(alpha) ;
    for h = hh 
    g = @(h) max(abs(eig(eye(2) + h*A + (h^2/2)*A^2 + (h^3/6)*A^3 + (h^4/24)*A^4))) - 1 ; 
    gg = [gg g(h)] ; 
    end
    GG = [GG ; gg] ;
    gg = [] ;
end

for i = 1 : size(GG,1)
    plot(hh, GG(i,:), 'LineWidth',1.2) ; % Plotting the functionals for different alphas
    grid on; hold on;
end

plot(linspace(-5,5,100), zeros(100,1),'--k','LineWidth',1.2) % g = 0 line
plot(-2.8,0,'s','LineWidth',5) ; plot(0,0,'s', 'LineWidth',5) ; plot(2.8,0,'s','LineWidth',5) ;
lgd = legend('$\alpha = 0$', '$\alpha = \pi/6$', '$\alpha = \pi/3$', ...
             '$\alpha = \pi/2$','$\alpha = 2/3\pi$','$\alpha = 5/6\pi$','$\alpha = \pi$', 'Interpreter','latex','FontSize',15); 
xlabel('$h$','Interpreter','latex','FontSize', 20);
ylabel('$g_{\alpha}(h)$','Interpreter','latex','FontSize', 20);

rectangle('Position',[2.3 -3 1 6],'LineStyle','--') ;
rectangle('Position',[-0.5 -3 1 6],'LineStyle','--') ;

h_guess = 2.7 ; % Result of previous plot
hh_alpha = [] ;
LAMBDA_alpha = [] ;
for alpha = 0 : pi/64 : pi
    A = @(alpha) [0 1; -1 2*cos(alpha)] ;
    A = A(alpha) ; 
    lambda_alpha = eig(A) ;
    LAMBDA_alpha = [LAMBDA_alpha; lambda_alpha] ; 
    g = @(h) max(abs(eig(eye(2) + h*A + (h^2/2)*A^2 + (h^3/6)*A^3 + (h^4/24)*A^4))) - 1 ;
    h_alpha = fzero(g,h_guess, options) ;
    hh_alpha = [hh_alpha; h_alpha] ;

end

figure(3) ;
j = 0 ; LAMBDA_h = [] ; 
for i = 1 : 2 : length(LAMBDA_alpha)
    j = j + 1;
    lambda_h = hh_alpha(j) * LAMBDA_alpha(i:i+1) ; 
    LAMBDA_h = [LAMBDA_h; lambda_h] ;    
end

% Plotting RK4 margin stability 
plot(real(LAMBDA_h(1 : 2 : length(LAMBDA_h))), imag(LAMBDA_h(1 : 2 : length(LAMBDA_h))), ' k ', 'LineWidth',1.5 ) ;  
hold on ; axis equal ; grid on ; 
plot(real(LAMBDA_h(2 : 2 : length(LAMBDA_h))), imag(LAMBDA_h(2 : 2 : length(LAMBDA_h))), ' k ', 'LineWidth',1.5 ) ;
ha = gca;
ha.XAxisLocation = 'origin';
ha.YAxisLocation = 'origin';
lbl = xlabel('${Re}(\lambda\cdot h)$','Interpreter', 'latex','FontSize',20) ;
lbl.HorizontalAlignment = 'center' ; lbl.Position = [-1 -3.4 -2.3485] ;
lbl = ylabel('${Im}(\lambda \cdot h)$','Interpreter', 'latex','FontSize',20) ;
set(get(gca,'YLabel'),'Rotation',90)
lbl.HorizontalAlignment = 'center'  ; lbl.Position = [-5.3 0 -2.3485] ; 

% Filling up the stable regions 
p3 = fill(real(LAMBDA_h(1:2:length(LAMBDA_h))), imag(LAMBDA_h(1:2:length(LAMBDA_h))),'g', 'FaceAlpha',0.1);
p3 = fill(real(LAMBDA_h(2:2:length(LAMBDA_h))), imag(LAMBDA_h(2:2:length(LAMBDA_h))),'g', 'FaceAlpha',0.1);


% Point h_i*lambda exercise 3;
h_ex3 = [0.5 0.2 0.05 0.01] ; 
f = @(t) -((t^2-1) / (t^2+1)) ; % Function definition x_dot = f*x

t = linspace(0, 2, 100) ; 

% Computing the @initial and final time
hh_lambdai = [h_ex3] * f(0) ;
hh_lambdaf = [h_ex3] * f(2) ;

% Plotting in the same graph
figure(3) ; 
plot(hh_lambdai, 0, 'r o', 'LineWidth',3) ;
plot(hh_lambdaf, 0, 'b o', 'LineWidth',3) ;
legend('RK2 margin stability','' ,'RK2 stable region','', ...
        'RK4 margin stability', '','RK4 stable region','','$\lambda h(t=0)$', ...
        '','','', '$\lambda h(t=2)$','Interpreter','latex','FontSize',10,'Location', 'northwest') ; 
title('Stability regions', 'FontSize', 15)


%% EXERCISE 5
clearvars; close all; clc;

% Data
xx0 = [1; 1] ; % Initial condition
t0  = 0; tf  = 1 ; % Initial and final time
Tol = [1e-3; 1e-4; 1e-5; 1e-6] ; % Tollerances vector
alpha_vec = [0: pi/180 : pi] ; % Alpha span vector 

%%% Exercise 5.1: RK1 
% Initialization for the tollerances requested 
h = [] ;
theta = 0 ; 
for i = 1 : length(Tol)
    A = [0, 1; -1, 2*cos(theta)] ;
    x_anal = expm(A*tf)*xx0; % Analytic solution        
    g = @(h) norm((x_anal-(eye(2)+h*(1/h-(floor(1/h)))*A)*((eye(2)+h*A)^(floor(1/h)))*xx0)) - Tol(i) ; 
    h(i) = fzero(g, [eps 1]) ; % First h guess
end

h_guess = zeros(length(alpha_vec)+1, 4) ; % Creating the matrix 
h_guess(1, :) = h; % Updating the first guess
lambda = zeros(length(alpha_vec), 2) ; 
lambda_h = zeros(length(alpha_vec), 2) ; % Creating the matrix
hh = zeros(length(alpha_vec), length(Tol)) ;

% setting properties
options = optimoptions('fsolve', 'Display', 'off', 'Tolfun', 1e-12) ;
graph = {'k', 'g', 'b', 'r'} ;  % Plot colors

% Cycling for each alpha between 0 and pi 
figure();
for i = 1 : length(Tol)
    for j = 1 : length(alpha_vec)
        A = [0, 1; -1, 2*cos(alpha_vec(j))] ;   
        x_anal = expm(A*tf)*xx0 ;     
        lambda(j, :) = eig(A)' ;
        g = @(h) max(abs(x_anal-(eye(2)+rem(1,h)*A)*((eye(2)+h*A)^(floor(1/h)))*xx0))-Tol(i) ; 
        hh(j, i) = fsolve(g, h_guess(j, i), options) ;
        h_guess(j+1, i) = hh(j, i) ;        
        lambda_h(j, :) = hh(j, i)*lambda(j, :) ; % The quantity we want to plot
    end 
    % Plot of RK1
    plot_vec = zeros(2*length(alpha_vec), 1) ;
    plot_vec(1: length(alpha_vec)) = lambda_h(:, 1) ;
    plot_vec(length(alpha_vec)+1: end) = flip(lambda_h(:, 2),1) ;
    plot(plot_vec, 'color', graph{i}, 'LineWidth', 1.7)
    hold on;
end

% Settings 
axis equal; 
grid on; 
xlabel('${Re}(\lambda\cdot h)$','Interpreter', 'latex','FontSize',20) ;
ylabel('${Im}(\lambda \cdot h)$','Interpreter', 'latex','FontSize',20) ;
legend('$Tol = 10^{-3}$', '$Tol = 10^{-4}$', '$Tol = 10^{-5}$', ...
       '$Tol = 10^{-6}$','Interpreter', 'Latex', 'fontsize', 20) ;

% Function evaluations 
RK1_eval = zeros(length(Tol), 1) ;
hh = zeros(length(Tol), 1) ;
theta = pi ; 
for i = 1: length(Tol)    
    A = [0, 1; -1, 2*cos(theta)] ;
    x_anal = expm(A*tf)*xx0 ;
    lambda = eig(A)' ; % Transposing the eig vector     
    g = @(h) norm((x_anal - (eye(2) + rem(1, h)*A)*((eye(2)+h*A)^(floor(1/h)))*xx0), inf) - Tol(i) ;
    [hh(i), ~, ~, ~] = fsolve(g, h(i), options) ;
    RK1_eval(i) = tf / hh(i) ;
end

%%% Exercise 5.3): RK2
% Initialize step guess for every tolerance
h = zeros(1, length(Tol));
for i = 1: length(Tol)
    A    = [0, 1; -1, 2*cos(0)] ;
    x_anal = expm(A*tf)*xx0 ;         
    g    = @(h) max(abs(x_anal-(eye(2)+h*(1/h-(floor(1/h)))*A+(h*(1/h-(floor(1/h)))*A)^2/2)* ...
                        ((eye(2)+h*A+(h*A)^2/2)^(floor(1/h)))*xx0))-Tol(i) ; 
    
    h(i) = fzero(g, [eps 1]) ;
end
hh = zeros(length(alpha_vec), length(Tol)) ;
lambda = zeros(length(alpha_vec), 2) ;
lambda_h = zeros(length(alpha_vec), 2) ;
h_guess = zeros(length(alpha_vec)+1, length(Tol)) ;
h_guess(1, :) = h ;

figure();
for i = 1: length(Tol)
    for j = 1:length(alpha_vec)
        
        A = [0, 1; -1, 2*cos(alpha_vec(j))] ;
        lambda(j, :) = eig(A)' ;
        x_anal = expm(A*tf)*xx0 ;    
        g = @(h) norm((x_anal-(eye(2)+rem(1, h)*A+(rem(1, h)*A)^2/2)* ...
                        ((eye(2)+h*A+(h*A)^2/2)^(floor(1/h)))*xx0), inf)-Tol(i);
        hh(j, i) = fsolve(g, h_guess(j, i), options) ;
        h_guess(j+1, i) = hh(j, i) ; 
        lambda_h(j, :) = hh(j, i)*lambda(j, :) ;  
    end 
    
    % Vector correction
    plot_vec = zeros(2*length(alpha_vec), 1) ;
    plot_vec(1: length(alpha_vec))= lambda_h(:, 1) ;
    plot_vec(length(alpha_vec)+1: end) = flip(lambda_h(:, 2)) ;
    plot(plot_vec, 'color', graph{i}, 'LineWidth', 1.7)
    hold on;
end
axis equal; grid on;
xlabel('${Re}(\lambda\cdot h)$','Interpreter', 'latex','FontSize',20) ;
ylabel('${Im}(\lambda \cdot h)$','Interpreter', 'latex','FontSize',20) ;
legend('$Tol = 10^{-3}$', '$Tol = 10^{-4}$', '$Tol = 10^{-5}$', ...
        '$Tol = 10^{-6}$', 'Interpreter', 'Latex', 'fontsize', 20);

% Function evaluations 
RK2_eval = zeros(length(Tol), 1) ; 
hh = zeros(length(Tol), 1) ;
for i = 1: length(Tol)    
    A = [0, 1; -1, 2*cos(pi)] ;
    lambda = eig(A)' ;
    x_anal = expm(A*tf)*xx0;     
    g = @(h) norm((x_anal-(eye(2)+rem(1, h)*A+(rem(1, h)*A)^2/2)* ...
                  ((eye(2)+h*A+(h*A)^2/2)^(floor(1/h)))*xx0), inf)-Tol(i);
    [hh(i), ~, ~, ~] = fsolve(g, h(i), options) ;
    RK2_eval(i) = 2* tf / hh(i) ;
end

%%% Exercise 5.3): RK4

h = zeros(1, length(Tol)) ;
theta = 0 ;
for i = 1: length(Tol)
    
    A = [0, 1; -1, 2*cos(theta)] ;
    x_anal = expm(A*tf)*xx0 ;        
    g = @(h) max(abs(x_anal-(eye(2)+rem(1, h)*A+(rem(1, h)*A)^2/2+(rem(1, h)*A)^3/6+(rem(1, h)*A)^4/24)*...
             ((eye(2)+h*A+(h*A)^2/2+(h*A)^3/6+(h*A)^4/24)^(floor(1/h)))*xx0))-Tol(i) ; 
    h(i) = fzero(g, [eps 1]) ;
end

hh = zeros(length(alpha_vec), length(Tol));
lambda = zeros(length(alpha_vec), 2) ;
lambda_h = zeros(length(alpha_vec), 2) ;
h_guess = zeros(length(alpha_vec)+1, 4) ;
h_guess(1, :) = h ;

figure();
for i = 1: length(Tol)
    for j = 1:length(alpha_vec)
        
        A = [0, 1; -1, 2*cos(alpha_vec(j))] ;
        lambda(j, :) = eig(A)' ;
        x_anal = expm(A*tf)*xx0;    
        g = @(h) max(abs(x_anal-(eye(2)+rem(1, h)*A+(rem(1, h)*A)^2/2+(rem(1, h)*A)^3/6+(rem(1, h)*A)^4/24)*...
                 ((eye(2)+h*A+(h*A)^2/2+(h*A)^3/6+(h*A)^4/24)^(floor(1/h)))*xx0))-Tol(i) ; 
        hh(j, i) = fsolve(g, h_guess(j, i), options) ;
        h_guess(j+1, i) = hh(j, i) ; 
        lambda_h(j, :) = hh(j, i)*lambda(j, :) ;  
    end    
    
    % Vector correction
    plot_vec = zeros(2*length(alpha_vec), 1) ;
    plot_vec(1: length(alpha_vec)) = lambda_h(:, 1) ;
    plot_vec(length(alpha_vec)+1: end) = flip(lambda_h(:, 2)) ;
    plot(plot_vec, 'color', graph{i}, 'LineWidth', 1.7)
    hold on;
end
% Plot settings 
axis equal ; grid on ; 
xlabel('${Re}(\lambda\cdot h)$','Interpreter', 'latex','FontSize',20) ;
ylabel('${Im}(\lambda \cdot h)$','Interpreter', 'latex','FontSize',20) ;
legend('$Tol = 10^{-3}$', '$Tol = 10^{-4}$', '$Tol = 10^{-5}$', ...
        '$Tol = 10^{-6}$', 'Interpreter', 'Latex', 'fontsize', 20) ;

% Function evaluations 
RK4_eval = zeros(length(Tol), 1) ;
hh = zeros(length(Tol), 1) ;
theta = pi ; 
for i = 1: length(Tol)    
    A = [0, 1; -1, 2*cos(theta)] ;
    lambda = eig(A)' ;
    x_anal = expm(A*tf)*xx0 ;     
    g = @(h) norm((x_anal-(eye(2)+h*(1/h-(floor(1/h)))*A+(h*(1/h-(floor(1/h)))*A)^2/2+ ...
                   (h*(1/h-(floor(1/h)))*A)^3/6+(h*(1/h-(floor(1/h)))*A)^4/24)* ...
                   ((eye(2)+h*A+(h*A)^2/2+(h*A)^3/6+(h*A)^4/24)^(floor(1/h)))*xx0), inf) - Tol(i);
    [hh(i), ~, ~, ~] = fsolve(g, h(i), options) ;
    RK4_eval(i) = 4 * tf /hh(i) ;
end

% Plotting the number of function evaluations wrt the tollerance
figure() ; 
loglog(Tol, RK1_eval,'r' ,'LineWidth', 1.5) ; % log scale as it's more interesting 
grid on; hold on; 
loglog(Tol, RK2_eval,'b','LineWidth', 1.5) ;
loglog(Tol, RK4_eval,'g','LineWidth', 1.5) ;

% Settings 
legend('RK1', 'RK2', 'RK4', 'Interpreter', 'Latex', 'fontsize', 20)
xlabel('Tol', 'Interpreter', 'Latex', 'fontsize', 20) ;
ylabel('Function Evaluations', 'Interpreter', 'Latex', 'fontsize', 20)



%% EXERCISE 6
clearvars; close all; clc;

%%% Exercise 6.1)
A = @(alpha) [0 1; -1 2*cos(alpha)];

B_BI2 = @(h, alpha) inv(eye(2) - A(alpha)*(1-theta)*h + 1/2*(A(alpha)*(1-theta)*h)^2 ) ...
                           *(eye(2) + A(alpha)*theta*h + 1/2*(A(alpha)*theta*h)^2 );

% Substituting the first given theta
theta = 0.4 ; 
B_BI2_04 = @(h, alpha) inv(eye(2) - A(alpha)*(1-theta)*h + 1/2*(A(alpha)*(1-theta)*h)^2 ) ...
                           *(eye(2) + A(alpha)*theta*h + 1/2*(A(alpha)*theta*h)^2 );

%%%Exercise 6.2) 
% Building the functional looking for initial h values 
gg = [] ; GG = [] ; 
hh = linspace(-15, 15, 1000) ;
for alpha = 0 : pi/6 : pi
    A = @(alpha) [0 1; -1 2*cos(alpha)] ; 
    A = A(alpha) ;
    for h = hh 
    g = @(h) max(abs(eig(inv(eye(2) - A*(1-theta)*h + 1/2*(A*(1-theta)*h)^2 ) ...
                           *(eye(2) + A*theta*h + 1/2*(A*theta*h)^2 )))) - 1 ; 
    gg = [gg g(h)] ; 
    end
    GG = [GG ; gg] ;
    gg = [] ;
end

% Plotting the functional function of h, for different alpha
for i = 1 : size(GG,1)
    plot(hh, GG(i,:), 'LineWidth',1.5) ;
    grid on; hold on;
end

% h inizial conditions zones 
rectangle('Position',[6 -0.5 6 1],'LineStyle','--') ;
rectangle('Position',[-2.5 -0.5 5 1],'LineStyle','--') ;

lgd = legend('$\alpha = 0$', '$\alpha = \pi/6$','$\alpha = \pi/3$','$\alpha = \pi/2$', ...
             '$\alpha = 2/3\pi$','$\alpha = 5/6\pi$','$\pi$','Interpreter','latex') ;
lgd.FontSize = 15 ;
set(gca,'FontSize',12)
xlabel('$h$','Interpreter','latex','FontSize',20);
ylabel('$g(h)$','Interpreter','latex','FontSize',20);

title('Functional for BI2_{0.4}','FontSize',15) ;

options = optimset('Display','off') ;

% Deriving al the h through the fzero function  
h_guess = [0, 8] ; % From previous considerations
hh_alpha = [] ;
LAMBDA_alpha = [] ;
for alpha = 0 : pi/60 : pi
    A = @(alpha) [0 1; -1 2*cos(alpha)] ;
    A = A(alpha) ; 
    lambda_alpha = eig(A) ;
    LAMBDA_alpha = [LAMBDA_alpha; lambda_alpha] ; 

    g = @(h) max(abs(eig(inv(eye(2) - A*(1-theta)*h + 1/2*(A*(1-theta)*h)^2 ) ...  % Functional 
                           *(eye(2) + A*theta*h + 1/2*(A*theta*h)^2 )))) - 1 ;


    if (alpha > (pi/2 - pi/60)) && (alpha < (pi/2 + pi/60)) % Near pi/2 the initial guess is 0
        h_alpha = fzero(g,h_guess(1), options) ;
    else
        h_alpha = fzero(g,h_guess(2), options) ;
    end
    hh_alpha = [hh_alpha; h_alpha] ; % Building the h vector 

end

% Lambda*h (intermediate step for the plots) 
j = 0 ; LAMBDA_h = [] ; 
for i = 1 : 2 : length(LAMBDA_alpha)
    j = j + 1;
    lambda_h = hh_alpha(j) * LAMBDA_alpha(i:i+1) ; 
    LAMBDA_h = [LAMBDA_h; lambda_h] ;    
end

% Plotting the stability margin for BI2_0.4
figure() ;
plot(real(LAMBDA_h(1 : 2 : length(LAMBDA_h))), imag(LAMBDA_h(1 : 2 : length(LAMBDA_h))), ' k ', 'LineWidth',1.5 ) ;  
hold on ; axis equal ; grid on ;
plot(real(LAMBDA_h(2 : 2 : length(LAMBDA_h))), imag(LAMBDA_h(2 : 2 : length(LAMBDA_h))), ' k ', 'LineWidth',1.5 ) ;
hold on ; axis equal ; grid on ; 

% Filling up the UNSTABLE region
p3 = fill(real(LAMBDA_h(1:2:length(LAMBDA_h))), imag(LAMBDA_h(1:2:length(LAMBDA_h))),'r', 'FaceAlpha',0.1);
p3 = fill(real(LAMBDA_h(2:2:length(LAMBDA_h))), imag(LAMBDA_h(2:2:length(LAMBDA_h))),'r', 'FaceAlpha',0.1);

legend('Stability marging of BI2_{0.4}','','UNSTABLE region','FontSize',12)
title('Stability domains of backinterpolation θ=0.4 method','FontSize',15)

% Plot settings
ha = gca;
ha.XAxisLocation = 'origin'; % Centering the x,y at zero point
ha.YAxisLocation = 'origin';
lbl = xlabel('${Re}(\lambda\cdot h)$','Interpreter', 'latex','FontSize',20) ;
lbl.HorizontalAlignment = 'center' ; lbl.Position = [5 -6.3 -2.3485] ;
lbl = ylabel('${Im}(\lambda \cdot h)$','Interpreter', 'latex','FontSize',20) ;
set(get(gca,'YLabel'),'Rotation',90)
set(gca,'FontSize',12)
lbl.HorizontalAlignment = 'center'  ; lbl.Position = [-2.8 0 -2.3485] ;

%%% Exercise 6.3) 
theta = [0.1 0.3 0.7 0.9] ; 

% Looking for initial h guesses by plotting functional 
figure()
gg = [] ; 
hh = linspace(-10, 10, 1000) ;
for i = 1 : length(theta)
    GG = [] ; 
    for alpha = 0 : pi/6 : pi % Note: we could have smaller stepsize but the graphs would not be clear anymore
        A = @(alpha) [0 1; -1 2*cos(alpha)] ; 
        A = A(alpha) ;
        for h = hh 
        g = @(h) max(abs(eig(inv(eye(2) - A*(1-theta(i))*h + 1/2*(A*(1-theta(i))*h)^2 ) ...
                               *(eye(2) + A*theta(i)*h + 1/2*(A*theta(i)*h)^2 )))) - 1 ; 
        gg = [gg g(h)] ; 
        end
        GG = [GG ; gg] ;
        gg = [] ; 
    end
    % Subplotting the four functional for the 7 alpha values
    for j = 1:size(GG,1)
    subplot(2,2,i)
    plot(hh, GG(j,:), 'LineWidth',1) ;
    grid on ; hold on ; 
    plot(linspace(-10,10,100),zeros(1,100),'-- k','LineWidth',1.3) ;
    end
    % Titling and labeling the subplots
    title(theta(i),'FontSize',15 );
    xlabel('$h$','Interpreter','latex','FontSize',15);
    ylabel('$g(h)$','Interpreter','latex','FontSize',15);
end

% Deriving all the h through the fzero function  
figure() ;
h_guess = [0 3; 0 5; 0 4; 0 2 ] ; % From previous considerations
color_chat = {'r', 'b', 'k', 'g'} ; % Settings for plot

% Big cycle: for each theta requested
for i = 1 : length(theta)
    hh_alpha = [] ;
    LAMBDA_alpha = [] ;
    for alpha = 0 : pi/70 : pi
        A = @(alpha) [0 1; -1 2*cos(alpha)] ;
        A = A(alpha) ; 
        lambda_alpha = eig(A) ;
        LAMBDA_alpha = [LAMBDA_alpha; lambda_alpha] ; 
    
        g = @(h) max(abs(eig(inv(eye(2) - A*(1-theta(i))*h + 1/2*(A*(1-theta(i))*h)^2 ) ...  % Functional 
                               *(eye(2) + A*theta(i)*h + 1/2*(A*theta(i)*h)^2 )))) - 1 ;

        if (alpha > (pi/2 - pi/65)) && (alpha < (pi/2 + pi/65)) % Near pi/2 the initial guess is 0
            h_alpha = fzero(g,h_guess(i,1), options) ;
        else
            h_alpha = fzero(g,h_guess(i,2), options) ;
        end
        hh_alpha = [hh_alpha; h_alpha] ; % Building the h vector
    end
    
    % Lambda*h (intermediate step for the plots) 
    j = 0 ; LAMBDA_h = [] ; 
    for k = 1 : 2 : length(LAMBDA_alpha)
        j = j + 1;
        lambda_h = hh_alpha(j) * LAMBDA_alpha(k:k+1) ; 
        LAMBDA_h = [LAMBDA_h; lambda_h] ;    
    end
    
    % Plotting the stability margin for BI2_0.4
    plot(real(LAMBDA_h(1 : 2 : length(LAMBDA_h))), imag(LAMBDA_h(1 : 2 : length(LAMBDA_h))),'Color', color_chat{i}, 'LineWidth',1.5 ) ;  
    hold on ; axis equal ; grid on ;
    plot(real(LAMBDA_h(2 : 2 : length(LAMBDA_h))), imag(LAMBDA_h(2 : 2 : length(LAMBDA_h))),'Color', color_chat{i},'LineWidth',1.5 ) ;
    hold on ; axis equal ; grid on ; 

end
legend('\theta = 0.1','','\theta = 0.3','','\theta = 0.7','','\theta = 0.9','FontSize', 15, 'LineWidth',1.5)
% Plot settings
ha = gca;
ha.XAxisLocation = 'origin'; % Centering the x,y at zero point
ha.YAxisLocation = 'origin';
lbl = xlabel('${Re}(\lambda\cdot h)$','Interpreter', 'latex','FontSize',20) ;
lbl.HorizontalAlignment = 'center' ; lbl.Position = [0 -4.6 -2.3485] ;
lbl = ylabel('${Im}(\lambda \cdot h)$','Interpreter', 'latex','FontSize',20) ;
set(get(gca,'YLabel'),'Rotation',90)
set(gca,'FontSize',12)
lbl.HorizontalAlignment = 'center'  ; lbl.Position = [-5.6 0 -2.3485] ;


%% EXERCISE 7
clearvars; close all; clc;

% Data 
B = [-180.5 219.5; 179.5 -220.5] ; % xdot = Bx
xx0 = [1; 1] ; % Initial condition x(0) = xx0
t0 = 0; tf = 5; % Initial and final time of integration

%%% Exercise 7.1)
h = 0.1 ; % Time step assigned 

[t_h, xx_RK4] = RK4M(B, xx0, t0, tf, h) ; % Integrating: RK4 (B-input version)

%%% Exercise 7.2) 
xx_IEX4 = IEX4(B, xx0, t0, tf, h) ; % Integrating: IEX4

%%% Exercise 7.3) 
xx_an = @(t) expm(B*t) * xx0 ; % Exact solution

% Evaluating the analytic solution at each node
k = 0 ;  
for i = t_h
    k = k + 1; 
    xx_an_nod(: , k) = xx_an(i) ; 
end

% Difference between the analytic solution and the approximated 
err_RK4 = abs(xx_an_nod - xx_RK4) ; 
err_IEX4 = abs(xx_an_nod - xx_IEX4) ;

% Plotting the the error
% Note: It has just been considered the first component of the error as 
%       the second component behaves in the exact way
subplot(3,1,1)
plot(t_h, err_RK4(1,:),'b','LineWidth',2) ;
grid on; hold on;  
set(gca, 'YScale', 'log')
set(gca,'FontSize',12)

title('RK4 vs Analytic','FontSize',15);

subplot(3,1,2) ;
plot(t_h, err_IEX4(1,:),'r','LineWidth',2) ; 
grid on;   
set(gca, 'YScale', 'log')
set(gca,'FontSize',12)
title('IEX4 vs Analytic','FontSize',15);
ylabel('$|x_{an}-x_{approx}|$', 'Interpreter','latex','FontSize',20)
 
subplot(3,1,3); 
plot(t_h, err_RK4(1,:),'b','LineWidth',2) ; 
grid on; hold on;  
set(gca, 'YScale', 'log')
set(gca,'FontSize',12)
plot(t_h, err_IEX4(1,:),'r','LineWidth',2) ; 
title('RK4 vs IEX4','FontSize',15);
xlabel('$t$', 'Interpreter','latex','FontSize',20, 'Position', ...
       [2.5 3.900984254120122e-73 -0.999999999999986]) ;

%%% Exercise 7.4

% Plotting the region of stability of IEX4
A = @(alpha) [0 1; -1 2*cos(alpha)] ; 
F_IEX4 = @(h) -1/6*inv(eye(2)-A(alpha)*h)+4*inv(eye(2)-A(alpha)*h/2)^2-27/2*inv(eye(2)- ...
            A(alpha)*h/3)^3+32/3*inv(eye(2)-A(alpha)*h/4)^4;


h_guess = [0, 12] ; % Guessed as in the previous cases
hh_alpha = [] ; 
LAMBDA_alpha = [] ;
options = optimset('Display','off') ;

for alpha = [0 : pi/80 : pi]
    A = @(alpha) [0 1; -1 2*cos(alpha)] ;
    A = A(alpha) ; 
    lambda_alpha = eig(A) ;
    LAMBDA_alpha = [LAMBDA_alpha; lambda_alpha] ; 

    g1 = @(h) max(abs(eig(-1/6*inv(eye(2)-A*h)+4*inv(eye(2)-A*h/2)^2-27/2*inv(eye(2)- ...
            A*h/3)^3+32/3*inv(eye(2)-A*h/4)^4))) - 1 ;
    
        if (alpha > (pi/2 - pi/38)) && (alpha < (pi/2 + pi/38)) % Near pi/2 the initial guess is 0
            h_alpha = fzero(g1,h_guess(1), options) ;
        else
            h_alpha = fzero(g1,h_guess(2), options) ;
        end
        hh_alpha = [hh_alpha; h_alpha] ; % Building the h vector 

end

% Intermediate step for the plot
figure() ;
j = 0 ; LAMBDA_h = [] ;
for i = 1 : 2 : length(LAMBDA_alpha)
    j = j + 1;
    lambda_h = hh_alpha(j) * LAMBDA_alpha(i:i+1) ; 
    LAMBDA_h = [LAMBDA_h; lambda_h] ;  
end

% Plotting IEX4 margin stability 
plot(real(LAMBDA_h(1 : 2 : length(LAMBDA_h))), imag(LAMBDA_h(1 : 2 : length(LAMBDA_h))), ' k ', 'LineWidth',1.5 ) ;  
hold on ; axis equal ; grid on ; 
plot(real(LAMBDA_h(2 : 2 : length(LAMBDA_h))), imag(LAMBDA_h(2 : 2 : length(LAMBDA_h))), ' k ', 'LineWidth',1.5 ) ;

% Filling up the unstable region
p3 = fill(real(LAMBDA_h(1:2:length(LAMBDA_h))), imag(LAMBDA_h(1:2:length(LAMBDA_h))),'r', 'FaceAlpha',0.1);
p3 = fill(real(LAMBDA_h(2:2:length(LAMBDA_h))), imag(LAMBDA_h(2:2:length(LAMBDA_h))),'r', 'FaceAlpha',0.1);


% Plotting RK4 margin stability
h_guess = 2.7 ;
hh_alpha = [] ;
LAMBDA_alpha = [] ;
for alpha = [0 : pi/64 : pi]
    A = @(alpha) [0 1; -1 2*cos(alpha)] ;
    A = A(alpha) ; 
    lambda_alpha = eig(A) ;
    LAMBDA_alpha = [LAMBDA_alpha; lambda_alpha] ; 

    g = @(h) max(abs(eig(eye(2) + h*A + (h^2/2)*A^2 + (h^3/6)*A^3 + (h^4/24)*A^4))) - 1 ;
    
    h_alpha = fzero(g,h_guess, options) ;

    hh_alpha = [hh_alpha; h_alpha] ;

end

j = 0 ; LAMBDA_h = [] ; 
for i = 1 : 2 : length(LAMBDA_alpha)
    j = j + 1;
    lambda_h = hh_alpha(j) * LAMBDA_alpha(i:i+1) ; 
    LAMBDA_h = [LAMBDA_h; lambda_h] ;    
end

plot(real(LAMBDA_h(1 : 2 : length(LAMBDA_h))), imag(LAMBDA_h(1 : 2 : length(LAMBDA_h))), ' b ', 'LineWidth',1.5 ) ;  
hold on ; axis equal ; grid on ; 
plot(real(LAMBDA_h(2 : 2 : length(LAMBDA_h))), imag(LAMBDA_h(2 : 2 : length(LAMBDA_h))), ' b ', 'LineWidth',1.5 ) ;
% Filling up RK4 stable region
p3 = fill(real(LAMBDA_h(1:2:length(LAMBDA_h))), imag(LAMBDA_h(1:2:length(LAMBDA_h))),'g', 'FaceAlpha',0.1);
p3 = fill(real(LAMBDA_h(2:2:length(LAMBDA_h))), imag(LAMBDA_h(2:2:length(LAMBDA_h))),'g', 'FaceAlpha',0.1);

% Computing and plotting in the lambda*h plane the eigen values 
lambda = eig(B) ; 
plot(lambda(1)*0.1, 0, 's b','MarkerSize',10, 'LineWidth',4) ; 
plot(lambda(2)*0.1, 0, 's b','MarkerSize',10, 'LineWidth',4) ; 

% Label and settings
set(gca,'FontSize',12)
xlabel('${Re}(\lambda\cdot h)$','Interpreter', 'latex','FontSize',20) ;
ylabel('${Im}(\lambda \cdot h)$','Interpreter', 'latex','FontSize',20) ;
lgd = legend('IEX4 margin stability','','IEX4 unstable region','', 'RK4 margin stability','', ...
             'RK4 stable region', '','$\lambda_1h, \lambda_2h$','Interpreter','latex');
lgd.Location = 'northwest' ; lgd.FontSize = 15 ;



%% EXERCISE 8
clc; clear; close all

%Data
xx0 = [1; 1]; % x(t0) and y(t0) --> Initial condition
t0 = 0; tf = 3; % Initial and final time
h = 0.1; % Given timestep 
t_h = 0 : h : tf ; % Vector of times 

% Analytic solution
%%%% METTERE SOLUZIONE ANALITICA E FARE UN CONFRONTO FATTO COME DIO COMANDA

%%% Exercise 8.1) AB3 solution 
xx_AB3 = AB3(@f_ex8, xx0, t0, tf, h); 

%%% Exercise 8.2) AMR, ABM3, BDF3 solutions 
xx_AM3 = AM3(@f_ex8, xx0, t0, tf, h); 
xx_ABM3 = ABM3(@f_ex8, xx0, t0, tf, h) ; 
xx_BDF3 = BDF3(@f_ex8, xx0, t0, tf, h) ;

% Plotting the y solutions to compare them 
subplot(2,2,1)
plot(t_h, xx_AB3(2,:),'b','LineWidth',2); hold on

title('AB3','FontSize',15); 
xlabel('$t$','Interpreter','latex','FontSize',15); ylabel('$y$','Interpreter','latex','FontSize',15);
grid on; 

subplot(2,2,2)
plot(t_h, xx_AM3(2,:),'b','LineWidth',2)
title('AM3','FontSize',15)
xlabel('$t$','Interpreter','latex','FontSize',15); ylabel('$y$','Interpreter','latex','FontSize',15);
grid on; 

subplot(2,2,3)
plot(t_h, xx_ABM3(2,:),'b','LineWidth',2)
title('ABM3','FontSize',15)
xlabel('$t$','Interpreter','latex','FontSize',15); ylabel('$y$','Interpreter','latex','FontSize',15);
grid on; 

subplot(2,2,4)
plot(t_h, xx_BDF3(2,:),'b','LineWidth',2)
title('BDF3','FontSize',15)
xlabel('$t$','Interpreter','latex','FontSize',15); ylabel('$y$','Interpreter','latex','FontSize',15);
grid on; 


% Plotting the x solutions to compare them 
figure()

subplot(2,2,1)
plot(t_h, xx_AB3(1,:),'r','LineWidth',2); hold on
title('AB3','FontSize',15); 
xlabel('$t$','Interpreter','latex','FontSize',15); ylabel('$x$','Interpreter','latex','FontSize',15);
grid on; 

subplot(2,2,2)
plot(t_h, xx_AM3(1,:),'r','LineWidth',2)
title('AM3','FontSize',15)
xlabel('$t$','Interpreter','latex','FontSize',15); ylabel('$x$','Interpreter','latex','FontSize',15);
grid on; 

subplot(2,2,3)
plot(t_h, xx_ABM3(1,:),'r','LineWidth',2)
title('ABM3','FontSize',15)
xlabel('$t$','Interpreter','latex','FontSize',15); ylabel('$x$','Interpreter','latex','FontSize',15);
grid on; 

subplot(2,2,4)
plot(t_h, xx_BDF3(1,:),'r','LineWidth',2)
title('BDF3','FontSize',15)
xlabel('$t$','Interpreter','latex','FontSize',15); ylabel('$x$','Interpreter','latex','FontSize',15);
grid on; 

% x component eigenvalues 

lambda = @(t) -5/2*(1+8*sin(t)); % Time varian eigenvalue
t = linspace(0,3,1000) ; 

figure()
plot(t, 0.1*lambda(t),'k','LineWidth',2) ; 
grid on ; hold on; 
set(gca,'FontSize',12)
xlabel('$t$','Interpreter','latex','FontSize',20); 
ylabel('$h\lambda(t)$','Interpreter','latex','FontSize',20);
legend('Time-variant eigenvalue','FontSize',15);

g = yline(-0.6,'-- r'); g.DisplayName = 'AB3 limit' ; 
% g = yline(-6); g.DisplayName = 'AB3 limit'
g = yline(-1.65, '-- b'); g.DisplayName = 'ABM3 limit';
% g = yline(-6.5); g.DisplayName = 'BDF3 limit'

%% Functions 

% Exercise 1.1)
function [xx,N_it, ACC] = Bisection(f, a, b, Tol, n)
%
% [xx,it] = Bisection(f, a, b, Tol) 
%
% Bisection method: root finding method which consists on repeatedly bisecting 
%                   the interval and then selecting the subinterval in which 
%                   the function changes sign, and therefore must contain a root.
%
% INPUT:
% f         Continous function 
% a,b       Extreme of the initial interval     [1x1]
% Tol       Tolerance                           [1x1]
% n         n-digit accuracy                    [1x1]
%
% OUTPUT:
% xx        Vector containing iterates          [(N_it+1)x1]
% N_it      Number of iterations                [1x1]
% ACC       accuracy errors for each iterat.    [1xn] 

N_it=-1; % First iteration compute x0 
xx=[];   % Allocation 
ACC = [];  
err = Tol + 1 ; % Initial error to enter the while cycle 
acc = 1 ;
nmax = ceil(log2((b - a)/Tol) -1 ) ;  % Theoretical Relation (explained in report)

% Check on the left and right extreme selection
if (f (a) * f (b) > 0)
    error ('Error in the guess selection') % Aborting function
end

% Computing xx 
while (N_it <= nmax-1 && err > Tol && acc >= 10^(-n)) 
    
    N_it = N_it + 1 ; % Iteration updating 
    c = (b+a)/2; % Zero estimation 
    fc = f(c); % Evaluation of the potential root     
    if (fc == 0)
       err = 0;
    else
       err = abs(fc) ; 
    end
    xx = [xx ; c] ; 
    
    % Check to have the 8-digit accuracy 
    if N_it > 0
        acc = abs(xx(end) - xx(end-1)); 
    end
        ACC = [ACC acc] ;
    
    % Selection of the new interval         
    if (fc * f(a) > 0)
          a = c; 
    else 
          b = c; 
    end          
end 
% Print of the results 
fprintf('BISECTION METHOD RESULTS\n')
fprintf('   Root : %-12.8f \n', xx(end)); % Printing the final root
fprintf('   N. iterations : %d \n', N_it); % Printing the final root
fprintf('-------------------------\n')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exercise 1.2)
function [xx, N_it, ACC] = Secant(f, x0, x1, nmax, Tol, n)
%
% [xx, N_it] = Secant(f, x0, x1, nmax, Tol, n) 
%
% Secant method: root-finding algorithm that uses a succession of roots
%                of secant lines to approximate a root of a function.
%
% INPUT:
% f         Continous function 
% x0, x1    Initial guesses                     [1x1]
% nmax      Maximum number of iteration         [1x1]
% Tol       Tolerance                           [1x1]
% n         n-digit accuracy                    [1x1]
%
% OUTPUT:
% xx        Vector containing iterates          [(N_it+1)x1]
% N_it      Number of iterations                [1x1]
% ACC       accuracy errors for each iterat.    [1xn] 
% 

N_it = -1 ; % First iteration compute x0

% Initial error and accuracy to enter in the cycle
err = Tol + 1 ; 
acc = 1 ;
% Initialization 
xx = []; 
ACC = [] ; 

while (N_it <= nmax-1 && err > Tol && acc >= 10^(-n))

    N_it = N_it + 1 ; % Iteration updating
    
    f0 = f(x0) ; % Evaluation of the function at x0
    f1 = f(x1) ; % Evaluation of the function at x1

    y = x1 - ((x1-x0) / (f1-f0)) * f1 ; % [x0,x1] is the interval of the root

    err = abs(f(y)) ; % Computing error
    xx = [xx; y] ; % Updating the vector

    % Check to have the 8-digit accuracy 
    if N_it > 0
        acc = abs(xx(end) - xx(end-1)); 
    end
    ACC = [ACC acc] ; 

    x0 = x1 ;
    x1 = y ;
end

% Print of the results 
fprintf('SECANT METHOD RESULTS\n')
fprintf('   Root : %-12.8f \n', xx(end)); % Printing the final root
fprintf('   N. iterations : %d \n', N_it); % Printing the final root
fprintf('-------------------------\n')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exercise 1.3) 
function [xx, N_it, ACC] = Regula_Falsi(f, x0, x1, nmax, Tol, n)
%
% [xx, N_it] = Regula_Falsi(f, x0, x1, nmax, Tol, n)
%
% Regula Falsi method: root finding metod consisting in a proper
%                      combination of bisection and secant methods.
%
% INPUT
% f         Function 
% x0, x1    Extreme of the initial interval     [1x1]
% nmax      Maximum number of iteration         [1x1]
% Tol       Tolerance                           [1x1]
% n         n-digit accuracy                    [1x1]
%
% OUTPUT
% xx        Vector containing iterates          [(N_it+1)x1]
% N_it      Number of iterations                [1x1]
% ACC       accuracy error for each iterat.    [1xn] 
%

% Initial error and accuracy to enter in the cycle
err = Tol + 1 ;
acc = 1 ;
xx = [] ; % Vector allocation 
ACC = [] ; 
N_it = - 1 ; % First iteration compute x0

while (N_it <= nmax-1 && err > Tol && acc >= 10^(-n)) 
    
    N_it = N_it + 1 ; % Updating the counter

    f0 = f(x0) ; %Calculating the value of function at x0
    f1 = f(x1) ; %Calculating the value of function at x1

    y = x1 - ((x1-x0) / (f1-f0)) * f1; % [x0,x1] is the interval of the root

    err = abs(f(y)) ; % Error computation 

    xx = [xx; y] ;
    
    % Check to have the 8-digit accuracy 
    if N_it > 0
        acc = abs(xx(end) - xx(end-1)); 
    end
    ACC = [ACC acc] ; 

    f2 = f(y); % Evaluation of the root 

    if (f1) * (f2) < 1
        % The new interval is [y,x1]
        x0 = y;  
        x1 = x1;
    else
        % The new interval is [x0,y]
        x0 = x0 ; 
        x1 = y ;
    end
end

% Print of the results 
fprintf('REGULA FALSI METHOD RESULTS\n')
fprintf('   Root : %-12.8f \n', xx(end)); % Printing the final root
fprintf('   N. iterations : %d \n', N_it); % Printing the final root
fprintf('-------------------------\n')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 2) 
function [XX, N_iter, ERR] = newton(f, J, xx0, Tol, kmax)
%
% [xx,N_iter] = newton(f, J, xx0, Tol, kmax)
%
% Newton's method: zero finding root method
%
% INPUT:
% f          function                               [2x1]
% J          Jacobian of the function               [2x2]
% xx0        Initial guess                          [2x1]
% Toll       Tollerace: arrest criteria             [1x1]
% kmax       Maximum number of iterations           [1x1]
%
% OUTPUT:
% xx         Vector containing the last iteration   [2x1]
% N_iter     Number of iterations                   [1x1]
%
% Inizialization
k = 0 ; 
err = Tol + 1 ; ERR = [] ; 
xx = xx0;
XX = [] ; 
XX(:,1) = xx0 ; % Setting as first iteration 

while err >= Tol && k < kmax
    JJ = J(xx) ; % Jacobian evaluation 
    F = f(xx) ; % Function evaluation 
    delta = - JJ\F ; % \ : inverse command (more efficient than inv() ) 
    xx = xx + delta ; % new iteration
    XX = [XX xx] ; 
    err = norm(delta); ERR = [ERR; err] ;   
    k = k + 1 ; 
end

res = norm(f(XX(:,end))) ; 
N_iter = k ;
end

% Exercise 2.2)
%
% [XX, N_iter, ERR] = newton_diff_fin(f, xx0, Tol, kmax, type)
%
% Newton's method: zero finding root method
%
% INPUT:
% f          function                               [2x1]
% xx0        Initial guess                          [2x1]
% Toll       Tollerace: arrest criteria             [1x1]
% kmax       Maximum number of iterations           [1x1]
% type       forward or centred
%
% OUTPUT:
% xx         Vector containing the last iteration   [2x1]
% N_iter     Number of iterations                   [1x1]
%
function [XX, N_iter, ERR] = newton_diff_fin(f, xx0, Tol, kmax, type)

k = 0 ; 
err = Tol + 1 ; ERR = [] ; 
x = xx0; XX = [] ; 

while err >= Tol && k < kmax
     k = k + 1 ; 
    epsilon_x = [max(sqrt(eps),sqrt(eps)*max(abs(x(1)))); 0] ; % theroteical definition 
    epsilon_y = [0; max(sqrt(eps),sqrt(eps)*max(abs(x(2))))] ;
    if type == 'forward'
        df_fin = [(f(x + epsilon_x)-f(x))./epsilon_x(1), (f(x + epsilon_y) - f(x))./epsilon_y(2)]; % Forward finite differece     
    elseif type ==  'centred'
        df_fin =[(f(x+epsilon_x) - f(x-epsilon_x))./(2*epsilon_x(1)),(f(x+epsilon_y) ...
                - f(x-epsilon_y))./(2*epsilon_y(2))]; % Centred finite difference
    end
 
    F = f(x) ; % Function evaluation 
    delta = - df_fin\F ; % \ : inverse command (more efficient than inv() ) 
    x = x + delta ; % Iteration 
    XX = [XX x] ; % Updating
    err = norm(delta); 
    ERR = [ERR; err] ;   
end

% res = norm(f(XX(:,end))) ; % computed just for some analysis  
N_iter = k ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exercise 3) 
function [t_h, xx] = Heun(f, x0, h, t0, tf)
%
% [t_h, xx] = Heun(f, x0, h, t0, tf)
% 
% Heun's method:  numerical procedure for solving ODEs  
%
% INPUT:
% f         function such that x_dot = f*x           
% x0        Initial value: x(t=0) = x0              [1x1]
% h         Time step                               [1x1]
% t0, tf    Initial and final time of the IVP       [1x1]
%
% OUTPUT:
% t_h       Discretized time vecotor                [1xn]
% xx        Discretized solution                    [1xn]
%

t_h = t0 : h : tf ; % Building the discrete time vector
xx = zeros(1, length(t_h)) ; % Allocation 

xx(1) = x0; % Imposing the initial condition 

for k = 1 : length(t_h)-1

    xp = xx(k)+ h * f(xx(k),t_h(k)) ; % Predictor  
    xx(k+1) = xx(k) + h/2 * (f(xx(k),t_h(k)) + f(xp, t_h(k+1))) ; % Corrector

end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exercise 3)
function [t_h, xx] = RK4(f, x0, h, t0, tf)
%
% [t_h, xx] = RK4(f, x0, h, t0, tf)
% 
% Fourth Order Runge-Kutta: iterative methods, used in temporal discretization 
%                           for the approximate solutions of nonlinear equations.
%
% INPUT:
% f         function         
% x0        Initial value: x(t=0) = x0              [1x1]
% h         Time step                               [1x1]
% t0, tf    Initial and final time of the IVP       [1x1]
%
% OUTPUT:
% t_h       Discretized time vecotor                [1xn]
% xx        Discretized solution                    [1xn]
%

t_h = t0 : h : tf ; % Time instant vector 
xx = zeros(1, length(t_h)); % Allocation 
xx(1) = x0; % Initial condition

% Iterative cycle to compute the solution 
for k = 1 : length(t_h) - 1 

    xp1 = xx(k) + 1/2 * h * f(xx(k), t_h(k)); % First predictor
    xp2 = xx(k) + 1/2 * h * f(xp1, t_h(k)+1/2*h); % Second predictor 
    xp3 = xx(k) + h * f(xp2, t_h(k)+1/2*h); % Third predictor 
    
    % Collector 
    xx(k+1) = xx(k) + h * (1/6*f(xx(k),t_h(k)) + 1/3*f(xp1, t_h(k)+1/2*h) ...
                           + 1/3* f(xp2, t_h(k)+1/2*h) + 1/6*f(xp3, t_h(k)+h));

end
end

% Exercise 7.1)
function [t_h, xx] = RK4M(A, xx0, t0, tf, h)
%
% [t_h, xx] = RK4M(f, x0, h, t0, tf)
% 
% Fourth Order Runge-Kutta: iterative methods, used in temporal discretization 
%                           for the approximate solutions of nonlinear equations.
%
% INPUT:
% A         A: xdot = Ax                            [2x2]         
% xx0        Initial value: x(t=0) = x0              [1x1]
% h         Time step                               [1x1]
% t0, tf    Initial and final time of the IVP       [1x1]
%
% OUTPUT:
% t_h       Discretized time vecotor                [1xn]
% xx        Discretized solution                    [1xn]
%
h2 = h/2; h6 = h/6;

t_h = t0 : h : tf ; % Time instant vectorxx(:,1) = xx0;
xx(:, 1) = xx0; % Initial condition

for i = 1 : length(t_h) - 1   
    k1 = A * xx(:,i);
    k2 = A * (xx(:,i) + h2*k1);
    k3 = A * (xx(:,i) + h2*k2);
    k4 = A * (xx(:,i) + h*k3);
    xx(:,i+1) = xx(:,i) + h6 * (k1 + 2*k2 + 2*k3 + k4); % Collector

end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xx = IEX4(A, xx0, t0, tf, h)
%
% xx = IEX4(A, xx0, t0, tf, h)
% 
% Fourth Order IMPLICIT extrapolation technique
%
% INPUT:
% A         A: xdot = Ax                            [2x2]         
% xx0       Initial value: x(t=0) = x0              [1x1]
% h         Time step                               [1x1]
% t0, tf    Initial and final time of the IVP       [1x1]
%
% OUTPUT:
% xx        Discretized solution                    [1xn]
%
xx(:,1) = xx0;
options = optimset('Display','off'); % To avoid messages due to the fsolve

for i = 1 : (tf-t0)/h
    
    k_1  = fsolve(@(k1) k1 - xx(:,i) - h * A * k1, xx(:,i), options); % First predictor

    k_2a = fsolve(@(k2a) k2a - xx(:,i) - h/2 * A * k2a, xx(:,i),options ); % Second predictor
    k_2  = fsolve(@(k2) k2 - k_2a - h/2 * A * k2, xx(:,i),options);
    
    k_3a = fsolve(@(k3a) k3a - xx(:,i) - h/3 * A * k3a, xx(:,i),options ); % Third predictor 
    k_3b = fsolve(@(k3b) k3b - k_3a - h/3 * A * k3b, xx(:,i), options);
    k_3  = fsolve(@(k3) k3 - k_3b - h/3 * A * k3, xx(:,i),options) ;

    k_4a = fsolve(@(k4a) k4a - xx(:,i) - h/4 * A * k4a, xx(:,i),options); % Fourth predicord
    k_4b = fsolve(@(k4b) k4b - k_4a - h/4 * A * k4b, xx(:,i),options);
    k_4c = fsolve(@(k4c) k4c - k_4b - h/4 * A * k4c, xx(:,i),options);
    k_4  = fsolve(@(k4) k4 - k_4c - h/4 * A * k4, xx(:,i),options);

    xx(:,i+1) = -1/6*k_1 + 4*k_2 - 27/2*k_3 + 32/3*k_4; % Corrector 

end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exercise 8) differential equation to be solved with: AB3, AM3, ABM3, BDF3
function dxxdt = f_ex8(t, xx)

dxxdt = zeros(2,1);

dxxdt(1) = -5/2 * (1+8*sin(t)) * xx(1);
dxxdt(2) = (1-xx(1)) * xx(2) + xx(1);

end


function xx = AB3(f, xx0, t0, tf, h)
%
% xx = AB3(f, xx0, t0, tf, h)
% 
% Third order EXPLICIT Adams-Bashforth formulae
%
% INPUT:
% f         function to be solved                                   
% xx0       Initial value: x(t=0) = x0              [1x1]
% h         Time step                               [1x1]
% t0, tf    Initial and final time of the IVP       [1x1]
%
% OUTPUT:
% xx        Discretized solution                    [1xn]
%

h2 = h/2; h6 = h/6;
t = [] ; xx = [] ;% Inizialization 
xx(:,1) = xx0; % Adding initial condition to the solution vector 

% t_h vector 
t(1) = t0; t(tf/h+1) = tf;
for i = 2 : tf/h
    t(i) = (i-1)*h;
end

% Computing first 2 samples using RK4 method 
for i = 1 : 2
    
    k1 = f(t(i),xx(:,i)) ;  
    k2 = f(t(i)+h2, xx(:,i)+h2*k1) ;
    k3 = f(t(i)+h2, xx(:,i)+h2*k2) ;
    k4 = f(t(i)+h, xx(:,i)+h*k3);
    
    xx(:,i+1) = xx(:,i) + h6 * (k1 + 2*k2 + 2*k3 + k4);

end

% Building the vector using AB3 method 
for k = 3 : (tf-t0)/h
    xx(:,k+1) = xx(:,k) + h/12 * (23*f(t(k), xx(:,k)) - 16*f(t(k-1), xx(:,k-1))+ 5*f(t(k-2), xx(:,k-2)));
end

end


function xx = AM3(f, xx0, t0, tf, h)
%
% xx = AM3(f, xx0, t0, tf, h)
% 
% Third order IMPLICIT Adams-Moulton formulae
%
% INPUT:
% f         function to be solved                                   
% xx0       Initial value: x(t=0) = x0              [1x1]
% h         Time step                               [1x1]
% t0, tf    Initial and final time of the IVP       [1x1]
%
% OUTPUT:
% xx        Discretized solution                    [1xn]
%
t = [] ;  xx = [] ; 
xx(:,1) = xx0;
t(1) = t0;
t(tf/h+1) = tf;
options = optimset('Display','off'); % To avoid messages due to the fsolve

for i = 2 : tf/h
    t(i) = (i-1)*h;
end

for i = 1
    k1 = f(t(i), xx(:,i));
    k2 = f(t(i)+h/2, xx(:,i)+h/2*k1);
    k3 = f(t(i)+h/2, xx(:,i)+h/2*k2);
    k4 = f(t(i)+h, xx(:,i)+h*k3);
    
    xx(:,i+1) = xx(:,i) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
end

% Implicit relation to be solved
% Note: fsolve as it's vectorial 
for k = 2 : (tf-t0)/h 
    xx(:, k+1) = fsolve(@(xkp1) xkp1 - xx(:, k) - h/12 * (5*f(t(k+1), xkp1) ...
                          + 8*f(t(k), xx(:,k)) - f(t(k-1), xx(:,k-1))), xx(:,k),options);
end

end


function xx = ABM3(f, xx0, t0, tf, h)
%
% xx = ABM3(f, xx0, t0, tf, h)
% 
% Third order EXPLICIT Adams-Bashforth-Moulton predictor-corrector
%
% INPUT:
% f         function to be solved                                   
% xx0       Initial value: x(t=0) = x0              [1x1]
% h         Time step                               [1x1]
% t0, tf    Initial and final time of the IVP       [1x1]
%
% OUTPUT:
% xx        Discretized solution                    [1xn]
%
t = [] ; xx = [] ; 
xx(:,1) = xx0;
t(1) = t0;
t(tf/h+1) = tf;

for i = 2 : tf/h
    t(i) = (i-1)*h;
end

% Approximating first steps with RK4
for i = 1 : 2
    k1 = f(t(i), xx(:,i));
    k2 = f(t(i)+h/2, xx(:,i)+h/2*k1);
    k3 = f(t(i)+h/2, xx(:,i)+h/2*k2);
    k4 = f(t(i)+h, xx(:,i)+h*k3);
    
    xx(:,i+1) = xx(:,i) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
end

% Continuing with ABM3 methid 
for k = 3 : (tf-t0)/h 
    % AB3 predictor stage
    xx_p = xx(:,k) + h/12 * (23*f(t(k), xx(:,k)) - 16*f(t(k-1), xx(:,k-1))+ 5*f(t(k-2), xx(:,k-2))); 
    % AM3 corrector stage
    xx(:, k+1) = xx(:,k) + h/12 * (5*f(t(k+1),xx_p) + 8*f(t(k),xx(:,k))- f(t(k-1), xx(:,k-1))) ; 

end

end


function xx = BDF3(f, xx0, t0, tf, h)
%
% xx = BDF3(f, xx0, t0, tf, h)
% 
% Third order IMPLICIT Backward Difference Formulae
%
% INPUT:
% f         function to be solved                                   
% xx0       Initial value: x(t=0) = x0              [1x1]
% h         Time step                               [1x1]
% t0, tf    Initial and final time of the IVP       [1x1]
%
% OUTPUT:
% xx        Discretized solution                    [1xn]
% 
h2 = h/2; h6 = h/6;
t = [] ; xx = [] ; 
xx(:,1) = xx0;
t(1) = t0;
t(tf/h+1) = tf;

% Discretized time vector
for i = 2 : tf/h
    t(i) = (i-1)*h;
end

options = optimset('Display','off'); % To avoid messages due to the fsolve

% Approximating first two steps by using RK4
for i = 1 : 2
    
    k1 = f(t(i),xx(:,i)) ;  
    k2 = f(t(i)+h2, xx(:,i)+h2*k1) ;
    k3 = f(t(i)+h2, xx(:,i)+h2*k2) ;
    k4 = f(t(i)+h, xx(:,i)+h*k3);
    
    xx(:,i+1) = xx(:,i) + h6 * (k1 + 2*k2 + 2*k3 + k4);

end

% BDF3 method 
for k = 3 : (tf-t0)/h 
    % Implicit equation to be solved using fsolve matlab function 
    xx(:, k+1) = fsolve(@(xkp1) xkp1 - 18/11 * xx(:,k) + 9/11 * xx(:,k-1) - 2/11 * xx(:,k-2) - 6/11 * h * f(t(k+1),xkp1), xx(:,k),options);
end

end





    


