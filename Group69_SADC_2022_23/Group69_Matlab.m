% Spacecraft Attitude and Dynamics, A.Y. 2022/2023
%
% Project Assignment: Attitude Determination and Control of a 3U CubeSat
%
% Group 69
% Authors: Lanza Davide, Larocca Rocco, Mascelloni Matteo, Shakeel Afaq
%
% Last Update: 01/07/2023
%
% Estimated time of running: half an hour, depending on the pc. If a faster
% running is desired, please reduce the relative tollerance in the Simulink
% file (normally set at 1e-12).
%% Orbit Characterization
% Initialization
clc; 
close all;
clear;
format long g;

% Constants
orbit.mu = astroConstants(13);       % Earth's gravitational parameter [km^3/s^2]
orbit.w_earth = deg2rad(15.04/3600); % Earth's angular velocity [rad/s]
orbit.Re = astroConstants(23);       % Earth's mean radius [km]

% Orbit of the satellite
a = 8000;                     % Semi-major axis [km]
e = 0.01;                     % Eccentricity [-]
i =deg2rad(0.1613);           % Inclination [rad]
OM =deg2rad(pi/6);               % RAAN [rad]
om =deg2rad(0);               % Argument of perigee [rad]
th =deg2rad(0);               % True anomaly [rad]

T = 2*pi*sqrt(a^3/orbit.mu);  % Period of the orbit [s]
n = sqrt(orbit.mu/a^3);       % Mean angular velocity [rad/s]

% Initial state vector              
orbit.kep = [a,e,i,OM,om,th];             % Keplerian coordinates
[r0,v0] = kep2car(orbit.kep,orbit.mu);   % Cartesian coordinates
n0 = sqrt(orbit.mu/norm(r0)^3);         % Initial mean angular velocity [rad/s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Environment 

% Earth's magnetic field: mathematical model
H0 = sqrt((-29615e-9)^2+(-1728e-9)^2+(5186e-9)^2);
R3H0 = orbit.Re^3*H0;

% Uncontrolled initial conditions
wx0 = 0.1; wy0 = 0.1; wz0 = 0.1;
w0 = [wx0,wy0,wz0]';             % Initial angular velocity vector [rad/s] 
A_NB_0 = eye(3);
A_BN_0 = A_NB_0';                % Initial DCM (inertial to body)

%% Components Characterization
% Initialization
clc; 
close all; 
format long g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensors

% Earth Horizon sensor: CubeSense N First Gen
hs.FOV = 90;                     % [+/- deg] 
hs.rate = 2;                     % [Hz]
hs.sampling_time = 1/hs.rate;
hs.accuracy = 0.2;               % < 0.2° 3Sigma

% Magnetic sensor: SpaceMag-Lite Three-axis Magnetometer
mag.rate = 1;                   % [Hz]  
mag.sampling_time = 1/mag.rate;
mag.saturation = 60e-6;         % [T]

a_x = [1,deg2rad(rand(1)),deg2rad(rand(1))];
a_y = [deg2rad(rand(1)),1,deg2rad(rand(1))];
a_z = [deg2rad(rand(1)),deg2rad(rand(1)),1];
mag.A_error = [a_x;a_y;a_z];     % Error matrix of the sensor

mag.noise = 50e-12*[1,1,1]; 

% Gyroscope: STIM300 
gyro.rate = 20;                    % [samples/sec] 
% (Real one: 2000 samples/s. Reduced for sppeding the running)
gyro.sampling_time = 1/gyro.rate;
gyro.ARW = 0.15;                   % [°/sqrt(h)]
gyro.bias = 0.3;                   % [°/h]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actuators

% 3x Reaction wheels: RocketLab-10mNms RW-0.01
rw.A = eye(3);
rw.momentum = 18e-3;   % [Nms]
rw.maxtorque = 1e-3;   % [Nm]
rw.mass = 120e-3;      % [kg]

% 1x Inertia wheel: CubeWheel Small CW0017 
inertia.rpm = 6000;                        % [rpm] 
inertia.w = [0,0,inertia.rpm*2*pi/60];
inertia.R = (28e-3)/2;                     % [m]
inertia.R_in = (28e-3)/2;                  % [m]
inertia.h = 26.2e-3;                       % [m]
inertia.V_in=pi*inertia.R_in^2*inertia.h;  % [m^3]
inertia.mass_in = 60e-3;                   % [kg]
inertia.rho=inertia.mass_in/inertia.V_in;


inertia.Ix = 1/12*inertia.mass_in*inertia.h^2+1/4*inertia.mass_in*inertia.R_in^2;
inertia.Iy = 1/12*inertia.mass_in*inertia.h^2+1/4*inertia.mass_in*inertia.R_in^2;
inertia.Iz = 1/2*inertia.mass_in*inertia.R_in^2;
inertia.I = diag([inertia.Ix,inertia.Iy,inertia.Iz]);

%% Spacecraft Characterization
% Initialization
clc; 
close all; 
format long g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spacecraft (S/C) and Solar Panels (SP) models
sc.m_sc =3.6;      % S/C Mass [kg]
sc.w = 0.1;        % Width    [m]
sc.d = 0.1;        % Depth    [m]
sc.h = 0.300;      % Height   [m]
m = [0.01 0.05 0.01]';    % Magnetic induttance [Am^2]

sc.msp = 0.400;    % SP Mass      [kg]
sc.th = 1.5e-3;    % SP Thickness [m]

sc.m = sc.m_sc + rw.mass*3 + inertia.mass_in; % Total Mass [kg]

sc.rcg = ([0,0,0]*sc.m + [sc.h/2-sc.th/2,sc.w/2+sc.h/2,0]*sc.msp/2 + ...
          +[sc.h/2-sc.th/2,-(sc.w/2+sc.h/2),0]*sc.msp/2)/(sc.m+sc.msp);
         % Centre of gravity position in body frame
sc.rcg_sp1 = [sc.h/2;sc.w/2+sc.h/2;0];    % CG SP1
sc.rcg_sp2 = [sc.h/2;-(sc.w/2+sc.h/2);0]; % CG SP2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inertia matrix computation

% Spacecraft components
sc.Ix = 1/12*sc.m*(sc.w^2+sc.d^2);
sc.Iy = 1/12*sc.m*(sc.d^2+sc.h^2);
sc.Iz = 1/12*sc.m*(sc.h^2+sc.w^2);

% Solar Panels components
sc.Ix_sp = 1/12*sc.msp/2*(sc.d^2+sc.h^2);
sc.Iy_sp = 1/12*sc.msp/2*(sc.d^2+sc.th^2);
sc.Iz_sp = 1/12*sc.msp/2*(sc.h^2+sc.th^2);
sc.I_sp1 = diag([sc.Ix_sp sc.Iy_sp sc.Iz_sp])+ ...
           + sc.msp/2*(norm(sc.rcg_sp1)*eye(3)*kron(sc.rcg_sp1,sc.rcg_sp1'));
sc.I_sp2 = diag([sc.Ix_sp sc.Iy_sp sc.Iz_sp])+ ...
           + sc.msp/2*(norm(sc.rcg_sp2)*eye(3)*kron(sc.rcg_sp2,sc.rcg_sp2'));

sc.I = diag([sc.Ix,sc.Iy,sc.Iz])+sc.I_sp1+sc.I_sp2; % Inertia matrix I
sc.I_inv = inv(sc.I);                               % Inverse of I

%% Mission Intervals
% Initialization
clc; 
close all; 
format long g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First phase: Detumbling
t.detumbling = 300; % [s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second phase: Slew manoeuvre
t.slew = 1200;      % [s]

%% Data from Simulation 
% Initialization
clc; 
close all; 
format long g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running the simulation and retreiving data

open("Group69_Simulink.slx"); 
out = sim('Group69_Simulink');

time = out.tout;           % Time vector

% Dynamics outputs
w = squeeze(out.w.Data);   % Angular velocity vector

% Kinematics outputs
A_BN = out.A_BN.Data;      % DCM (from inertial to body)
A_NL = out.A_NL.Data;      % DCM (from inertial to LVLH)

% Orbit propagation outputs
orbit.r = out.r.Data;      % Position vector in LVLH frame

% Environment outputs
M_GG = squeeze(out.T_gravity.Data);           % Gravity Gradient Torque
M_mag = squeeze(out.T_magnetic.Data);         % Magnetic Torque

% Control outputs
A_BN_m = out.A_BN_m.Data;       % Measured DCM
angerr = out.ang_err.Data;      % Pointing Error 
T_c = out.T_c.Data;             % Reaction Wheel Torque

 
%% Plots Execution
% Initialization
clc; 
close all; 
format long g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET-UP
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesLabelFontSize', 1.4);
set(groot,'defaultAxesFontSize', 18);
set(groot,'DefaultAxesTitleFontSize',1.4)
set(groot,'defaultLegendFontSize',1.3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find phases times in the vector time
dtm = find(time==t.detumbling);
slw = find(time==t.slew);

% Actual vs Measured DCM (to run manually!)
% res11 = zeros(1,length(time)-slw);
% res12 = zeros(1,length(time)-slw);
% res13 = zeros(1,length(time)-slw);
% res21 = zeros(1,length(time)-slw);
% res22 = zeros(1,length(time)-slw);
% res23 = zeros(1,length(time)-slw);
% res31 = zeros(1,length(time)-slw);
% res32 = zeros(1,length(time)-slw);
% res33 = zeros(1,length(time)-slw);
% 
% for i = 1:(length(time)-slw)
%     res11(i) = A_BN_m(1,1,i+slw)-A_BN(1,1,i+slw);
%     res12(i) = A_BN_m(1,2,i+slw)-A_BN(1,2,i+slw);
%     res13(i) = A_BN_m(1,3,i+slw)-A_BN(1,3,i+slw);
%     res21(i) = A_BN_m(2,1,i+slw)-A_BN(2,1,i+slw);
%     res22(i) = A_BN_m(2,2,i+slw)-A_BN(2,2,i+slw);
%     res23(i) = A_BN_m(2,3,i+slw)-A_BN(2,3,i+slw);
%     res31(i) = A_BN_m(3,1,i+slw)-A_BN(3,1,i+slw);
%     res32(i) = A_BN_m(3,2,i+slw)-A_BN(3,2,i+slw);
%     res33(i) = A_BN_m(3,3,i+slw)-A_BN(3,3,i+slw);
% 
%     plot(time((slw+1):end),res11,'LineWidth',1.5)
%     hold on
%     plot(time((slw+1):end),res12,'LineWidth',1.5)
%     hold on
%     plot(time((slw+1):end),res13,'LineWidth',1.5)
%     hold on
%     plot(time((slw+1):end),res21,'LineWidth',1.5)
%     hold on
%     plot(time((slw+1):end),res22,'LineWidth',1.5)
%     hold on
%     plot(time((slw+1):end),res23,'LineWidth',1.5)
%     hold on
%     plot(time((slw+1):end),res31,'LineWidth',1.5)
%     hold on
%     plot(time((slw+1):end),res32,'LineWidth',1.5)
%     hold on
%     plot(time((slw+1):end),res33,'LineWidth',1.5)
%     grid on
%     grid minor
%     xlabel('Time [s]')
%     ylabel('Estimation Error [-]')
%     
% end

% Pointing error
figure(1)
plot(time(slw:end),angerr(slw:end),'LineWidth',1)
grid on
grid minor
xlabel('Time [s]')
ylabel('Pointing Error [deg]')
xlim([time(slw) time(end)])
ylim([-0.05 0.1])

% Angular Velocity
figure(2)
subplot(2,3,[1 2 3]) 
plot(time,w(1,:),time,w(2,:),time,w(3,:),'LineWidth',1.5)
xlim([0 time(end)])
grid on
grid minor
xlabel('Time [s]')
ylabel('Angular velocity [rad/s]')
legend('$\omega_x$','$\omega_y$','$\omega_z$')
subplot(2,3,4)
plot(time(1:dtm),w(1,1:dtm),time(1:dtm),w(2,1:dtm),time(1:dtm),w(3,1:dtm),'LineWidth',1.5)
xlim([time(1) time(dtm)])
xlabel('Time [s]')
ylabel('Angular velocity [rad/s]')
legend('$\omega_x$','$\omega_y$','$\omega_z$')
% title('De-tumbling')
grid on
grid minor
subplot(2,3,5)
plot(time(dtm:slw),w(1,dtm:slw),time(dtm:slw),w(2,dtm:slw),time(dtm:slw),w(3,dtm:slw),'LineWidth',1.5)
xlim([time(dtm) time(slw)])
xlabel('Time [s]')
ylabel('Angular velocity [rad/s]')
legend('$\omega_x$','$\omega_y$','$\omega_z$')
% title('Slew maneuvre')
grid on
grid minor
subplot(2,3,6)
plot(time(slw:end),w(1,slw:end),time(slw:end),w(2,slw:end),time(slw:end),w(3,slw:end),'LineWidth',1.5)
xlim([time(slw) time(end)])
xlabel('Time [s]')
ylabel('Angular velocity [rad/s]')
legend('$\omega_x$','$\omega_y$','$\omega_z$')
% title('Pointing')
grid on
grid minor

% Environment Torques
figure(3)
plot(time,vecnorm(M_GG),'LineWidth',1.5)
hold on
plot(time,vecnorm(M_mag),'LineWidth',1.5)
xlim([0 time(end)])
legend('$M_{GG}$','$M_{mag}$')
grid on
grid minor
xlabel('Time [s]')
ylabel('Environmental Torques [Nm]')

% Actuator torque
figure(4)
subplot(2,3,[1 2 3]) 
plot(time,T_c(:,1),time,T_c(:,2),time,T_c(:,3),'LineWidth',1.5)
xlim([0 time(end)])
grid on
grid minor
xlabel('Time [s]')
ylabel('Torque [Nm]')
legend('$T_x$','$T_y$','$T_z$')
subplot(2,3,4) 
plot(time(1:dtm-10),T_c(1:dtm-10,1),time(1:dtm-10),T_c(1:dtm-10,2),time(1:dtm-10),T_c(1:dtm-10,3),'LineWidth',1.5)
xlim([time(1) time(dtm)])
xlabel('Time [s]')
ylabel('Torque [Nm]')
legend('$T_x$','$T_y$','$T_z$')
% title('De-tumbling')
grid on
grid minor
subplot(2,3,5) 
plot(time(dtm:slw),T_c(dtm:slw,1),time(dtm:slw),T_c(dtm:slw,2),time(dtm:slw),T_c(dtm:slw,3),'LineWidth',1.5)
xlim([time(dtm) time(slw)])
xlabel('Time [s]')
ylabel('Torque [Nm]')
legend('$T_x$','$T_y$','$T_z$')
% title('Slew maneuvre')
grid on
grid minor
subplot(2,3,6) 
plot(time(slw:end),T_c(slw:end,1),time(slw:end),T_c(slw:end,2),time(slw:end),T_c(slw:end,3),'LineWidth',1)
xlim([time(slw) time(end)])
xlabel('Time [s]')
ylabel('Torque[Nm]')
legend('$T_x$','$T_y$','$T_z$')
% title('Pointing')
grid on


%% Functions

function [r,v]=kep2car(k,mu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% kep2car:  function useful to convert from Keplerian elements to Cartesian
%           coordinates
%
% INPUTS:
%
%  a         [1x1] Semi-majopr axis        [km]
%  e         [1x1] Eccentricity            [-]
%  i         [1x1] Inclination             [rad]
%  OM        [1x1] RAAN                    [rad]
%  om        [1x1] Pericentre anomaly      [rad]
%  theta     [1x1] True anomaly            [rad]
%  mu        [1x1] Gravitational parameter [km^3/s^2]
%
% OUTPUTS:
%
%  r         [3x1] Position vector         [km]
%  v         [3x1] Velocity vector         [km/s]
%
% First Edition: 01/07/2022
%
% Authors: Rocco Larocca, Matteo Mascelloni, Davide Lanza, Afaq Shakeel

a=k(1);
e=k(2);
i=k(3);
OM=k(4);
om=k(5);
theta=k(6);


p=a*(1-(e^2));
r_norm=p/(1+(e*cos(theta)));
r_PF=r_norm*[cos(theta); sin(theta); 0];
v_PF=(sqrt(mu/p))*[-sin(theta); (e+cos(theta)); 0];
R_3_OM=[cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];
R_1_i=[1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
R_3_om=[cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];
T_PF=(R_3_OM)'*(R_1_i)'*(R_3_om)';
r=T_PF*r_PF;
v=T_PF*v_PF;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = astroConstants(in)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% astroConstants.m - Returns astrodynamic-related physical constants.
%
% PROTOTYPE:
%   out = astro_constants(in)
%
% DESCRIPTION:
%   Returns a row vector of constants, in which there is the corresponding
%   constant for each element of the input vector.
%
%   List of identifiers:
%       Generic astronomical constants:
%           1   Universal gravity constant (G) (from DITAN and Horizon) [km^3/(kg*s^2)]
%           2   Astronomical Unit (AU) (from DE405) [km]
%               Note:  The value for 1 au is from the IAU 2012 Resolution B1.
%       Sun related:
%           3   Sun mean radius (from DITAN) [km]
%           4   Sun planetary constant (mu = mass * G) (from DE405)
%               [km^3/s^2]
%           31  Energy flux density of the Sun (from Wertz,SMAD)
%               [W/m2 at 1 AU]
%       Other:
%           5   Speed of light in the vacuum (definition in the SI and Horizon) [km/s]
%           6   Standard free fall (the acceleration due to gravity on the
%               Earth's surface at sea level) (from Wertz,SMAD) [m/s^2]
%           7   Mean distance Earth-Moon (from Wertz,SMAD) [km]
%           8   Obliquity (angle) of the ecliptic at Epoch 2000 (from
%               Horizon) [rad]
%           9   Gravitatonal field constant of the Earth (from Wertz,SMAD,
%               taken from JGM-2). This should be used in conjunction to
%               Earth radius = 6378.1363 km
%           32  Days in a Julian year y = 365.25 d  (from Horizon)
%       Planetary constants of the planets (mu = mass * G) [km^3/s^2]:
%           11  Me      (from DE405)
%           12  V       (from DE405)
%           13  E       (from DE405)
%           14  Ma      (from DE405)
%           15  J       (from DE405)
%           16  S       (from DE405)
%           17  U       (from DE405)
%           18  N       (from DE405)
%           19  P       (from DE405)
%           20  Moon    (from DE405)
%       Mean radius of the planets [km]:
%           21  Me      (from Horizon)
%           22  V       (from Horizon)
%           23  E       (from Horizon)
%           24  Ma      (from Horizon)
%           25  J       (from Horizon)
%           26  S       (from Horizon)
%           27  U       (from Horizon)
%           28  N       (from Horizon)
%           29  P       (from Horizon)
%           30  Moon    (from Horizon)
%
%   Notes for upgrading this function:
%       It is possible to add new constants.
%       - DO NOT change the structure of the function, as well as its
%           prototype.
%       - DO NOT change the identifiers of the constants that have already
%           been defined in this function. If you want to add a new
%           constant, use an unused identifier.
%       - DO NOT add constants that can be easily computed starting form
%           other ones (avoid redundancy).
%       Contact the author for modifications.
%
% INPUT:
%   in      Vector of identifiers of required constants.
%
% OUTPUT:
%   out     Vector of constants.
%
% EXAMPLE:
%   astroConstants([2, 4, 26])
%      Returns a row vector in which there is the value of the AU, the Sun
%      planetary constant and the mean radius of Saturn.
%
%   astroConstants(10 + [1:9])
%      Returns a row vector with the planetary constant of each planet.
%
% REFERENCES:
%   - DITAN (Direct Interplanetary Trajectory Analysis), Massimiliano
%       Vasile, 2006.
%	- Wertz J. R., Larson W. J., "Space Mission Analysis and Design", Third
%       Edition, Space Technology Library 2003.
%   [A]   DE405 - http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf
%   [B]   Explanatory Supplement to the Astronomical Almanac. 1992. K. P.
%         Seidelmann, Ed., p.706 (Table 15.8) and p.316 (Table 5.8.1),
%         University Science Books, Mill Valley, California. 
%   [C]   Tholen, D.J. and Buie, M.W. 1990. "Further Analysis of
%         Pluto-Charon Mutual Event Observations" BAAS 22(3):1129.
%   [D]   Seidelmann, P.K. et al. 2007. "Report of the IAU/IAG Working
%         Group on cartographic coordinates and rotational elements: 2006"
%         Celestial Mech. Dyn. Astr. 98:155-180. 
%   [F]   Anderson, J.D., et al. 1987. "The mass, gravity field, and
%         ephemeris of Mercury" Icarus 71:337-349.
%   [G]   Konopliv, A.S., et al. 1999. "Venus gravity: 180th degree and
%         order model" Icarus 139:3-18.
%   [H]   Folkner, W.M. and Williams, J.G. 2008. "Mass parameters and
%         uncertainties in planetary ephemeris DE421." Interoffice Memo.
%         343R-08-004 (internal document), Jet Propulsion Laboratory,
%         Pasadena, CA. 
%   [I]   Jacobson, R.A. 2008. "Ephemerides of the Martian Satellites -
%         MAR080" Interoffice Memo. 343R-08-006 (internal document),
%         Jet Propulsion Laboratory, Pasadena, CA. 
%   [J]   Jacobson, R.A. 2005. "Jovian Satellite ephemeris - JUP230"
%         private communication. 
%   [K]   Jacobson, R.A., et al. 2006. "The gravity field of the Saturnian
%         system from satellite observations and spacecraft tracking data"
%         AJ 132(6):2520-2526. 
%   [L]   Jacobson, R.A. 2007. "The gravity field of the Uranian system and
%         the orbits of the Uranian satellites and rings" BAAS 39(3):453. 
%   [M]   Jacobson, R.A. 2008. "The orbits of the Neptunian satellites and
%         the orientation of the pole of Neptune" BAAS 40(2):296. 
%   [N]   Jacobson, R.A. 2007. "The orbits of the satellites of Pluto -
%         Ephemeris PLU017" private communication.
%   [W1]  http://ssd.jpl.nasa.gov/?planet_phys_par Last retrieved
%         20/03/2013
%   [W2]  http://ssd.jpl.nasa.gov/?sat_phys_par Last retrieved
%         20/03/2013
%   [W3]  http://ssd.jpl.nasa.gov/horizons.cgi Last retrieved
%         20/03/2013
%   [M1]  Bills, B.G. and Ferrari, A.J. 1977. ``A Harmonic Analysis of
%         Lunar Topography'', Icarus 31, 244-259.
%   [M2]  Standish, E. M. 1998. JPL Planetary and Lunar Ephemerides,
%         DE405/LE405.
%   [M3]  Lunar Constants and Models Document, Ralph B. Roncoli, 23 Sept 2005,
%         JPL Technical Document D-32296 
%
%
% CALLED FUNCTIONS:
%   (none)
%
% AUTHOR:
%   Matteo Ceriotti, 2006, MATLAB, astroConstants.m
%
% PREVIOUS VERSION:
%   Matteo Ceriotti, 2006, MATLAB, astro_constants.m, Ver. 1.2
%       - Header and function name in accordance with guidlines.
%
% CHANGELOG:
%   26/10/2006, Camilla Colombo: Updated.
%   22/10/2007, Camilla Colombo: astroConstants(8) added (Obliquity (angle)
%       of the ecliptic at Epoch 2000).
%   02/10/2009, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%   12/11/2010, Camilla Colombo: astroConstants(9) added (J2) Note: the
%       present value of J2 is not consistent with the value of the Earth
%       radius. This value of J2 should be used in conjunction to Earth
%       radius = 6378.1363 km
%   19/03/2013, Camilla Colombo: constants updated to NASA JPL website.
%       References added.
%   20/03/2013, REVISION, Francesca Letizia.
%   22/03/2013, Francesca Letizia: all GM from DE405.
%
% -------------------------------------------------------------------------

% 9: J2
% 32: 365.25

out = zeros(1,length(in));
for i=1:length(in)
    switch in(i)
        case 1
            out(i)=6.67259e-20; % From DITAN and Horizon
        case 2
            out(i)=149597870.691; % From DE405
        case 3
            % out(i)=700000; % From DITAN
            out(i)=6.955*10^5; % From Horizon [W3]
        case 4
            % out(i)=0.19891000000000E+31*6.67259e-20; % From DITAN
            out(i)=1.32712440017987E+11; % From DE405 [A]
        case 5
            out(i)=299792.458; % Definition in the SI, Horizon, DE405
        case 6
            out(i)=9.80665; % Definition in Wertz, SMAD
        case 7
            % out(i)=384401; % Definition in Wertz, SMAD
            out(i)=384400; % From Horizon [W3]
        case 8
            % out(i)=23.43928111*pi/180; % Definition in Wertz, SMAD
            out(i)=84381.412/3600*pi/180; % Definition in Horizon
            % obliquity of ecliptic (J2000)    epsilon = 84381.412 (± 0.005) arcsec 
        case 9
            out(i)=0.1082626925638815e-2; % Definition in Wertz, SMAD
        case 11
            % out(i)=0.33020000000000E+24*6.67259e-20; % From DITAN
            %out(i)=0.330104E+24*6.67259e-20;    % From Horizon [F]
            out(i)=2.203208E+4;    % From DE405
        case 12
            % out(i)=0.48685000000000E+25*6.67259e-20; % From DITAN
            %out(i)=4.86732E+24*6.67259e-20;     % From Horizon [G]
            out(i)=3.24858599E+5; % From DE405
        case 13
            % out(i)=0.59736990612667E+25*6.67259e-20; % From DITAN
            % out(i)=5.97219E+24*6.67259e-20;     % From Horizon [H]
            out(i) = 3.98600433e+5; % From DE405
        case 14
            % out(i)=0.64184999247389E+24*6.67259e-20; % From DITAN
            %out(i)=0.641693E+24*6.67259e-20; 	% From Horizon [I]
            out(i) = 4.2828314E+4; %Frome DE405
        case 15
            % out(i)=0.18986000000000E+28*6.67259e-20; % From DITAN
            %out(i)=1898.13E+24*6.67259e-20; 	% From Horizon [J]
            out(i) = 1.26712767863E+08; % From DE405
        case 16
            % out(i)=0.56846000000000E+27*6.67259e-20; % From DITAN
            % out(i)=568.319E+24*6.67259e-20;     % From Horizon [k]
            out(i) = 3.79406260630E+07; % From DE405
        case 17
            % out(i)=0.86832000000000E+26*6.67259e-20; % From DITAN
            % out(i)=86.8103E+24*6.67259e-20;     % From Horizon [L]
            out(i)= 5.79454900700E+06; % From DE405
        case 18
            % out(i)=0.10243000000000E+27*6.67259e-20; % From DITAN
            % out(i)=102.410E+24*6.67259e-20;     % From Horizon [M]
            out(i) = 6.83653406400E+06; % From DE405
        case 19
            % out(i)=0.14120000000000E+23*6.67259e-20; % From DITAN
            %out(i)=.01309E+24*6.67259e-20;     % From Horizon [N]
            out(i) = 9.81601000000E+02; %From DE405
        case 20
            % out(i)=0.73476418263373E+23*6.67259e-20; % From DITAN
             out(i)=4902.801;                 % From Horizon  [M2]
            %out(i)=4902.801076;                % From Horizon  [M3]
        case 21
            % out(i)=0.24400000000000E+04; % From DITAN
            out(i)=2439.7; % From Horizon [D]
        case 22
            % out(i)=0.60518000000000E+04; % From DITAN
            out(i)=6051.8; % From Horizon [D]
        case 23
            % out(i)=0.63781600000000E+04; % From DITAN
            % out(i)=6371.00; % From Horizon [B]
            out(i)=6371.01; % From Horizon [W3]
        case 24
            % out(i)=0.33899200000000E+04; % From DITAN
            % out(i)=3389.50; % From Horizon [D]
            out(i)=3389.9; % From Horizon [W3]            
        case 25
            % out(i)=0.69911000000000E+05; % From DITAN
            out(i)=69911;   % From Horizon [D]
        case 26
            % out(i)=0.58232000000000E+05; % From DITAN
            out(i)=58232;   % From Horizon [D]
        case 27
            % out(i)=0.25362000000000E+05; % From DITAN
            out(i)=25362;   % From Horizon [D]
        case 28
            % out(i)=0.24624000000000E+05; % From DITAN
            % out(i)=24622;   % From Horizon [D]
            out(i)= 24624; % From Horizon [W3]            
        case 29
            % out(i)=0.11510000000000E+04; % From DITAN
            out(i)=1151; 	% From Horizon [C]
        case 30
            % out(i)=0.17380000000000E+04; % From DITAN
            % out(i)=1737.5;  % From Horizon [M1]
            out(i)=1738.0;    % From Horizon  [M3]
        case 31
            out(i)=1367; % From Wertz, SMAD
            % out(i)=1367.6;  % From Horizon  [W3]
        case 32
            out(i)=365.25; % From Horizon
        % Add an identifier and constant here. Prototype:
        % case $identifier$
        %     out(i)=$constant_value$;
        otherwise
            warning('Constant identifier %d is not defined!',in(i));
            out(i)=0;
    end
end

end