% Physical Parameters
g0=9.81*0.01; g=g0;
labper = 24; % lab scale rotation period, try 120
f = 2*2*pi/labper; % lab scale
cd = 0.00025; % drag coefficient
Lx = 5; Ly = 5;
H = 0.2; c0 = sqrt(g*H); 

% Grid Parameters
Nx = 128*4; %zonal resolution
Ny = 128*4; %meridional resolution

% Timestepping Parameters
numsteps = 500;
numouts = 10;
t = 0;
h = 0.01;

% Filter Parameters
f_cutoff = 0.3;
f_order = 4;
f_strength = 1e-1;

% Overhead automatically creates wavenumbers, grid and filter
overhead

