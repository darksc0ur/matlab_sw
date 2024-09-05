% create the grid and truncate
x1 = linspace(-1,1,Nx+1)*Lx;  x = x1(1:end-1);
y1 = linspace(-1,1,Ny+1)*Ly;  y = y1(1:end-1);
[xx,yy] = meshgrid(x,y);
% grid spacing and CFL condition
dx = abs(xx(2,2)-xx(1,1)); dy = abs(yy(2,2)-yy(1,1));

cfl = min([dx dy])/c0;
dt = min([0.05*cfl h]);
if dt>(0.5*cfl)
    disp('dt is getting too close to CFL condition')
end

[k, l] = calculate_wavenums(Nx, Ny, Lx, Ly, USE_GPU);
k2 = k.*k; l2 = l.*l; kl = k.*l;

Dx = sqrt(-1)*k; Dy = sqrt(-1)*l;

myfilt = sbfilter(k, l, Nx, Ny, f_cutoff, f_order, f_strength, USE_GPU);

dump1 = H*H/6;
a11 = 1+dump1*k2;
a12 = (dump1*kl);
a21 = (dump1*kl);
a22 = (1+dump1*l2);
mydet = a11.*a22-a12.*a21;
b11 = a22./mydet;
b22 = a11./mydet;
b12 = -a12./mydet;
b21 = -a21./mydet;
