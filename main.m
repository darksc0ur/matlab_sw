% SW equations in doubly periodic domain in physical formulation
% not in conservative form including bottom and surface stress
% inertia free particles optional
clear all,close all
%DEVICE = gpuDevice(1);
USE_GPU = 1;
SHOW_PLOTS = 1;

% Load physical parameters
configuration

% Load Initial Conditions
ee0 = 0.4*H*sech((yy-0.1*Lx*sin(4*pi*xx/Lx))/(0.05*Lx)).^8;
uu0 = zeros(size(ee0));
vv0 = zeros(size(ee0));

% a band of tracer (dye)

cc0 = exp(-(xx/(0.1*Lx)).^2);
dd0 = exp(-(yy/(0.1*Ly)).^2);

if USE_GPU == 1
    ee0 = gpuArray(ee0);
    uu0 = gpuArray(uu0);
    vv0 = gpuArray(vv0);
    cc0 = gpuArray(cc0);
    dd0 = gpuArray(dd0);
  
end

%% One Euler Step Before Leapfrog
% Calculate some derivatives
Dux = fft_deriv(uu0, k); % X derivative of u
Duy = fft_deriv(uu0, l); % y derivative of u
Dvx = fft_deriv(vv0, k); % x derivative of v
Dvy = fft_deriv(vv0, l); % y derivative of v
Dex = fft_deriv(ee0, k); % x derivative of eta
Dey = fft_deriv(ee0, l); % y derivative of eta
Dheux = fft_deriv((H + ee0).*uu0, k); % x derivative of (H + eta)*u
Dhevy = fft_deriv((H + ee0).*vv0, l); % y derivative of (H + eta)*v
Dcx = fft_deriv(cc0, k); % x derivative of tracer
Dcy = fft_deriv(cc0, l); % y derivative of tracer
Ddx = fft_deriv(dd0, k);
Ddy = fft_deriv(dd0, l);

% Timestep
uu1 = uu0 - dt*(uu0.*Dux + vv0.*Duy + g*Dex - f*vv0);
vv1 = vv0 - dt*(uu0.*Dvx + vv0.*Dvy + g*Dey + f*uu0);
ee1 = ee0 - dt*(Dheux + Dhevy);
cc1 = cc0 - dt*(uu0.*Dcx + vv0.*Dcy);
dd1 = dd0 - dt*(uu0.*Ddx + vv0.*Ddy);

%% Leapfrog Iterations
tic
for ii=1:numouts
    for jj = 1:numsteps
    t=t+dt;
    funf = fft2(f.*uu1);
    fvnf = fft2(f.*vv1);
    umf = fft2(uu0);
    vmf = fft2(vv0);
    emf = fft2(ee0);

    uux = fft2(uu1.*real(fft_deriv(uu1, k)));
    uvx = fft2(uu1.*real(fft_deriv(vv1, k)));
    vuy = fft2(vv1.*real(fft_deriv(uu1, l)));
    vvy = fft2(vv1.*real(fft_deriv(vv1, l)));
    eex = Dx.*fft2(ee1);
    eey = Dy.*fft2(ee1);

    epc = Dx.*fft2((H+ee1).*uu1) + Dy.*fft2((H+ee1).*vv1);
    ccx = Dx.*fft2(cc1);
    ccy = Dy.*fft2(cc1);
    cpc = uu1.*real(ifft2(ccx))+vv1.*real(ifft2(ccy));
    ddx = Dx.*fft2(dd1);
    ddy = Dy.*fft2(dd1);
    dpc = uu1.*real(ifft2(ddx)) + vv1.*real(ifft2(ddy));

    spd=sqrt(uu1.^2 + vv1.^2);

    % The bottom drag
    taubx=-fft2(cd*spd.*uu1);
    tauby=-fft2(cd*spd.*vv1);
    phx=-g*real(ifft2(eex));phy=-g*real(ifft2(eey));
    rhu=a11.*umf+a12.*vmf-2*dt*(-fvnf+uux+vuy+g*eex-taubx);
    rhv=a21.*umf+a22.*vmf-2*dt*(funf+uvx+vvy+g*eey-tauby);
    uu2=b11.*rhu+b12.*rhv;
    vv2=b12.*rhu+b22.*rhv;
    ee2=emf-2*dt*epc;
    ccpf = fft2(cc0-2*dt*cpc);
    ddpf = fft2(dd0 - 2*dt*dpc);

    % hyperviscosity or filter
     uu2=real(ifft2(uu2.*myfilt));
     vv2=real(ifft2(vv2.*myfilt));
     ee2=real(ifft2(ee2.*myfilt));
     cc2 = real(ifft2(ccpf.*myfilt));
     dd2 = real(ifft2(ddpf.*myfilt));

    % Cycle Arrays
    uu0 = uu1;
    uu1 = uu2;
    vv0 = vv1;
    vv1 = vv2;
    ee0 = ee1;
    ee1 = ee2;
    cc0 = cc1;
    cc1 = cc2;
    dd0 = dd1;
    dd1 = dd2;
    end
    if (any(isnan(uu1(:)) | isnan(vv1(:))))
    disp('NaNs detected');
    break;
    end
if SHOW_PLOTS == 1
     figure(1)
     clf
     % Plotting Defaults
     colormap(darkjet);
     set(gcf,'DefaultLineLineWidth',2,'DefaultTextFontSize',12,...
         'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
         'DefaultAxesFontWeight','bold');
    [etaplot, vortplot, keplot, tracplot] = plotfuncs();
    subplot(2, 2, 1);
    etaplot(xx, yy, ee1, H, t, labper, Lx, Ly);
    subplot(2, 2, 2);
    vortplot(xx, yy, vmf, umf, Dx, Dy, Lx, Ly);
    subplot(2, 2, 3);
    keplot(xx, yy, uu1, vv1, t, labper, Lx, Ly);
    subplot(2, 2, 4);
    tracplot(xx, yy, cc1, Lx, Ly);
    drawnow
 end
end
toc