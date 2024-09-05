function[etap, vortp, kep, tracp] = plotfuncs()
    etap = @etaplot;
    vortp = @vortplot;
    kep = @keplot;
    tracp = @tracplot;
end



function etapp = etaplot(xx, yy, ee1, H, t, labper, Lx, Ly)
     %uu=uu1; vv=vv1; % could use this to get subdomain for extended cases
     %figure(1)
     %clf
     etapp = pcolor(xx,yy,ee1/H);
     shading interp;
     title(['eta/H at t = ' num2str(t/labper,3) ' rotation periods'],'fontweight','bold','fontsize',12)
     axis([-Lx Lx -Ly Ly])
     ylabel('y (m)')
     clim([-1 1]*0.1)
end

function vortpp = vortplot(xx, yy, vmf, umf, Dx, Dy, Lx, Ly)
     vort = real(ifft2(Dx.*vmf - Dy.*umf));
     vortpp = pcolor(xx,yy,vort);
     shading interp;
     clim([-1 1]*0.05)
     axis([-Lx Lx -Ly Ly])
     title('Vorticity')
end

function kepp = keplot(xx, yy, uu1, vv1, t, labper, Lx, Ly)
      ke=0.5*(uu1.^2+vv1.^2);kemax=max(ke(:));
      kepp = pcolor(xx,yy,ke/kemax);
      shading interp;
      title(['scaled KE at t = ' num2str(t/labper,3) ' rotation periods'],'fontweight','bold','fontsize',12)
      axis([-Lx Lx -Ly Ly])
      ylabel('y (m)')
      xlabel('x (m)')
end

function tracpp = tracplot(xx, yy, tdata, Lx, Ly)
     tracpp = pcolor(xx,yy,tdata);
     shading interp;
     %title(['scaled KE at t = ' num2str(t/labper,3) ' rotation periods'],'fontweight','bold','fontsize',12)
     axis([-Lx Lx -Ly Ly])
     clim([0 1])
     xlabel('x (m)')
     title('Tracer')
end