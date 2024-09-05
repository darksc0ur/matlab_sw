function [my_k, my_l] = calculate_wavenums(Nx, Ny, Lx, Ly, USE_GPU)
    dk=pi/Lx;
    dl=pi/Ly;

    if USE_GPU == 0
        ksvec(1) = 0; ksvec(Nx/2+1) = 0;
        lsvec(1) = 0; lsvec(Ny/2+1) = 0;
    else
        ksvec(1) = gpuArray(0); ksvec(Nx/2+1) = 0;
        lsvec(1) = gpuArray(0); lsvec(Ny/2+1) = 0;
    end


    for ii=2:(Nx/2)
       ksvec(ii)=ii-1;
       ksvec(Nx/2+ii)=-Nx/2 + ii -1;   
    end
    for ii=2:(Ny/2)
       lsvec(ii)=ii-1;
       lsvec(Ny/2+ii)=-Ny/2 + ii -1;   
    end
    ksvec=ksvec*dk; lsvec=lsvec*dl;
    [my_k, my_l]=meshgrid(ksvec,lsvec);
end
