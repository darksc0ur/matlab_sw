function myfilt = sbfilter(k, l, Nx, Ny, f_cutoff, f_order, f_strength, USE_GPU)
    %%Build filter
    kmax = max(k(:));
    lmax = max(l(:));
    %% SUBICH
    kcrit = kmax*f_cutoff;
    lcrit = lmax*f_cutoff;

    if USE_GPU == 0
        myfilt = ones(size(k));
    else
        myfilt = ones(size(k), 'gpuArray');
    end
    %filter in x
    mymask = (abs(k)<kcrit);
    myfilt = myfilt.*(mymask + ... % Less than cutoff, no filtering
               (1-mymask).* ... % exponential filter
               exp(log(f_strength)*((abs(k)-kcrit)./(max(k(:))-kcrit)).^f_order));
    
    %filter in y
    mymask = (abs(l)<lcrit);
    myfilt = myfilt.*(mymask + (1-mymask).* ...
		      exp(log(f_strength).*((abs(l)-lcrit)./(max(l(:))-lcrit)).^f_order));
    
    % remove the Nyquist frequency entirely
    myfilt(:,floor(Nx/2+1)) = 0;
    myfilt(floor(Ny/2+1),:) = 0;
end