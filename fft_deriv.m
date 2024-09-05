function ans = fft_deriv(in_data, wavenums)
% This function uses FFT to calculate the derivative.

% The wavenumber that you pass to the function implies the direction
% that you are calculating the derivative in.

% If you want to use parallel computation pass in_data and wavenums as a
% gpuArray.
    i = sqrt(-1);
    ans = ifft2(i*wavenums.*fft2(in_data));
end