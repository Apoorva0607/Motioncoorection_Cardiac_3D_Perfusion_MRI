% --------------------------------------------------------------------------------
% F=jFFT(f)
% Jason Mendes - UCAIR - University of Utah - 2012
% --------------------------------------------------------------------------------
% Symmetric FFT.
% --------------------------------------------------------------------------------
function F=jFFT(f)

if (ndims(f)>3)
   error('jFFT(): Data must be 1D, 2D or 3D');
elseif (ndims(f)==3) % 3D
   F=fftshift(fft(fft(fft(ifftshift(f),[],1),[],2),[],3))/sqrt(numel(f));
elseif (min(size(f))>1) % 2D
   F=fftshift(fft2(ifftshift(f)))/sqrt(numel(f));
else % 1D
   F=fftshift(fft(ifftshift(f)))/sqrt(numel(f));
end
