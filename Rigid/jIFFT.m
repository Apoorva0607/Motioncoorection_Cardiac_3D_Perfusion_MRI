% --------------------------------------------------------------------------------
% f=jIFFT(F)
% Jason Mendes - UCAIR - University of Utah - 2012
% --------------------------------------------------------------------------------
% Symmetric inverse FFT.
% --------------------------------------------------------------------------------
function f=jIFFT(F)

if (ndims(F)>3)
   error('jFFT(): Data must be 1D, 2D or 3D');
elseif (ndims(F)==3) % 3D
   f=fftshift(ifft(ifft(ifft(ifftshift(F),[],1),[],2),[],3))*sqrt(numel(F));
elseif (min(size(F))>1) % 2D
   f=fftshift(ifft2(ifftshift(F)))*sqrt(numel(F));
else % 1D
   f=fftshift(ifft(ifftshift(F)))*sqrt(numel(F));
end
