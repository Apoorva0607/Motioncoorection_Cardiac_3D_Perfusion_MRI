function [I2,x,y,z,E] = Reg_GS_3D(I1,I0,step)
[sx,sy,sz] = size(I0);
Maxiter = 110;

a = 20;
s = 2;

kx = cos(2*pi*(0:sx-1)/sx);
ky = cos(2*pi*(0:sy-1)/sy);
kz = cos(2*pi*(0:sz/(sz-1):sz));
W = 2*(kx+ky.'+permute(kz,[3 1 2])-3);
W = (1-a*W).^-s;
 
[y,x,z] = meshgrid(1:sx,1:sy,1:sz);

for iter=1:Maxiter
    
    I2 = interp3(I1,y,x,z,'cubic');
    I2(isnan(I2))=0;
    figure(1)
    imagesc(I2(:,:,4))
    colormap gray
    axis image
    drawnow
    
    ddx = 0.5*(I2(3:end,:,:) - I2(1:end-2,:,:));
    ddy = 0.5*(I2(:,3:end,:) - I2(:,1:end-2,:));
    ddz = 0.5*(I2(:,:,3:end) - I2(:,:,1:end-2));
    ddx = cat(1,I2(2,:,:) - I2(1,:,:),ddx,I2(end,:,:) - I2(end-1,:,:));
    ddy = cat(2,I2(:,2,:) - I2(:,1,:),ddy,I2(:,end,:) - I2(:,end-1,:));
    ddz = cat(3,I2(:,:,2) - I2(:,:,1),ddz,I2(:,:,end) - I2(:,:,end-1));
    
    dI = I2-I0;
    
    dx = W.*fft3(ddx.*dI);
    dx = -real(ifft3(dx));

    dy = W.*fft3(ddy.*dI);
    dy = -real(ifft3(dy));
    
    dz = W.*fft3(ddz.*dI);
    dz = -real(ifft3(dz));
    
    d = sqrt(dx.^2+dy.^2+dz.^2);
    md = max(d(:));
    
    dx = dx/md;
    dy = dy/md;
    dz = dz/md;

    E(iter) = sum(d(:));
    figure(2)
    plot(E)
    drawnow
    
    x = x + step*dx;
    y = y + step*dy;
    z = z + step*dz;
    
    x(x<1) = 1;
    x(x>sx) = sx;
    y(y<1) = 1;
    y(y>sy) = sy;
    z(z<1) = 1;
    z(z>sz) = sz;
   
end

