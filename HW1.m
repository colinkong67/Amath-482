clear all; close all; clc

load subdata.mat % Imports the data as the 262144x49 (space by time) matrix called subdata
L = 10; % spatial domain
n = 64; % Fourier modes
x2 = linspace(-L,L,n+1); x = x2(1:n); y =x; z = x;
k = (2*pi/(2*L))*[0:(n/2 - 1) -n/2:-1]; ks = fftshift(k);

[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

for j=1:49
    Un(:,:,:)=reshape(subdata(:,j),n,n,n);
    M = max(abs(Un),[],'all');
    close all, isosurface(X,Y,Z,abs(Un)/M,0.7)
    axis([-20 20 -20 20 -20 20]), grid on, drawnow
    pause(1)
end

% Average the spectrum with normalizing fourier transformed data
ave = zeros(1,n);
for j = 1:49
    Un(:,:,:)=reshape(subdata(:,j),n,n,n);
    Utn = fftn(Un);
    ave = ave + Utn;
end
ave = abs(fftshift(ave))/49;
figure(1);
close all, isosurface(Kx,Ky,Kz,ave./max(ave(:)),0.7);
axis([-15 15 -15 15 -15 15]), grid on, drawnow
title("Averaging of spectrum ");
xlabel("Kx"),ylabel("Ky"),zlabel("Kz");

% Center Frequency
[x,y,z] = ind2sub([n,n,n], find(ave == max(abs(ave(:)))));
centFreq = [Kx(x,y,z), Ky(x,y,z), Kz(x,y,z)];
kx_0 = Kx(x,y,z); ky_0 = Ky(x,y,z); kz_0 = Kz(x,y,z);

% Gaussian filter
tau = 0.2;
filter = exp(-tau*((Kx-kx_0).^2+(Ky-ky_0).^2+(Kz-kz_0).^2));

% Determinethe path of the submarine.
positions = zeros(49,3);
for j = 1:49
    Un_f(:,:,:)=reshape(subdata(:,j),n,n,n);
    Utn_f = fftn(Un_f);
    Unft = filter.*fftshift(Utn_f);
    Unf = ifftn(Unft);
    Max = max(abs(Unf(:)));
    [positionsX,positionsY,positionsZ] = ind2sub([n,n,n], find(abs(Unf)==Max));
    positions(j,1) = X(positionsX,positionsY,positionsZ);
    positions(j,2) = Y(positionsX,positionsY,positionsZ);
    positions(j,3) = Z(positionsX,positionsY,positionsZ);
end
figure(2); 
plot3(positions(:,1),positions(:,2),positions(:,3),".-","Linewidth",1.5,"Markersize", 15);
axis([-10 10 -10 10 -10 10]);
xlabel("x"),ylabel("y"),zlabel("z"), grid on;
title("Path of the Submarine");
hold on
% Final positions for P-8 Poseidon subtracking aircraft
plot3(positions(49,1),positions(49,2),positions(49,3),"r*","Markersize",15);






