clear all; close all; clc
%% Load Video 1
v1 = VideoReader("monte_carlo_low.mp4");
vid1_dt = 1/v1.Framerate;
vid1_t = 0:1:v1.Duration;
vid1Frames = read(v1);
[height1, width1, RGB1, numFrames1] = size(vid1Frames);

v2 = VideoReader("ski_drop_low.mp4");
vid2dt = 1/v2.Framerate;
vid2_t = 0:1:v2.Duration;
vid2Frames = read(v2);
[height2, width2, RGB2, numFrames2] = size(vid2Frames);

for i = 1:numFrames1
    X = vid1Frames(:,:,:,i);
%    imshow(X);drawnow;
end
%% Set the size of the frame
numRows = 500-49;
numCols = 600-299;
gray_vid1 = zeros(numRows, numCols, numFrames1);

for j=1:numFrames1
    gimage = rgb2gray(vid1Frames(50:500,300:600,:,j));
    gray_vid1(:,:,j) = abs(255-gimage);
    %imshow(abs(255-gimage));drawnow
end    

X = reshape(gray_vid1, numRows*numCols, numFrames1);
height = numRows;
width = numCols;

X1 = X(:,1:end-1);
X2 = X(:,2:end);
%% Setup DMD
[U,S,V] = svd(X1,'econ');
r = 2;
U_r = U(:, 1:r);
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);
A_tilde = U_r' * X2 * V_r / S_r;
[W_r,D] = eig(A_tilde);
Phi = X2 * V_r / S_r * W_r;

lambda = diag(D);
omega = log(lambda)/vid1_dt;

%thresh = 0.001;
%bg = find(abs(omega) < thresh);
%omega_bg = omega(bg);
%plot(real(omega_bg),imag(omega_bg),'rx');

x1 =X1(:, 1);
b = Phi \ x1;
%% DMD Reconstruction
mm1 = size(X1,2);
time_dynamics = zeros(r,mm1);
t = (0:mm1-1)*vid1_dt;
for iter = 1:mm1
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));    
end
X_bg = Phi * time_dynamics;

for k = 1:numFrames1 - 1
    frame = reshape(X_bg(:,k),height,width);
    frame = uint8(real(frame));
    figure(2);
    %imshow(frame); drawnow
end
%% Set the foreground and the background
X_fg = X1 - abs(X_bg);
ind = find(X_fg < 0);
X_rec = X_fg+X_bg;
X_bgr = X_bg;
X_bgr(ind) = X_bg(ind) + X_fg(ind);
X_fgr = X_fg;
X_fgr(ind) = 0;

for i = 1:numFrames1 - 1
    frame = reshape(X_fg(:,i),height,width); 
    frame = uint8(real(frame));
    figure(3);
    %imshow(frame); drawnow
end
%% Frame out
figure(4)
subplot(2,3,1)
frame = reshape(X1(:,100),height,width); 
frame = uint8(real(frame));
imshow(reshape(frame,[height,width]));
title("Original Video");
subplot(2,3,2)
frame = reshape(X_rec(:,100),height,width); 
frame = uint8(real(frame));
imshow(reshape(frame,[height,width]));
title("Sum with Reconstruction");
subplot(2,3,3)
frame = reshape(X_bg(:,100),height,width); 
frame = uint8(real(frame));
imshow(reshape(frame,[height,width]));
title("Background without r");
subplot(2,3,4)
frame = reshape(X_bgr(:,100),height,width); 
frame = uint8(real(frame));
imshow(reshape(frame,[height,width]));
title("Background with r");
subplot(2,3,5)
frame = reshape(X_fg(:,100),height,width); 
frame = uint8(real(frame));
imshow(reshape(frame,[height,width]));
title("Foreground without r");
subplot(2,3,6)
frame = reshape(X_fgr(:,100),height,width); 
frame = uint8(real(frame));
imshow(reshape(frame,[height,width]));
title("Foreground with r");
