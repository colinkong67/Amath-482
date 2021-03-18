clear all;close all;clc;

load('cam1_1.mat');
load('cam2_1.mat');
load('cam3_1.mat');
[height1_1, width1_1, rgb1_1, numframes1_1] = size(vidFrames1_1);
[height2_1, width2_1, rgb2_1, numframes2_1] = size(vidFrames2_1);
[height3_1, width3_1, rgb3_1, numframes3_1] = size(vidFrames3_1);
%for j=1:num_frames
%    X=vidFrames1_1(:,:,:,j);
%    imshow(X); drawnow
%end
%filter = zeros(480,640);
%filter() = 1;
%% Test 1
X1_1 = [];Y1_1 = [];
X2_1 = [];Y2_1 = [];
X3_1 = [];Y3_1 = [];
for j = 1:numframes1_1
    X = vidFrames1_1(:,:,:,j); 
    X_gray = double(rgb2gray(X));
    X_gray(:,1:320) = 0;
    X_gray(:,380:end) = 0;
    X_gray(1:200,:) = 0;
   % X_modified = X_gray .* filter;
    [M,I] = max(X_gray(:));
    [x1_1,y1_1] = ind2sub(size(X_gray),I);
    X1_1 = [X1_1, x1_1];
    Y1_1 = [Y1_1, y1_1];
end
for j = 1:numframes2_1
    X = vidFrames2_1(:,:,:,j); 
    X_gray = double(rgb2gray(X));
    X_gray(:,1:260) = 0;
    X_gray(:,320:end) = 0;

    [M,I] = max(X_gray(:));
    [x2_1,y2_1] = ind2sub(size(X_gray),I);
    X2_1 = [X2_1, x2_1];
    Y2_1 = [Y2_1, y2_1];
end
for j = 1:numframes3_1
    X = vidFrames3_1(:,:,:,j); 
    X_gray = double(rgb2gray(X));
    X_gray(:,1:250) = 0;
    X_gray(310:end,:) = 0;
    X_gray(:,1:260) = 0;
    [M,I] = max(X_gray(:));
    [x3_1,y3_1] = ind2sub(size(X_gray),I);
    X3_1 = [X3_1, x3_1];
    Y3_1 = [Y3_1, y3_1];
end
[min, I] = min(X1_1(1:50));
X1_1 = X1_1(I:I+200);Y1_1 = Y1_1(I:I+200);
X2_1 = X2_1(I:I+200);Y2_1 = Y2_1(I:I+200);
X3_1 = X3_1(I:I+200);Y3_1 = Y3_1(I:I+200);
X1_1fft = fft(X1_1);Y1_1fft = fft(Y1_1);
X2_1fft = fft(X2_1);Y2_1fft = fft(Y2_1);
X3_1fft = fft(X3_1);Y3_1fft = fft(Y3_1);
X1_1fft(10:end-10) = 0;Y1_1fft(10:end-10) = 0;
X2_1fft(10:end-10) = 0;Y2_1fft(10:end-10) = 0;
X3_1fft(10:end-10) = 0;Y3_1fft(10:end-10) = 0;
X1_1ifft = abs(ifft(X1_1fft));Y1_1ifft = abs(ifft(Y1_1fft));
X2_1ifft = abs(ifft(X2_1fft));Y2_1ifft = abs(ifft(Y2_1fft));
X3_1ifft = abs(ifft(X3_1fft));Y3_1ifft = abs(ifft(Y3_1fft));
figure(1)
subplot(3,2,1)
plot(X1_1,'b','Linewidth',1.5); hold on;
plot(X1_1ifft,'r','Linewidth',1);
axis([0 200 0 480]);
xlabel("Time(Frames)");xlabel("Position of X");
title("Camera 1")
subplot(3,2,2)
plot(Y1_1,'b','Linewidth',1.5); hold on;
plot(Y1_1ifft,'r','Linewidth',1);
axis([0 200 0 640]);
xlabel("Time(Frames)");xlabel("Position of Y");
title("Camera 1")
subplot(3,2,3)
plot(X2_1,'b','Linewidth',1.5); hold on;
plot(X2_1ifft,'r','Linewidth',1);
axis([0 200 0 480]);
xlabel("Time(Frames)");xlabel("Position of X");
title("Camera 2")
subplot(3,2,4)
plot(Y2_1,'b','Linewidth',1.5); hold on;
plot(Y2_1ifft,'r','Linewidth',1);
axis([0 200 0 640]);
xlabel("Time(Frames)");xlabel("Position of Y");
title("Camera 2")
subplot(3,2,5)
plot(X3_1,'b','Linewidth',1.5); hold on;
plot(X3_1ifft,'r','Linewidth',1);
axis([0 200 0 480]);
xlabel("Time(Frames)");xlabel("Position of X");
title("Camera 3")
subplot(3,2,6)
plot(Y3_1,'b','Linewidth',1.5); hold on;
plot(Y3_1ifft,'r','Linewidth',1);
axis([0 200 0 640]);
xlabel("Time(Frames)");xlabel("Position of Y");
title("Camera 3")

%% Test 2
load('cam1_2.mat');
load('cam2_2.mat');
load('cam3_2.mat');
[height1_2, width1_2, rgb1_2, numframes1_2] = size(vidFrames1_2);
[height2_2, width2_2, rgb2_2, numframes2_2] = size(vidFrames2_2);
[height3_2, width3_2, rgb3_2, numframes3_2] = size(vidFrames3_2);
X1_2 = [];Y1_2 = [];
X2_2 = [];Y2_2 = [];
X3_2 = [];Y3_2 = [];
for j = 1:numframes1_2
    X = vidFrames1_2(:,:,:,j); 
    X_gray = double(rgb2gray(X));
    X_gray(:,1:450) = 0;
    X_gray(:,500:end) = 0;
    X_gray(1:313,:) = 0;
   % X_modified = X_gray .* filter;
    [M,I] = max(X_gray(:));
    [x1_2,y1_2] = ind2sub(size(X_gray),I);
    X1_2 = [X1_2, x1_2];
    Y1_2 = [Y1_2, y1_2];
end
for j = 1:numframes2_2
    X = vidFrames2_2(:,:,:,j); 
    X_gray = double(rgb2gray(X));
    X_gray(:,1:260) = 0;
    X_gray(:,320:end) = 0;
    [M,I] = max(X_gray(:));
    [x2_2,y2_2] = ind2sub(size(X_gray),I);
    X2_2 = [X2_2, x2_2];
    Y2_2 = [Y2_2, y2_2];
end
for j = 1:numframes3_2
    X = vidFrames3_2(:,:,:,j); 
    X_gray = double(rgb2gray(X));
    X_gray(:,1:250) = 0;
    X_gray(310:end,:) = 0;
    X_gray(:,1:260) = 0;
    [M,I] = max(X_gray(:));
    [x3_2,y3_2] = ind2sub(size(X_gray),I);
    X3_2 = [X3_2, x3_2];
    Y3_2 = [Y3_2, y3_2];
end
[max, I] = max(X1_2(1:50));
X1_2 = X1_2(I:I+200);Y1_2 = Y1_2(I:I+200);
X2_2 = X2_2(I:I+200);Y2_2 = Y2_2(I:I+200);
X3_2 = X3_2(I:I+200);Y3_2 = Y3_2(I:I+200);
Xsvd = [X1_2;Y1_2;X2_2;Y2_2;X3_2;Y3_2];
[U,S,V] = svd(Xsvd);
Xsvd = U(:,1:2)*S(1:2,1:2)*V(:,1:2)';
X1_2 = Xsvd(1,:);Y1_2 = Xsvd(2,:);
X2_2 = Xsvd(3,:);Y2_2 = Xsvd(4,:);
X3_2 = Xsvd(5,:);Y3_2 = Xsvd(6,:);
X1_2fft = fft(X1_2);Y1_2fft = fft(Y1_2);
X2_2fft = fft(X2_2);Y2_2fft = fft(Y2_2);
X3_2fft = fft(X3_2);Y3_2fft = fft(Y3_2);
X1_2fft(10:end-10) = 0;Y1_2fft(10:end-10) = 0;
X2_2fft(10:end-10) = 0;Y2_2fft(10:end-10) = 0;
X3_2fft(10:end-10) = 0;Y3_2fft(10:end-10) = 0;
X1_2ifft = abs(ifft(X1_2fft));Y1_2ifft = abs(ifft(Y1_2fft));
X2_2ifft = abs(ifft(X2_2fft));Y2_2ifft = abs(ifft(Y2_2fft));
X3_2ifft = abs(ifft(X3_2fft));Y3_2ifft = abs(ifft(Y3_2fft));
Xisvd = [X1_2ifft;Y1_2ifft;X2_2ifft;Y2_2ifft;X3_2ifft;Y3_2ifft];
[Ui,Si,Vi] = svd(Xisvd);
Xisvd = Ui(:,1:2)*Si(1:2,1:2)*Vi(:,1:2)';
X1_2ifft = Xisvd(1,:);Y1_2ifft = Xisvd(2,:);
X2_2ifft = Xisvd(3,:);Y2_2ifft = Xisvd(4,:);
X3_2ifft = Xisvd(5,:);Y3_2ifft = Xisvd(6,:);
figure(2)
subplot(3,2,1)
plot(X1_2,'b','Linewidth',1.5); hold on;
plot(X1_2ifft,'r','Linewidth',1);
axis([0 200 0 480]);
xlabel("Time(Frames)");xlabel("Position of X");
title("Camera 1")
subplot(3,2,2)
plot(Y1_2,'b','Linewidth',1.5); hold on;
plot(Y1_2ifft,'r','Linewidth',1);
axis([0 200 0 640]);
xlabel("Time(Frames)");xlabel("Position of Y");
title("Camera 1")
subplot(3,2,3)
plot(X2_2,'b','Linewidth',1.5); hold on;
plot(X2_2ifft,'r','Linewidth',1);
axis([0 200 0 480]);
xlabel("Time(Frames)");xlabel("Position of X");
title("Camera 2")
subplot(3,2,4)
plot(Y2_2,'b','Linewidth',1.5); hold on;
plot(Y2_2ifft,'r','Linewidth',1);
axis([0 200 0 640]);
xlabel("Time(Frames)");xlabel("Position of Y");
title("Camera 2")
subplot(3,2,5)
plot(X3_2,'b','Linewidth',1.5); hold on;
plot(X3_2ifft,'r','Linewidth',1);
axis([0 200 0 480]);
xlabel("Time(Frames)");xlabel("Position of X");
title("Camera 3")
subplot(3,2,6)
plot(Y3_2,'b','Linewidth',1.5); hold on;
plot(Y3_2ifft,'r','Linewidth',1);
axis([0 200 0 640]);
xlabel("Time(Frames)");xlabel("Position of Y");
title("Camera 3")

