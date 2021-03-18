%% Start code
figure(1)
[y, Fs] = audioread('GNR.m4a');
y = y.';
trgnr = length(y)/Fs; % record time in seconds
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Sweet Child O''Mine');
%p8 = audioplayer(y,Fs); playblocking(p8);
%%
n = length(y);
L = trgnr;
t2 = linspace(0,L,n+1); 
t = t2(1:n);
k = (1/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);
a = 500;
tau = 0:0.1:L;
Gnrguitar = [];
for j = 1:length(tau)
    g = exp(-a*(t - tau(j)).^2);
    Sg = g.*y;
    Sgt = fft(Sg);
    Sgt_spec(:,j) = fftshift(abs(Sgt));
    [M,I] = max(Sgt);
    Gnrguitar = [Gnrguitar; abs(k(I))];
end
figure(2) %% spectrogram of the guitar part of GNR.clip
pcolor(tau,ks,abs(Sgt_spec))
shading interp
set(gca,'ylim',[50 1000],'Fontsize',16)
colormap(hot)
xlabel('time (t)-seconds'), ylabel('frequency (k)-Hz')
title('Spectrogram of Sweet Child O'' Mine','Fontsize',16)
figure(3) %% Music scoree sheet
plot(tau, Gnrguitar, 'ko', 'MarkerFaceColor', 'r')
xlabel('Time (sec)'), ylabel('Frequency (Hz)')
yticks([277.18, 311.13, 369.99, 415.30, 554.37, 698.46, 739.99])
yticklabels({'C4#', 'D4#', 'F4#', 'G4#', 'C5#', 'F5', 'F5#'})
set(gca,'Fontsize',12,'ylim',[200 800])
title('Music score for the guitar in GNR clip')
%%
figure(4)
[y, Fs] = audioread('Floyd.m4a');
trfloyd = length(y)/Fs; % record time in seconds
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Comfortably Numb');
L = 15;
n = L*Fs;
y = y(1:n).';
t2 = linspace(0,L,n+1); 
t = t2(1:n);
k = (1/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);
a = 200;
tau = 0:0.1:L;
Floydbass = [];
for j = 1:length(tau)
    g = exp(-a*(t - tau(j)).^2);
    Sg = g.*y;
    Sgt = fft(Sg);
    Sgt_spec2(:,j) = fftshift(abs(Sgt));
    [M,I] = max(Sgt);
    Floydbass = [Floydbass; abs(k(I))];
end
figure(5) %% spectrogram of the bass part of Floyd.clip
pcolor(tau,ks,abs(Sgt_spec2)) 
shading interp
set(gca,'ylim',[60 200],'Fontsize',16)
colormap(hot)
xlabel('time (t)-seconds'), ylabel('frequency (k)-Hz')
title('Spectrogram of Comfortably Numb','Fontsize',16)
figure(6) %% Music scoree sheet
plot(tau, Floydbass, 'ko', 'MarkerFaceColor', 'r')
xlabel('Time (sec)'), ylabel('Frequency (Hz)')
yticks([82.407, 92.499, 97.999, 110.00, 123.47])
yticklabels({'E2', 'F2#', 'G2', 'A2', 'B2'})
set(gca,'Fontsize',12,'ylim',[60 160])
title('Music score for the bass in Floyd clip')
%% 
[y, Fs] = audioread('Floyd.m4a');
trfloyd=length(y)/Fs;
y = y(1:end-1);
y = y.';
n = length(y);
L = trfloyd;
k = (1/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);
yt = fft(y);
[Ny,Nx] = size(yt);
wx = 10000;
wy = 10000;
filter = ones(size(y));
filter(wy+1:Nx-wy+1) = zeros(1,Nx-2*wy+1); %% apply Shannon filter
yf = yt.*filter;
ynf = ifft2(yf);
figure(7)
plot(ks,abs(ynf),'Linewidth',2);

t2 = linspace(0,L,n+1); 
t = t2(1:n);
a = 500;
tau = 0:0.1:L;
for j = 1:length(tau)
    g = exp(-a*(t - tau(j)).^2);
    Sg = g.*ynf;
    Sgt = fft(Sg);
    Sgt_spec3(:,j) = abs(Sgt);
end
figure(8)
pcolor(tau,ks,abs(Sgt_spec3)) 
shading interp
set(gca,'ylim',[50 1000],'Fontsize',16)
colormap(hot)
xlabel('time (t)-seconds'), ylabel('frequency (k)-Hz')
title('Spectrogram of Comfortably Numb','Fontsize',16)




