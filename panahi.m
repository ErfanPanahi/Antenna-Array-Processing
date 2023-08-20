 %%% HW4 - BSS - Erfan Panahi 810198369
clc
clear
fprintf("HW4 - BSS - Erfan Panahi 810198369\n");

%% Definitions:

fs = 1e6;
T = 1e-3;
t = 0:1/fs:T;
f1 = 20e3;
s1 = exp(1i*2*pi*f1*t);
f2 = 10e3;
s2 = exp(1i*2*pi*f2*t);
M = 10;
fc = 150e6;
d = 1;
D = (1:M) - d; % Distance between antennas and first antenna 
c = 3e8;
k = 2*pi*fc/c;
theta1 = 10;
theta2 = 20;
a_theta1 = exp(-1i*k*D'*sind(theta1));
a_theta2 = exp(-1i*k*D'*sind(theta2));
Noise = randn(M,length(s1));
Y = a_theta1*s1 + a_theta2*s2 + Noise;
[U,G,V] = svd(Y);


%% Part b

Theta = 0:0.5:90;
aTheta = exp(-1i*k*D'*sind(Theta));
aUsig = aTheta'*U(:,1:2);
f_beamforming = sqrt(sum(aUsig.*conj(aUsig),2));
figure(1)
plot(Theta,f_beamforming);
title('max[$f(\theta)$] Beamforming Method','Interpreter','latex');
ylabel('$f(\theta)$','Interpreter','latex');
xlabel('$degree(grad)$','Interpreter','latex');
[p,loc] = findpeaks(f_beamforming,'MinPeakHeight',1);
theta1_beam = Theta(loc(1));
theta2_beam = Theta(loc(2));
fprintf("\nPart b:\n Beamforming Method: %s1 = %f , %s2 = %f\n",char(952),theta1_beam,char(952),theta2_beam);

%% Part c

aUnoise = aTheta'*U(:,3:end);
f_MUSIC = 1./sqrt(sum(aUnoise.*conj(aUnoise),2));
figure(2)
plot(Theta,f_MUSIC);
title('max[$f(\theta)$] MUSIC Method','Interpreter','latex');
ylabel('$f(\theta)$','Interpreter','latex');
xlabel('$degree(grad)$','Interpreter','latex');
[p,loc] = findpeaks(f_MUSIC,'MinPeakHeight',1);
theta1_music = Theta(loc(1));
theta2_music = Theta(loc(2));
fprintf("\nPart c:\n MUSIC Method: %s1 = %f , %s2 = %f\n",char(952),theta1_music,char(952),theta2_music);

%% Part d - method1

F = 0:10:5e4;
SF = exp(1i*2*pi*F'*t);
SVsig = SF*V(:,1:2);
g_beamforming_f = sqrt(sum(SVsig.*conj(SVsig),2));
figure(3)
plot(F,g_beamforming_f);
title('max[$g(f)$] Beamforming Method','Interpreter','latex');
ylabel('$g(f)$','Interpreter','latex');
xlabel('$f$','Interpreter','latex');
[p,loc] = findpeaks(g_beamforming_f,'MinPeakHeight',15);
f2_beam = F(loc(1));
f1_beam = F(loc(2));
fprintf("\nPart d:\n Beamforming Method: f1 = %f , f2 = %f\n",f1_beam,f2_beam);

%% Part e - method1

SVnoise = SF*V(:,3:end);
g_music_f = 1./sqrt(sum(SVnoise.*conj(SVnoise),2));
figure(4)
plot(F,g_music_f);
title('max[$g(f)$] MUSIC Method','Interpreter','latex');
ylabel('$g(f)$','Interpreter','latex');
xlabel('$f$','Interpreter','latex');
[p,loc] = findpeaks(g_music_f,'MinPeakHeight',0.05);
f2_music = F(loc(1));
f1_music = F(loc(2));
fprintf("\nPart e:\n MUSIC Method: f1 = %f , f2 = %f\n",f1_music,f2_music);

%% Part d & e - method 2

%beamforming

a1b = exp(-1i*k*D'*sind(theta1_beam));
a2b = exp(-1i*k*D'*sind(theta2_beam));
A_beam = [a1b,a2b]; 
PinvA_beam = pinv(A_beam);
S_hat_beam = PinvA_beam*Y;
f = fs*(-1/2:T:1/2);
S1_f_beam = abs(fftshift(fft(S_hat_beam(1,:))));
S2_f_beam = abs(fftshift(fft(S_hat_beam(2,:))));
figure(5)
plot(f,S1_f_beam,f,S2_f_beam);
legend('fft(s1)','fft(s2)');
title('fft[s(t)] - Beamforming Method','Interpreter','latex');
ylabel('$|S(f)|$','Interpreter','latex');
xlabel('$frequency(Hz)$','Interpreter','latex');
[p1,is1b] = max(S1_f_beam);
[p2,is2b] = max(S2_f_beam);
fprintf("\nPart d-2:\n Beamforming Method (fft): f1 = %f , f2 = %f \n",f(is1b),f(is2b));

%music 

a1m = exp(-1i*k*D'*sind(theta1_music));
a2m = exp(-1i*k*D'*sind(theta2_music));
A_music = [a1m,a2m]; 
PinvA_music = pinv(A_music);
S_hat_music = PinvA_music*Y;
S1_f_music = abs(fftshift(fft(S_hat_music(1,:))));
S2_f_music = abs(fftshift(fft(S_hat_music(2,:))));
figure(6)
plot(f,S1_f_music,f,S2_f_music);
legend('fft(s1)','fft(s2)');
title('fft[s(t)] - MUSIC Method','Interpreter','latex');
ylabel('$|S(f)|$','Interpreter','latex');
xlabel('$frequency(Hz)$','Interpreter','latex');
[p1_,is1m] = max(S1_f_music);
[p2_,is2m] = max(S2_f_music);
fprintf("\nPart e-2:\n MUSIC Method (fft): f1 = %f , f2 = %f \n",f(is1m),f(is2m));






