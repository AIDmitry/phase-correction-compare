clc;
close all;
clear all;

Fs = 1.2e6;
df01=-3096.8; % оценка частотной отстройки для каждого из сигналов, данные от СТЦ
df02=-176037;
debug=1;
duration = 5;

fileID = fopen('/home/dmirty/MatlabProjects/GenSig/PSP/ml_psp0_1200000_1000000_5sec.bin', 'r');
psp0 = fread(fileID, 2*duration*Fs, "short");
psp= complex(psp0(1:2:end), psp0(2:2:end));% чтение сигналов из файлов, приведение их в коплексный вид

fileID = fopen('/home/dmirty/MatlabProjects/GenSig/Raw_Eutelsat21B_10990V_2022_10_05_15_33_09_000/Data.bin', 'r');
signal10 = fread(fileID, 2*duration*Fs, "short");
signal1= complex(signal10(1:2:end), signal10(2:2:end));
t= (1:length(signal1))'/Fs;
signal1 = signal1.*(exp(1i*2*pi*(df01)*t));

fileID = fopen('/home/dmirty/MatlabProjects/GenSig/Raw_Eutelsat10A_10990V_2022_10_05_15_33_09_000/Data.bin', 'r');
signal20 = fread(fileID, 2*duration*Fs, "short");
signal2= complex(signal20(1:2:end), signal20(2:2:end));
t= (1:length(signal2))'/Fs;
signal2 = signal2.*(exp(1i*2*pi*(df02)*t));

% figure;
% title(["Корреляция сигнала sat21"]);
% [~, S]= max(abs(xcorr(psp, signal1)));
% % subplot(2,1,1);
% % plot(abs(xcorr(psp, signal1)));
% % subplot(2,1,2);
% % plot(abs(xcorr(psp, psp)));
% [~, S0]= max(abs(xcorr(psp, psp))); % вычисление рассогласование спутниковых сигналов и псп
% lag1=S0-S;
% 
% % figure;
% % title(["Корреляция сигнала sat10"]);
% [~, S]= max(abs(xcorr(psp, signal2)));
% % subplot(2,1,1);
% % plot(abs(xcorr(psp, signal2)));
% % subplot(2,1,2);
% % plot(abs(xcorr(psp, psp)));
% [~, S0]= max(abs(xcorr(psp, psp)));
% lag2=S0-S;

lag1=304000;
lag2=306357;

signal2=signal2(1+lag2:end); %компенсация рассогласований
signal1=signal1(1+lag1:end);

%% without correction sat21
time=[];
for i=1:7
psp1=psp(1:i*Fs/2);

sig1=signal1(1:i*Fs/2);

Corr(i)=max(abs(xcorr(psp1, sig1)))/mean(abs(xcorr(psp1, sig1)));
time(i)=i/2;
clear sig1;
clear psp1;
end

figure;
subplot(2, 1, 1);
plot(time, Corr);
title(["Корреляция от времени сигнала sat21 на psp до и после побработки"]);
clear sign1;
clear psp1;

%% Correction sat21

[ErrorSig1, ErrorSig1_, df01, Fs1] = PhaseFluctuationsEstim2(signal1,psp(1:end-lag1),Fs); % оценка искажений

sig_est1=ErrorSig1_;
s1 = signal1;

t= (1:length(s1))'/Fs;
s1= s1.*exp(1i*2*pi*(-df01)*t);

L = Fs/Fs1;
sig_est_interp1 = repmat(sig_est1,1, L); % интерполяция сигнала ошибки на ЧД сигналов
sig_est_interp1 = sig_est_interp1.';
sig_est_interp1 = sig_est_interp1(:);



L = min(length(sig_est_interp1), length(s1));

sig_est_interp1 = sig_est_interp1(1:L);
s1=s1(1:L);

s1=s1.*conj(sig_est_interp1); % применеине коррекции

%% with correction sat21

L = min(L, length(psp));

for i=1:7

psp1=psp(1:i*Fs/2);

sig1=s1(1:i*Fs/2);

C=abs(xcorr(psp1, sig1));
%plot(C);
Max = max(C);
M=mean(C);
Corr1(i) = Max/M;

time(i)=i/2;

end

subplot(2, 1, 2);
plot(time, Corr1);

%% whithout correction sat10
time=[];
for i=1:7
psp1=psp(1:i*Fs/2);

sig2=signal2(1:i*Fs/2);

Corr(i)=max(abs(xcorr(psp1, sig2)))/mean(abs(xcorr(psp1, sig2)));
time(i)=i/2;
clear sig2;
clear psp1;
end

figure;
subplot(2, 1, 1);
plot(time, Corr);
title(["Корреляция от времени сигнала sat10 на psp до и после обработки"]);
clear sig2;

%% Correction sat10

[ErrorSig2, ErrorSig2_, df02, Fs2] = PhaseFluctuationsEstim2(signal2,psp(1:end-lag2),Fs);

sig_est2=ErrorSig2_;

s2 = signal2;

t= (1:length(s2))'/Fs;
s2= s2.*exp(1i*2*pi*(-df02)*t);

L = Fs/Fs2;
sig_est_interp2 = repmat(sig_est2,1, L);
sig_est_interp2 = sig_est_interp2.';
sig_est_interp2 = sig_est_interp2(:);


L = min(length(sig_est_interp2), length(s2));

sig_est_interp2 = sig_est_interp2(1:L);
s2=s2(1:L);

s2=s2.*conj(sig_est_interp2);

%% with correction sat10

L = min(L, length(psp));

for i=1:7

psp1=psp(1:i*Fs/2);

sig2=s2(1:i*Fs/2);
C=xcorr(psp1, sig2);
Corr2(i)=max(abs(C))/mean(abs(C));
time(i)=i/2;

end

subplot(2, 1, 2);
plot(time, Corr2);
%% crosscorr

% figure;
% subplot(2, 1, 1);
% plot(abs(xcorr(signal1, signal2)));
% subplot(2,1,2);fileID
% plot(abs(xcorr(s1, s2)));
% 
df0=250.32;
t= (1:length(signal1))'/Fs;
signal1 = signal1.*(exp(1i*2*pi*df0*t));

t= (1:length(s1))'/Fs;
s1 = s1.*(exp(1i*2*pi*df0*t));

for i= 1:7
sig1befor=signal1(1:i*Fs/2);
sig2befor=signal2(1:i*Fs/2);
sig1after=s1(1:i*Fs/2);
sig2after=s2(1:i*Fs/2);
crosscoerr1(i)=max(abs(xcorr(sig1befor, sig2befor)))/mean(abs(xcorr(sig1befor, sig2befor)));
crosscoerr2(i)=max(abs(xcorr(sig1after, sig2after)))/mean(abs(xcorr(sig1after, sig2after)));
%figure; plot(abs(xcorr(sig1after, sig2after)));
time(i);
end




figure;
title(["Взаимняа корреляция"]);
subplot(2,1,1);
plot(crosscoerr1);
subplot(2,1,2);
plot(crosscoerr2);

figure; plot(unwrap(angle(sig_est_interp1)))
figure; plot(unwrap(angle(sig_est_interp2)))

sw1=[];
for i=1:length(s1)
sw1(2*i-1)=real(s1(i));
sw1(2*i)=imag(s1(i));
end

sw2=[];
for i=1:length(s2)
sw2(2*i-1)=real(s2(i));
sw2(2*i)=imag(s2(i));
end

fid=fopen("/home/dmirty/MatlabProjects/GenSig/Equalized Signals/Eu21_eq.bin", "w");
fwrite(fid, sw1, "float");
fclose(fid);
fid=fopen("/home/dmirty/MatlabProjects/GenSig/Equalized Signals/Eu10_eq.bin", "w");
fwrite(fid, sw2, "float");
fclose(fid);

% fid= fopen("/home/dmirty/MatlabProjects/GenSig/Equalized Signals/Eu21_eq.bin", "r");
% Data21_1=fread(fid, "float");
% Data21_1=complex(Data21_1(1:2:end), Data21_1(2:2:end));
% fclose(fid);


figure; plot(unwrap(angle(signal1(1:L).*conj(s1(1:L)))));
% figure; plot(unwrap(angle(Data21_1(1:L).*conj(s1(1:L)))));