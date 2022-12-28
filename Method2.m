clc;
close all;
clear all;

Fs = 1.2e6;
df01=-3096.7; % оценка частотной отстройки для каждого из сигналов, данные от СТЦ
df02=-176037;

fileID = fopen('psp.bin', 'r'); 
psp = fread(fileID, 2*5*Fs, "short");
psp= complex(psp(1:2:end), psp(2:2:end));% чтение сигналов из файлов, приведение их в коплексный вид

fileID = fopen('Data21.bin', 'r'); 
signal1 = fread(fileID, 2*5*Fs, "short");
signal1= complex(signal1(1:2:end), signal1(2:2:end));
t= (1:length(signal1))'/Fs;
signal1 = signal1.*(exp(1i*2*pi*(df01)*t));

fileID = fopen('Data10.bin', 'r'); 
signal2 = fread(fileID, 2*5*Fs, "short");
signal2= complex(signal2(1:2:end), signal2(2:2:end));
t= (1:length(signal2))'/Fs;
signal2 = signal2.*(exp(1i*2*pi*(df02)*t));

% figure;
% title(["Корреляция сигнала sat21"]);
[~, S]= max(abs(xcorr(psp, signal1)));
%  subplot(2,1,1);
%  plot(abs(xcorr(psp, signal1)));
%  subplot(2,1,2);
%  plot(abs(xcorr(psp, psp)));
[~, S0]= max(abs(xcorr(psp, psp))); % вычисление рассогласование спутниковых сигналов и псп
lag1=S0-S;

% figure;
% title(["Корреляция сигнала sat10"]);
[~, S]= max(abs(xcorr(psp, signal2)));
%  subplot(2,1,1);
%  plot(abs(xcorr(psp, signal2)));
%  subplot(2,1,2);
%  plot(abs(xcorr(psp, psp)));
[~, S0]= max(abs(xcorr(psp, psp)));
lag2=S0-S;


signal2=signal2(1+lag2:end); %компенсация рассогласований
signal1=signal1(1+lag1:end); 


%%
tic
N=100; % кол-во кусочков

d=length(signal1)/N;% длина кусочков
d=floor(d);
A1=[];
pos1=[];
phases1=[];
for i=1:N
 corr=xcorr(psp(d*(i-1)+1:d*i),signal1);% корреляция кусочка ПСП и сигнала 
[A1(i), pos1(i)]=max(abs(corr)); % отыскание корреляционного пика 
phases1(i)=angle(corr(pos1(i)));% вычисление фазы пика
end

%phases1=rot90(rot90(phases1));% разворот массива фаз; пики идут в обратном порядке, см массив pos, взаимная корреляция без разворота хуже


phases1=repmat(phases1, d, 1); %интерполяция сигнала корреккции на ЧД сигналов 
phases1=phases1(:);

L=min(length(signal1), length(phases1));
phases1=phases1(1:L);
signal1=signal1(1:L);
phases1=phases1-phases1(1);

sig1=signal1.*exp(-j*phases1); % коррекция

figure; 
subplot(2,1,1);
plot(abs(xcorr(signal1, psp)));
title("Корреляция сигнала sat21 до коррекции");
subplot(2,1,2)
plot(abs(xcorr(sig1, psp)));
title("Корреляция сигнала sat21 после коррекции");

d=length(signal2)/N;% длина кусочков
d=floor(d);
A2=[];
pos2=[];
phases2=[];
for i=1:N
 corr=xcorr(psp(d*(i-1)+1:d*i),signal2);
[A2(i), pos2(i)]=max(abs(corr));
phases2(i)=angle(corr(pos2(i)));
end


%phases2=rot90(rot90(phases2));

phases2=repmat(phases2,d, 1);
phases2=phases2(:);

L=min(length(signal2), length(phases2));
phases2=phases2(1:L);
signal2=signal2(1:L);
phases2=phases2-phases2(2);

sig2=signal2.*exp(-j*phases2);

figure;
subplot(2,1,1)
plot(abs(xcorr(signal1, psp)));
title("Корреляция сигнала sat10 до коррекции");
subplot(2,1,2);
plot(abs(xcorr(sig1, psp)));
title("Корреляция сигнала sat10 после коррекции");

figure;
subplot(2,1,1)
crosscorr1=xcorr(signal1, signal2);
plot(abs(crosscorr1));
title("Взаимная корреляция до коррекции");
subplot(2,1,2);
crosscorr2=xcorr(sig1, sig2);
plot(abs(crosscorr2));
title("Взаимная корреляция после коррекции");

P1=max(abs(crosscorr1))/mean(abs(crosscorr1));
P2=max(abs(crosscorr2))/mean(abs(crosscorr2));
N=P2/P1;
%%

res1=[];
res2=[];
for i=1:9
c1=abs(xcorr(signal1(1:Fs*0.5*i), signal2(1:Fs*0.5*i)));
c2=abs(xcorr(sig1(1:Fs*0.5*i), sig2(1:Fs*0.5*i)));
res1(i)=max(c1)/mean(c1);
res2(i)=max(c2)/mean(c2);
end
figure;
subplot(2,1,1)
plot(res1);
title("Корреляция от длительности записи до коррекции");
subplot(2,1,2)
plot(res2);
title("Корреляция от длительности записи после коррекции");

figure; plot(abs(xcorr(psp(1:d), signal2)))
figure; plot(unwrap(-1*phases1))
figure; plot(unwrap(-1*phases2))
toc
