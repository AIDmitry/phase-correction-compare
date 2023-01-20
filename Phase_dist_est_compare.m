% Generate signal

clc;
clear all;
close all;

SymbolRate = 1e6; % символьная скорость исходного информационного сигнала
Duration = 5; % длительность сигнала
Beta = 0.25;
Modulation=4; % количество точек модуляционного созвездия 
SNR_db=-2; %% ОСШ информационного сигнала

Span =40; % длина ИХ RRC фильтра в символах
Sps = 12;% количество отсчетов в одном символе ИХ 
R=2/9; % кодовая скорость

SymNumbers = randi(Modulation, 1, Duration*SymbolRate); % формирование символов информационного сигнала
ModSymbols = qammod(SymNumbers-1, Modulation);% модуляция сформированных символов созвездием 
%scatterplot(Signal);

Pulse = rcosdesign(Beta, Span, Sps); % формирование ИХ RRC фильтра с дадаными выше коэффициентами

FiltredSig = upfirdn(ModSymbols, Pulse, Sps); % прорежение модуляционных символов с коэффициентом SPS и последующая фильтрация

FiltredSig=FiltredSig(Sps*(Span-1) + 1 : end-Sps*(Span-1)); % убираем лишние модуляционные символы

ResultSignal= resample(FiltredSig, 1, 6); % понижаем ЧД до 2е6 

L=length(ResultSignal);

ResultInfoSig=ResultSignal; % запомниим информационный сигнал для вычисления корреляции 

Es=sum(abs(ResultSignal).^2)/L; % вычисляем энергию сигнала
% SNR = Eb/N0 N0 - спектральная плотность мощности АБГШ Eb - энергия, приходящаяся
% на передачу одного полезного бита
% D=N0/2 * Fs - дисперсия одной квадратуры шума, 2D=N0*Fs -суммарняа дисперсия
% N0=2D/Fs
% Eb=E/R/Fs где E - средний квадрат точек модуляционного
% созвездия
% SNR = E/(R*2D)
% отсюда D=E/(2*R*10^(SNR_db/10)

E=mean(qammod([0:Modulation-1], 4).*conj(qammod([0:Modulation-1], 4))); % средний квадрат точек модуляционного
% созвездия
D=E/(2*R*(10^((SNR_db)/10))); % общая дисперсия
sigma=sqrt(D/2); % СКО для одной квадратуры
noise = (randn(1, L) + j*randn(1, L))*sigma; 
ResultSignal=ResultSignal+noise; % добавление шума в сигнал


%% generate psp
T=0.01;
psp=main_psp_gen(T, SymbolRate); % генерация PSP длинной T с символьной скоростью, аналогичной информационному сигналу
SNR_PSP_db=-10; % ОСШ для PSP 

N=length(ModSymbols)/length(psp);%  повторение PSP до длительности информационного сигнала
N=round(N);                      %
fullPSP=repmat(psp, 1, N);       %

FiltredPSP = upfirdn(fullPSP, Pulse, Sps); % фильтрация с передискеитезицией на 2е6
FiltredPSP=FiltredPSP(Sps*(Span-1) + 1 : end-Sps*(Span-1));
FiltredPSP= resample(FiltredPSP, 1, 6);
p=1;
phaseShifts=cumsum((2*pi*rand(1, N/p) - pi)/5);% генерация скачков фазы

% phasenoise=phaseShifts;
% phasenoise=repmat(phasenoise, 2*length(psp), 1);
% phasenoise=phasenoise(:);

phasenoise=resample(phaseShifts, 2*length(psp)*p, 1);

L=length(FiltredPSP);

Espsp=sum(abs(FiltredPSP).^2)/L; % вычисление энергии сигнала

ResultPSP=FiltredPSP*sqrt(D/Espsp*10^(SNR_PSP_db/10)); % изменение энергии psp под требуемое осш

% Espsp=sum(abs(ResultPSP).^2)/L;
% snr=Espsp/D;

L=min(length(ResultPSP), length(ResultSignal));

ResultSignal_clear=ResultSignal+ResultPSP(1:L); % сохранение копии сигнала без фазовых искажений

ResultSignal=ResultSignal_clear.*exp(j*phasenoise(1:L)); % добавение фазовых искажений




%%

ResultSignal_clear = upfirdn(ResultSignal_clear, Pulse, Sps);
 
ResultSignal_clear=ResultSignal_clear(Sps*(Span-1) + 1 : end-Sps*(Span-1));% фильтрация сигналов чистого и с искажениями фазы 

ResultSignal_clear= resample(ResultSignal_clear, 1, 12);

ResultSignal = upfirdn(ResultSignal, Pulse, Sps);
 
ResultSignal=ResultSignal(Sps*(Span-1) + 1 : end-Sps*(Span-1));

ResultSignal= resample(ResultSignal, 1, 12);

 for i=1:length(psp) % прорежение psp для повышения ЧД до 2е6
     PSP0(2*i)=0;
     PSP0(2*i-1)=psp(i);
 end

%PSP0=psp;

corr=xcorr(PSP0, ResultSignal); % вычисление корреляции на psp 


 fullPSP0=repmat(PSP0, 1, N);
 fullPSP0=fullPSP0(1:length(ResultSignal));

  plot(abs(xcorr(ResultInfoSig(1:length(ResultSignal_clear)), ResultSignal))) % корреляция с готовым сигналом
  title("Корреляция с искаженным сигналом")
 figure;
 plot(abs(xcorr(ResultInfoSig(1:length(ResultSignal_clear)), ResultSignal_clear))) % корреляция с сигналом без фазовых искажений
 title("Корреляция с сигналом без искажений");
%%

Fs=2e6;

[ErrorSig2, ErrorSig2_, df, Fs2] = PhaseFluctuationsEstim(ResultSignal,fullPSP0,Fs,1); % оценка искажений между сигналом с искажениями и psp

sig_est2=ErrorSig2;

L = Fs/Fs2;
sig_est_interp2 = repmat(sig_est2,L,1);
sig_est_interp2 = sig_est_interp2(:);  %интерполяция сигнала на 2e6


ResultSignal_corrected1=ResultSignal;
 
L = min(length(sig_est_interp2), length(ResultSignal_corrected1));

sig_est_interp2 = sig_est_interp2(1:L);
ResultSignal_corrected1=ResultSignal_corrected1(1:L);

ResultSignal_corrected1=ResultSignal_corrected1.*conj(sig_est_interp2'); % коррекция 

corr_corrected1=xcorr(ResultInfoSig, ResultSignal_corrected1);

figure;
plot(abs(corr_corrected1))
title("После коррекции 1-ым способом");



%%

x=zeros(N, 2);


for i=1:N % нахождение корреляционных пиков как максимум корреляции на участке, где происходти наложение psp на 
[x(i, 1), ~]=max(abs(corr(length(psp)+2*length(psp)*(i-1):length(psp)+2*length(psp)*i))); % очередную свою копию суммарном сигнале
[~, x(i, 2)]=max(abs(corr(length(psp)+2*length(psp)*(i-1):length(psp)+2*length(psp)*i))); % получение позиции пика на очередном участке
x(i,2)=x(i,2)+length(psp)+2*length(psp)*(i-1)-1;%  получение позиции пика на длинне всей корреляции 
end

% phaseShifts_est=[];
% for i=1:N
% phaseShifts_est=[phaseShifts_est angle(corr(x(i, 2)))]; % создание массива фаз пиков
% end
phaseShifts_est=angle(corr(x(:, 2)));
phaseShifts_est=-1*rot90(rot90(phaseShifts_est));

phase_dist_est=[];
for i=1:N
   phase_dist_est=[phase_dist_est repmat(phaseShifts_est(i), 1, 2*length(psp))]; % создания сигнала оценки фазовых искажений
end

correction_signal=exp(-j*phase_dist_est(1:length(ResultSignal))); % создание сигнала корреккции

ResultSignal_corrected=ResultSignal.*correction_signal;% коррекция

corr_corrected=xcorr(ResultInfoSig, ResultSignal_corrected);

figure;
plot(abs(corr_corrected));

title("После коррекции 2-ым способом");

figure; plot((unwrap(phase_dist_est)));
title("оценка фазы методом 2");
figure;  plot((unwrap(angle(sig_est_interp2))))
title("оценка фазы методом 1");
figure;  plot((phasenoise))
title("сгенерированная фаза");