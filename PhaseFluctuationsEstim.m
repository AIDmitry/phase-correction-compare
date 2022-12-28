function [ErrorSig, ErrorSig1, df, Fs2] = PhaseFluctuationsEstim(s1,s2,Fs,debug)
% Оценка фазовой нестабильности
% s1 и s2 - скомпенсированные на сдвиг по частоте, пропущенные через ФНЧ и синхронные
% во врмени сигналы
% Fs - частота дискретизации, Гц
% ErrorSig - exp(j*phi), phi - оценненая нестабильность фазы
% ErrorSig1 - exp(j*phi1), phi1 - оценненая нестабильность фазы за вычетом
% линейного набега
% df - оцененный сдвиг по частоте

% константы
Fs_out = 100; % Частота оценки фазовой нестабильности, Гц

Ts = 1/Fs;

%% сигнал рассогласования
e = s1.*conj(s2);
te = (1:length(e))'/Fs;

%% Понижение частоты дискретизации и его фильтрация
Q = 100;
Fs2 = Fs/Q;
e2 = resample(e,1,Q);
te2 = (1:length(e2))'/Fs2;

Q = 60;
Fs2 = Fs2/Q;
e2 = resample(e2,1,Q);
te2 = (1:length(e2))'/Fs2;

% второй этап понижения ЧД
Q = round(Fs2/Fs_out);
Fs2 = Fs2/Q;
Ts2 = 1/Fs2;

e2 = resample(e2,1,Q);
te2 = (1:length(e2))'/Fs2;

% фильтрация сигнала рассогласования
Fpass = 10; % Passband Frequency
Fstop = 30; % Stopband Frequency
% Fpass = 1; % Passband Frequency
% Fstop = 2; % Stopband Frequency
Dpass = 0.0057501127785; % Passband Ripple
Dstop = 0.0001; % Stopband Attenuation
dens = 40; % Density Factor
% minimum-order lowpass filter
% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs2/2), [1 0], [Dpass, Dstop]); % Оценка порядка оптимального КИХ-фильтра Parks-McClellan

% Calculate the coefficients using the FIRPM function.
b = firpm(N, Fo, Ao, W, {dens}); % Оценка порядка оптимального КИХ-фильтра Parks-McClellan
Hd1 = dfilt.dffir(b);
% fvtool(Hd1);

e2 = filtfilt(Hd1.Numerator,1,e2);% фильтрация

%% извлечение ошибки фазы

% Фаза сигнала рассогласования
phi = unwrap(angle(e2))*180/pi;
phi = phi-phi(1);

% сглаживание фазы
Nf = round(0.1*Fs2);
phif = filtfilt(ones(Nf,1)/Nf,1,phi);


phif(1:Nf-1) = phif(Nf);
phif(end-Nf+1:end) = phif(end-Nf);
phif = phif - phif(1);
% сигнал коррекции
ErrorSig = exp(1i*phif*pi/180);

phi_id = unwrap(angle(ErrorSig))*180/pi;
phi_id = phi_id-phi_id(1);

k = (0:length(phi_id)-1)'; order = 1;
px1 = polyfit(k,phi_id,order);
est1 = polyval(px1,k);
phi_id1 = phi_id-est1;
phi_id1 = phi_id1-phi_id1(1);
% оценка сдвига частоты
df = (est1(end)-est1(1))*(pi/180)/( (length(est1)-1)*Ts2*2*pi);
% сигнал коррекции без линейного набега
ErrorSig1 = exp(1i*phi_id1*pi/180); %% сигнал коррекции