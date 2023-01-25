%% 使用加窗傅里叶变换以及apFFT对锯齿波幅度相位谱进行分析
close all;
clc;
clear all;
%% 加窗傅里叶变换和全相位傅里叶变换分析
figure(1);
fv = 20;  % 锯齿波频率
N = 100;
fs = 200;  % 采样频率
% t = (0:N-1)/fs;
% t = -N + 1: 2 * N -1;
t = (1:2*N-1)/fs; % 传入两倍的信号-1长度，就是实际产生的是两倍信号减一的长度，N还是N
y = sawtooth(2*pi*fv*t);
y1=y(N:2*N-1);  % N长的数据送入FFT处理（这就是两倍长度信号的后一段，）
y2=y(1:2*N-1);  % 2N-1长的数据送入APFFT处理（这就是两倍长度-1信号的整段）
plot(y);  

%% 普通加窗傅里叶变换
figure(2);
subplot(4,2,[1,3]);
w = hamming(N);
f = fs / N * (0 : 1 : N-1);  % Fs/N就是这个频谱中的最小频率间隔！！！！！所以N越大，分辨率会越高
y1_ = fft(w'.* y1, N);  % 加窗傅里叶变换，具体为何这样写 我也不懂
amp1 = abs(y1_);
plot(f(1:N/2), amp1(1:N/2));  % fft算法默认是双边谱,通常我们只取一半
title('windowFFT-amp');

subplot(4,2,[5,7]);
% pha1 = angle(y1_) * 180 / pi;
pha1 = angle(y1_);
plot(f(1:N/2), pha1(1:N/2));  % fft算法默认是双边谱,通常我们只取一半
title('windowFFT-pha');


%% apFFT预处理
w1 = hanning(N)';  % 1.构建一个N点的汉宁窗
w2 = conv(w1,w1);  % 2.汉宁窗对自己求卷积，得到2N-1点的卷积窗
w2 = w2 / sum(w2);  % 3.对卷积窗进行归一化
% y2 = w2 .* y(1:2*N-1);  % 4.将数据的2N-1项与和归一化卷积窗相乘
y2 = y2 .* w2;
y2 = y2(N:end) + [0, y2(1:N-1)];  % 5.固定步骤，得到经过全相预处理的N点序列

%% apFFT预处理结束
subplot(4,2,[2,4]);
f = fs / N * (0 : 1 : N-1);  % Fs/N就是这个频谱中的最小频率间隔！！！！！所以N越大，分辨率会越高
y2_ = fft(y2,N);
amp2 = abs(y2_);
plot(f(1:N/2), amp2(1:N/2));  % fft算法默认是双边谱,通常我们只取一半
title('apFFT-amp');

subplot(4,2,[6,8]);
% pha2 =mod(angle(y2_)*180/pi,360);
pha2 = angle(y2_);
plot(f(1:N/2), pha2(1:N/2));  % fft算法默认是双边谱,通常我们只取一半
title('apFFT-pha');