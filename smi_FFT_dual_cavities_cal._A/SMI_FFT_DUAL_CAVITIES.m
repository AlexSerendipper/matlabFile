%% 频谱滤波
clc;
clear all;
close all;
%% 产生自混合信号
figure(1);
fs = 200000;  % 采样率，即fs(s)采一个。
N = 50000;  % 采样点——采样点可能会报错，原因在于某些采样点让direction出错了，如4000。原因未知
fv1 = 200;  % 震动频率
fv2 = 150;
C = 0.7;
alpha = 3;
[t, lambda, L01, L02, Lt1, Lt2, phi01, p1, p2] = SMI_API(fs, N, fv1, fv2, C, alpha);
% p = p + randn(size(p));  % 添加噪声
subplot(3,1,1);
p1 = p1 - mean(p1);  % 消除直流分量
plot(p1);
title(['自混合信号,fv1= ', num2str(fv1)]);
subplot(3,1,2);
p2 = p2 - mean(p2);  % 消除直流分量
plot(p2);
title(['自混合信号,fv2= ', num2str(fv2)]);

%% 双波长傅里叶频谱分析（主频是振动频率,主谐波频率可以计算振幅A）
w = hamming(N);
f = fs / N * (0 : 1 : N-1);  % Fs/N就是这个频谱中的最小频率间隔！！！！！所以N越大，分辨率会越高
figure(2);
p1_ = fft(w'.* p1, N) * 2;  % 加窗傅里叶变换，具体为何这样写 我也不懂
amp1 = abs(p1_) * 2 / N ;
plot(f(1:N/2), amp1(1:N/2));  % fft算法默认是双边谱,通常我们只取一半
hold on;

p2_ = fft(w'.* p2, N) * 2;
amp2 = abs(p2_) * 2 / N ;
plot(f(1:N/2), amp2(1:N/2));
title('自混合信号的频谱'); 


% 根据主谐波阶计算振幅, 若当前设置100为第一个，2000为第20个
A1 = 50 * lambda / (4 * pi);  % 幅值真实值
A2 = 40 * lambda / (4 * pi);
A0 = @ (nd) (1 / 0.95) * ((lambda/(4 * pi))) * (nd + 1.2);  % 预测值
A0_2 = @ (nd) (1 / 0.96) * ((lambda/(4 * pi))) * (nd + 1.25);  % 预测值




