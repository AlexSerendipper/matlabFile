%% 频谱滤波
clc;
clear all;
close all;
%% 产生自混合信号
figure(1);
fs = 200000;  % 采样率，即fs(s)采一个点。
N = 50000;  % 采样点——采样点可能会报错，原因在于某些采样点让direction出错了，如4000。原因未知
fv = 200;  % 震动频率
C = 3;
alpha = 4.6;
[t, lambda, L0, Lt, phi0, p] = SMI_API(fs, N, fv, C, alpha);
% p = p + randn(size(p));  % 添加噪声
subplot(3,1,1);
p = p - mean(p);  % 消除直流分量
plot(p);
title(['加噪声的自混合信号,fv= ', num2str(fv)]);


%% 得到重构所需的相关信息
% [top_p, loc_p, top_v, loc_v, top_r, loc_r, direction] = SMI_FRINGE(p,N);
% direction = -direction;  % 如果初始震动用的cos，那方向是反的，一定要注意!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%% 傅里叶频谱分析（主频是振动频率,主谐波频率可以计算振幅A）
% figure(2);
% f = fs / N * (0 : 1 : N-1);  % Fs/N就是这个频谱中的最小频率间隔！！！！！所以N越大，分辨率会越高
% p_ = fft(p, N);
% amp1 = abs(p_) * 2 / N;
% plot(f(1:N/2), amp1(1:N/2));  % fft算法默认是双边谱,通常我们只取一半
% title('自混合信号的幅度谱');

figure(2);
w = hamming(N);
f = fs / N * (0 : 1 : N-1);  % Fs/N就是这个频谱中的最小频率间隔！！！！！所以N越大，分辨率会越高
p_ = fft(w'.* p, N) * 2;  % 加窗傅里叶变换，具体为何这样写 我也不懂
amp1 = abs(p_) * 2 / N ;
plot(f(1:N/2), amp1(1:N/2));  % fft算法默认是双边谱,通常我们只取一半
title('自混合信号的幅度谱');

% figure(2)
% p_ = fft(p);
% fshift = (-N/2:N/2-1)*(fs/N);  % 平移后信号的频域范围
% p_ = fftshift(p_);  % fftshift将零频分量移动到数组中心，重新排列
% amp1 = abs(p_) * 2 / N ;
% plot(fshift,amp1);
% title("取出的一次谐波成分（更改了频域范围）");

% 根据主谐波阶计算振幅, 若当前设置100为第一个，2000为第20个
A = 50 * lambda / (4 * pi);  % 真实值
A0 = @ (nd) (1 / 0.95) * ((lambda/(4 * pi))) * (nd + 1.2);  % 预测值
A0_2 = @ (nd) (1 / 0.96) * ((lambda/(4 * pi))) * (nd + 1.25);  % 预测值




