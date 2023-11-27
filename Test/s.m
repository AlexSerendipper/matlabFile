% addpath(genpath('..\smi_api'));  % 添加
% rmpath(genpath('..\smi_api'));

subplot(7,1,1);
fs = 200000; % 采样率
N = 4000;  % 采样点
fv = 100;  % 震动频率
C = [0.5]; 
c=C(1);
alpha = 4;
[t, lambda, L0, Lt, phi0, p] = MOVE_API_HARMONIC(fs, N, fv, C, alpha);  % 1 简谐振动的自混合信号
% [t, lambda, L0, Lt, phi0, p] = MOVE_API_COS(fs, N, C, alpha);  % 2 余弦调制信号的自混合信号
% cut = 200;  % cut降采样，输入一个能被N整除的数，将N分为N/cut段  % 3
% [t, lambda, L0, Lt, phi0, p] = MOVE_API_ALEATORY(fs, N, cut, C, alpha);  % 3 产生随机振动时，方向×负
% [t, lambda, L0, Lt, phi0, p] = MOVE_API_ALEATORY_LOAD(fs, N, C, alpha);  % 4 加载储存好的随机振动，方向×负
plot(p);
% p = awgn(p,40);  % 10db，加高斯白噪声
% p = p .* (1+0.2*cos(2*pi*75*t));  % 给自混合信号加包络，加了一个幅值为0.2，频率为75的包络
hold on;
title(['自混合信号，C=', num2str(c)]);