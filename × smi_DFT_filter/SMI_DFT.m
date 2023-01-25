clc;
clear all;
close all;
%% DFT滤波，虽然有滤波的效果，但是效果是真差啊
%% 产生自混合信号
figure(1);
fs = 204800;  % 采样率，即fs(s)采一个点。
N = 4096;  % 采样点——采样点可能会报错，原因在于某些采样点让direction出错了，如4000。原因未知
fv = fs/256;  % 震动频率
C = 0.7;
alpha = 3;
[t, lambda, L0, Lt, phi0, p] = SMI_API(fs, N, fv, C, alpha);
% p = p + randn(size(p));  % 添加噪声
plot(p);
title("加噪声的自混合信号");


%% 得到重构所需的相关信息
% [top_p, loc_p, top_v, loc_v, top_r, loc_r, direction] = SMI_FRINGE(p,N);
% direction = -direction;  % 如果初始震动用的cos，那方向是反的，一定要注意!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%% 傅里叶频谱滤波
p_ = fft(p, N);
amp1 = abs(p_) * 2 / N;

f = fs / N * (0 : 1 : N-1);  % Fs/N就是这个频谱中的最小频率间隔！！！！！所以N越大，分辨率会越高
f_norm = f;  % 储存归一化后的f
f_max = max(f);
f_min = min(f);
for i =1:length(f)
    f_norm(i) =  0 + (f(i) - f_min)/(f_max - f_min) * (1 - 0);  % 把数据归一化至0~1之间！！
end

figure(2);
plot(f_norm(1:N/2), amp1(1:N/2));  % fft算法默认是双边谱,通常我们只取一半
title("自混合信号的频谱");

y = @(x) 2 * x * (f_max - f_min) + f_min;
pp_ = zeros(1,length(p_));  % 用来储存滤波后的信号

% for i = 0 : N-1
%     if ((fs / N * i) > y(0.2777)) && ((fs / N * i) < y(1 - 0.2777)) % || (fs / N * i) > (fs - 40) && (fs / N * i) < (fs - 20)
%         pp_(i+1) = 0;
%     else
%         pp_(i+1) = p_(i+1);
%     end
% end

for i = 0: length(f_norm)-1
    if f_norm(i+1) > 0.28 && f_norm(i+1) < 1 - 0.28 %|| f_norm(i+1) < 0.0014
        pp_(i+1) = 0;
    else 
        pp_(i+1) = p_(i+1);
    end
end

amp2 = abs(pp_) * 2 / N;
figure(3);
plot(f_norm(1:N/2), amp2(1:N/2));

ppp_ = real(ifft(pp_, N));
figure(4);
plot(ppp_);

RMSe = sqrt(sum((p - ppp_).^2)/N);



