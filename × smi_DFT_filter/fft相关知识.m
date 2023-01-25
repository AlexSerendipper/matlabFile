%% 对该信号继续傅里叶频谱分析:  信号为x = 5 + 7*cos(2*pi*15*t - 30*pi/180) + 3*cos(2*pi*40*t - 90*pi/180)的简单分析
clc;
clear all;
close all;
fs = 128;  % 采样频率
L = 256;  % 信号长度
t = (0:L-1)/fs;  % 时间
x = 5 + 7*cos(2*pi*15*t - 30*pi/180) + 3*cos(2*pi*40*t - 90*pi/180);
y = x + randn(size(t));  % 为信号添加噪声
subplot(3, 1, 1);
plot(t, y)
title("加了噪声的信号")

N = 2^nextpow2(L);  % 确定采样点数
Y = fft(y) ;  % fft后的赋值被压缩，乘 2/N 后恢复真实幅值
amp = abs(Y) * 2 / N;
f = fs/N * (0:1:N-1);  % Fs/N就是这个频谱中的最小频率间隔！！！！！所以N越大，分辨率会越高
subplot(3, 1, 2);
plot(f(1:N/2), amp(1:N/2));  % fft算法默认是双边谱,通常我们只取一半
title("幅值谱");

pha = angle(Y);
subplot(3, 1, 3);
plot(f(1:N/2), pha(1:N/2)); 
title("相位谱")

%% x = 5*sin(2*pi*10*t)+cos(2*pi*30*t),将信号频谱中20-40hz的波滤除，采样间隔dt=0.01,给出滤波前后的振幅谱和时域信号
clc;
clear all;
close all;
N = 512;
dt = 0.01;  % 即采样频率为fs = 100;
fs = 100;
t = (0:N-1)*dt;
x = 5*sin(2*pi*10*t)+cos(2*pi*30*t);

subplot(4, 1, 1);
plot(t, x);
title("滤波前时域信号");

subplot(4, 1, 2);
Y = fft(x);
amp = abs(Y) * 2 / N;  % 好像要乘在这个位置
f = fs/N * (0:1:N-1);  % Fs/N就是这个频谱中的最小频率间隔！！！！！所以N越大，分辨率会越高
plot(f(1:N/2), amp(1:N/2));
title("滤波前幅值谱");

YY = zeros(1,length(Y));  % 用来储存滤波后的信号

for i = 0 : N-1
    if (fs / N * i) > 20 && (fs / N * i) < 40 % || (fs / N * i) > (fs - 40) && (fs / N * i) < (fs - 20)  % 后半部分是那奎斯特频率的东西
        YY(i+1) = 0;
    else
        YY(i+1) = Y(i+1);
    end
end
subplot(4, 1, 3);
amp = abs(YY) * 2 / N;  % 这里一定要×啊！！！！！！！！！！
plot(f(1:N/2), amp(1:N/2));
title("滤波后幅值谱");

xx = real(ifft(YY));
subplot(4, 1, 4);
plot(t, xx);
title("滤波后逆变换回时域信号")

