%% 该程序用于调试，其实也是非常好用且重要的
%% 据说能从频域中直接判断出方向

%% 全局变量
clc;
clear all;
close all;
fs = 20000;  % 采样率，即fs(s)采一个点。
N = 100000;  
fv = 0.4;  % 震动频
alpha = 5;

%% 产生自混合信号
figure(1);
subplot(7,1,1);
t = (0:N-1)/fs;  % 采样时间，设N=10, fs=200，即采样了0.05s，t为[0...0.045]
lambda = 650e-9;  % 波长
A = 3*10^-4;  % 幅值（✔）  % 初始值为2 * lambda
L0 = 20 * lambda;  % 外腔距离（✔） 
Lt = A.* sin(2*pi*fv*t);  % L0为标称位置，外腔长度
beta = 1;  % the amplitude of selfmixing signal
phi0 = 4*pi*(L0+Lt)/lambda;

p = zeros(1,N);
C = [0.1, 0.1];
 % C的变化是一个正弦曲线，不能随机数！
C_lower = C(1);
C_upper = C(2);
% 这个乘和加保证了c的上下限，这里可以设置变换的周期！！但是这个变换周期需要长一点否则会报错！！
x = linspace(0, 3*pi, N);
c = (C_upper-C_lower)/2 * cos(x) + (C_upper - (C_upper-C_lower)/2);
% plot(x,c);

for i = 1:N 
    C = c(i);
    p(i) = beta * cos(solve_phiF(C, phi0(i), alpha));  % 遍历所有的phi0
end
% p = awgn(p,10);  % 10db，加高斯噪声
% p = p .* (1+0.2*cos(2*pi*75*t));  % 给自混合信号加包络，加了一个幅值为0.2，频率为75的包络
plot(p);
hold on;
title("自混合信号")


%% 傅里叶频谱分析（主频是振动频率,主谐波频率可以计算振幅A）
Fringe = p(6000:6000+255);  % 取出低条纹，结果大于0
% Fringe = p(6000:1856+512);      
N = length(Fringe);  % 用的是拓延前的N，记住
subplot(7,1,2);
plot(Lt);
title("Lt");
subplot(7,1,3);
plot(Fringe);

%% 数据拓延至2N-1
subplot(7,1,4);
% padding = length(Fringe)/2;
% Fringe = wextend('1D','per',Fringe,padding);  % 数据延拓，一定要是周期拓延
% Fringe = wextend('1D','sym',Fringe,padding);  % 数据延拓，一定要是周期拓延
% Fringe(end) = [];  % 这个真不确定
% N = length(Fringe);


Fringe = per_wextend(Fringe);
Fringe = per_wextend(Fringe);
Fringe = Fringe(1:511);
plot(Fringe);
y2 = Fringe;

%% apFFT预处理
figure(2);
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














%% 周期拓延函数
function per = per_wextend(sig)  % 传入输入的信号，进行信号拓延
    sig_init = sig;
    mmin = min([sig_init(1),sig_init(2),sig_init(3)]);  % 三个点可能一样，所以多放几个点
    mmax = max([sig_init(1),sig_init(2),sig_init(3)]);
    ind = find(sig_init>mmin & sig_init<mmax);
    sig(ind(end-2):end) = [];
    
    per = [sig,sig_init];
end











%% subfunction （selfmixing-power）
function phiF = solve_phiF(C,phi0,alpha)  % 求解出每一个phi0对应的phiF
    if C<=1  % 每个phi0的解区间
        [phiF_min, phiF_max] = bounds1(C,phi0);
    else
        [phiF_min, phiF_max] = bounds2(C,phi0,alpha);
    end
    
    excessphaze_equation = @(phiF)phiF-phi0+C*sin(phiF+atan(alpha));
    
    if (excessphaze_equation (phiF_min)>0)  % 文章的解释为,phiF_min值可能刚好比零点大一点点点，这时候取phiF_min为近似零点
        phiF = phiF_min;
    elseif (excessphaze_equation (phiF_max)<0)
    	phiF = phiF_max;
    else
    	phiF = fzero(excessphaze_equation,[phiF_min,phiF_max]);
    end  
end

%---C < 1时的解区间函数------------------------------------------
function [phiF_min, phiF_max] = bounds1(C,phi0) 
    phiF_min = phi0-C;
    phiF_max = phi0+C;
end

%---C > 1时的解区间函数------------------------------------------
function [phiF_min, phiF_max] = bounds2(C,phi0,alpha)
persistent m; 
if isempty (m); m = 0; end
mlower = ceil ((phi0 + atan (alpha) + acos (1/C)- sqrt (C*C- 1))/(2*pi)- 1.5);
mupper = floor ((phi0 + atan (alpha)- acos (1/C) + sqrt (C*C- 1))/(2*pi)- 0.5);
if (m < mlower); m = mlower; end
if (m > mupper); m = mupper; end
phiF_min = (2*m+1)*pi + acos (1/C)- atan (alpha); 
phiF_max = (2*m+3)* pi- acos (1/C)- atan (alpha); 
end


%% subfunction2(reconstruction_T_N间转换)
% 该函数实现，将N转换为对应的t，也就是从T=1，N=10这些简单的时候推出来N(i)与t的关系
% 其中N为需要转换的点，N为采样点数，T为模拟时间
function  temp = N2T(temp1,N,T)
temp = zeros(1,length(temp1));    
    for i = 1:length(temp1) 
        temp(i) = T*(temp1(i)-1)/N;
    end
end

function temp = T2N(temp1,N,T)
temp = zeros(1,length(temp1));    
    for i = 1:length(temp1) 
        temp(i) = N*temp1(i)/T + 1;
    end
end