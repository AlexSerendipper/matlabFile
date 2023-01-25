%% 该程序用于调试，其实也是非常好用且重要的
%% 据说能从频域中直接判断出方向,但是apFFT并不好用，我要试试内插FFT，提高一下频谱精度试试

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

subplot(7,1,2);
plot(Lt);
title("Lt");

%% 傅里叶频谱分析（主频是振动频率,主谐波频率可以计算振幅A）
% Fringe = p(1:1+255);  % 取出低条纹，结果大于0
Fringe = p(1:1+255);      
N = length(Fringe);  % 用的是拓延前的N，记住
subplot(7,1,3);
plot(Fringe);

figure(2);
f = fs / N * (0 : 1 : N-1);
Fringe_ = fft(Fringe);
subplot(4,1,1);
amp2 = abs(Fringe_);
plot(f(1:N/2), amp2(1:N/2));
subplot(4,1,2);
pha2 = angle(Fringe_);
plot(f(1:N/2), pha2(1:N/2)); 

direction  = solve_pha(Fringe);

%% 判断方向函数我试试
function direction = solve_pha(fringe)
    N = length(fringe);
    % f = fs / N * (0 : 1 : N-1);
    fringe_ = fft(fringe);
    amp2 = abs(fringe_);
    pha2 = angle(fringe_);
    [top, ind] = findpeaks(amp2(1:N/2),'minpeakheight',1) ;
    if length(ind) >= 2
        if 2 * pha2(ind(1)) - pha2(ind(2))  > 0
            direction(1:N) = 1;
        else
            direction(1:N) = -1;
        end   
    else
        direction(1:N) = 0;
    end
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