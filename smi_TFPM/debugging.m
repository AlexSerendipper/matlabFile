%% 
% 理论上已经推证了。使用未归一化信号直接重构的效果一样（甚至更好）,所以代码实现用归一化信号进行重构，用归一化的信号来判断方向
% 对于V，C=1.7以下设置为0.9都合适。当2>C>1.7,则设置0.5合适。我们发现椭圆加杠表示上抬的条纹，越抬会被当作分量被抑制
%% 全局变量
clc;
clear all;
close all;

fs = 200000;  % 采样率，即fs(s)采一个点。
N = 4000;  
fv = 200;  % 震动频
C = [.5,.5];
alpha = 4.6;

%% 全局变量
windowLength = 128; % 窗长
overlapLength = floor(windowLength * 0.9);  % OverlapLength后为指定的重叠长度
% overlapLength = floor(50);  % OverlapLength后为指定的重叠长度
window = hamming(windowLength, "periodic");  % 使用汉明窗作为滑动的窗口
% window = kaiser(windowLength,5);
fftLength = 5*windowLength;
% fftLength = 512;
M = C(1)/10; % 抑制因子
V = 0.9; % 抑制因子2

%% 产生自混合信号
figure(1);
subplot(7,2,1);
t = (0:N-1)/fs;  % 采样时间，设N=10, fs=200，即采样了0.05s，t为[0...0.045]
lambda = 650e-9;  % 波长
A = 2 * lambda;  % 幅值（✔）
L0 = 20 * lambda;  % 外腔距离（✔） 
Lt = A.* sin(2*pi*fv*t);  % L0为标称位置，外腔长度
beta = 1;  % the amplitude of selfmixing signal
phi0 = 4*pi*(L0+Lt)/lambda;

p = zeros(1,N);


 % C的变化是一个正弦曲线，不能随机数！
C_upper = C(1);
C_top = C(2);
% 这个乘和加保证了c的上下限
x = linspace(0, 3*pi, N);
c = (C(2)-C(1))/2 * cos(x) + (C_top - (C(2)-C(1))/2);
plot(x,c);
for i = 1:N 
    CC = c(i);
    p(i) = beta * cos(solve_phiF(CC, phi0(i), alpha));  % 遍历所有的phi0
end

% p = awgn(p,10); % 10db，加噪声在抑制是可以的，问题是加了噪声没法从原信号判断方向了。
plot(p);
p_init = p;
title(['自混合信号,C=',num2str(C)]);


%% 打padding
judge = true; 
padding = windowLength;  % 逆变换不可避免的会让信号变短，要padding
while judge
    if mod(((2 * padding + N - overlapLength) / (windowLength - overlapLength)), 1) == 0
        judge = false;
    else
        padding = padding - 1;
    end
end
p = [zeros(1,padding) p zeros(1,padding)];

%% 时频分析
% stft(p, fs, 'Window', window, 'OverlapLength', overlapLength, 'FFTLength', 5*windowLength);  % 如果不返回值就是直接画图
TF = stft(p, fs, 'Window', window, 'OverlapLength', overlapLength, 'FFTLength', fftLength);  % （wf每一列就是每一个时间维度的频率）
weight1 = abs(TF)./max(abs(TF));  % 要把这个矩阵当作权值使用，所以按列先归一化

% weight1(find(weight1 == 1)) = weight1(find(weight1 == 1)) - 1e-11;  % 是找到weight1中等于1的位置（理论上每列都有），并都减去num极小值，因为这个后边要用做分母，所以减去一个无穷小量，不要让他产生无穷大
% weight2 = (weight1.^64)./(1 - weight1.^64); % 算法1：构建一个函数，让原来频率分量大的更大，小的更小，windowLenght此处为抑制因子
% weight2 = weight1 .^ 2;

weight2 = weight1;
% weight2( find(weight2 < (std(weight2)).^(M)) ) = 0;  % 删除weight2中每列小于其标准差的数，这步意义不大。✔此处抑制效果越小抑制效果越好，越大抑制效果越差！
weight2( find(weight2 < V) ) = 0; 
weight2 = weight2 ./ max(weight2);
TF_curb = TF .* weight2;  % 获得抑制后的信号

size = size(TF);
[x,y] = meshgrid(1:size(2),1:size(1));
subplot(7,2,[2,4,6]);
mesh(x,y,abs(TF));
% mesh(x,y,angle(TF));
view(0,90)
title('抑制前')

subplot(7,2,[10,12,14]);
mesh(x,y,abs(TF_curb));
% mesh(x,y,angle(TF_curb));
view(0,90);
title('抑制后');

%% 逆变化
p = (istft(TF_curb, fs, 'Window', window, 'OverlapLength', overlapLength, 'FFTLength', 5*windowLength))';
p_init2 = p;  % 用p_init2存储未归一化的信号
p_init2 = p_init2(padding+1:end-padding);  % ✔去掉之前的padding
p_init2 = sgolayfilt(p_init2,4,31);  % ✔ 超参数，平滑滤波器, 这里的参数也需要细调！！！c小的时候可以设置4，31   c大的时候可以设置，2，21 

figure(1);
subplot(7,2,3);
plot(p_init2);
title('未归一化');


% subplot(7,2,5);
% p = 2 * (p - min(p))./(max(p) - min(p)) - 1;
% plot(p);
% title('归一化1');

subplot(7,2,7);
p = 2 * (p - -abs(hilbert(-p)))./(abs(hilbert(p)) - -abs(hilbert(-p))) - 1;
% p = p./abs(p);  % ×这个有问题，我自己乱推的
p = p(padding+1:end-padding);  % ✔去掉之前的padding
p = sgolayfilt(p,4,31);  % ✔ 超参数，平滑滤波器, 这里的参数也需要细调！！！c小的时候可以设置4，31   c大的时候可以设置，2，21 
plot(p);
title('归一化2');

%% 得到重构所需的方向信息
% [top_p, loc_p, top_v, loc_v, top_r, loc_r, direction] = SMI_API_FRINGE(p_init,N);  % 法1求方向，求的是初始信号的方向，所以都能用
direction = SMI_TFPM_DIRECTION(p,N);  % 法2求方向，必须用归一化2的信号

% direction = -direction;  % 如果初始震动用的cos，那方向是反的，一定要注意
hold on;
% plot(direction);

%% hilbert transform
inverse_hb = (imag(hilbert(p_init2))) .* direction;  % 得到sin ，使用未归一化的信号，精度更高
% inverse_hb = (imag(hilbert(p))) .* direction;
% plot(inverse_hb,"b");
% hold on;
% plot(p,'r');
% title("blue,inverse-HT and orginal-HT","red,self-mixing signal");

subplot(7, 2, 9);
phiF_wrapped = atan(inverse_hb./p_init2);  %  使用未归一化的信号，精度更高
plot(phiF_wrapped);
hold on;
title("phiF-wrapped")

%% arctan相位解包裹
for i = 2:N
    if (phiF_wrapped(i) - phiF_wrapped(i-1) > pi/2)
        phiF_wrapped(i:end) = phiF_wrapped(i:end) -  pi;
    elseif (phiF_wrapped(i) - phiF_wrapped(i-1) < -pi/2)  % 找峰值点，加pi
        phiF_wrapped(i:end) = phiF_wrapped(i:end) +  pi;            
    end
end

phiF_reconstruct = phiF_wrapped;
% subplot(7, 1, 6);
% plot(phiF_reconstruct)
% title("重构后的phiF")

%% 计算一下抑制后，C的值
[C_reconstruct, alpha_reconstruct] = SMI_API_ESTIMATE_C(phiF_reconstruct); 

%% 重构
subplot(7, 2, 11);
phi0_reconstrut = phiF_reconstruct + C_reconstruct*sin(phiF_reconstruct+atan(alpha));
Lt_reconstruct = phi0_reconstrut * lambda / (4 * pi);
Lt_reconstruct = Lt_reconstruct - mean(Lt_reconstruct);
plot(Lt_reconstruct,'b');
hold on;
plot(Lt, 'r');
title(['重构后的Lt,C=',num2str(C_reconstruct)]);

subplot(7, 2, 13);
plot(Lt-Lt_reconstruct)
title("误差")























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