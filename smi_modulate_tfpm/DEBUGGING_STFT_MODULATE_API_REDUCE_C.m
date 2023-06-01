%% 主要想利用该方法解决两个问题，一个是拓展调制的C的范围（力图使调制能工作在C比较大的时候）
%% 第二是解决了TFPM算法中无法判断方向的问题
%% 调制+TFPM的东西，应该是和参数设置有关系，下一步是实现自动化（画图时根据T,F来画）。
%% 第二步就是把降C的工作引入到其中。
%% 第三步就是去直流，现在好像就差 去直流了，引入就结束了！！！ 

%% 经过多次仿真，fs/fm需要为整数才能实现频谱的正常搬移，但是一定要注意搬移距离至少要大到有足够间隔
%   如果频谱搬运正常，但是逆变换后的图像很乱，可能是窗长设置问题，窗长设置为整数比较好
%   如果频谱搬运正常，但是逆变换后的图像是一条直线，可能是取实部或虚部的问题
%   如果出现一次谐波清楚，二次谐波淡，可能是调制深度的问题
%   如果出现时频谱没有按照调制深度进行搬移，可能是调制深度的问题
%   重构的误差（可能重构的Lt是波折的），和调制深度有很大关系
%% 全局变量
clc;
clear all;
close all;
fs = 100000;  % 采样率
N = 8000;  
fv = 30;  % 震动频
alpha = 5;
dir = 1; % 方向

h = 300; % 调制深度
fm = 10000;  % 调制频率
gamma = 0;  % 调制初相位

windowLength = 500; % 窗长
%% 产生自混合信号
figure(1);
subplot(7,1,1);
t = (0:N-1)/fs;  % 采样时间，设N=10, fs=200，即采样了0.05s，t为[0...0.045]
lambda = 650e-9;  % 波长
% A = 15*97.5*10^-9;  % 幅值（✔）默认幅值为0.3*lambda
A = 3*lambda;  % 这个振幅还不能设置那么小，否则时频谱看着有点怪
L0 = 20 * lambda;  % 外腔距离（✔） 
Lt = A.* sin(2*pi*fv*t);  % L0为标称位置，外腔长度
beta = 1;  % the amplitude of selfmixing signal
phi0 = 4*pi*(L0+Lt)/lambda;

p = zeros(1,N);
C = [2.2,2.2];
 % C的变化是一个正弦曲线，不能随机数！
C_lower = C(1);
C_upper = C(2);
% 这个乘和加保证了c的上下限，这里可以设置变换的周期！！但是这个变换周期需要长一点否则会报错！！
x = linspace(0, 3*pi, N);
c = (C_upper-C_lower)/2 * cos(x) + (C_upper - (C_upper-C_lower)/2);
% plot(x,c);

for i = 1:N 
    C = c(i);
    phiF(i) = solve_phiF(C, phi0(i), alpha);
    p(i) = beta * cos(solve_phiF(C, phi0(i), alpha));  % 遍历所有的phi0
end
p_init = p;
[p,h] = SMI_API_MODULATE(beta,phiF,h,fm,gamma,t);  % 调制深度/调制频率/调制信号初始相位
% [phi0,h] = SMI_API_MODULATE2(phi0,300,80000,0,t);
% for i = 1:N 
%     C = c(i);
%     p(i) = beta * cos(solve_phiF(C, phi0(i), alpha));  % 遍历所有的phi0
% end
% p = awgn(p,30);  % 10db，加高斯噪声
% p = p .* (1+0.2*cos(2*pi*75*t));  % 给自混合信号加包络，加了一个幅值为0.2，频率为75的包络
plot(p);
hold on;
title("自混合信号");

%% 傅里叶变换看频谱
figure(1);
subplot(7,1,2);
% w = hamming(N);
f = fs / N * (0 : 1 : N-1);  % Fs/N就是这个频谱中的最小频率间隔！！！！！所以N越大，分辨率会越高
%----------------------------
p_ = fft(p);
% 平移频域信号
% fshift = (-N/2:N/2-1)*(fs/N);  % 平移后信号的频域范围
% p_ = fftshift(p_);  % fftshift将零频分量移动到数组中心，重新排列
% amp1 = abs(p_) * 2 / N ;
amp1 = abs(p_);
plot(amp1);
title("平移后频域信号（未更改频域范围）");
% subplot(5,1,3);
% plot(fshift,amp1);
% title("平移后频域信号（更改频域范围）");
f2N = @(x) N/fs * x + 1;  % 映射了从频域到N的对应关系

p_([f2N(fm),f2N(2*fm),f2N(3*fm)])=0;
p_([f2N(fs-fm),f2N(fs-2*fm),f2N(fs-3*fm)])=0;
amp2 = abs(p_);
subplot(7,1,3);
plot(amp2)
p = ifft(p_);

%% 全局变量
overlapLength = floor(windowLength * 0.9);  % OverlapLength后为指定的重叠长度
window = hamming(windowLength, "periodic");  % 使用汉明窗作为滑动的窗口
fftLength = 5*windowLength;  % 每个时刻傅里叶变换的长度
V = 0.65; % 抑制因子2

%% padding,逆变换不可避免的会让信号变短，要padding
judge = true; 
padding = windowLength;  % 逆变换不可避免的会让信号变短，要padding
while judge
    if mod(((2 * padding + N - overlapLength) / (windowLength - overlapLength)), 1) == 0
        judge = false;
    else
        padding = padding - 1;
    end
end
p = wextend('1D','sym',p,padding);  % 数据延拓
% p = [zeros(1,padding) p zeros(1,padding)];

%% 先不进行抑制，看看时频谱
figure(2);
[TF,F,T] = stft(p, fs, 'Window', window, 'OverlapLength', overlapLength, 'FFTLength', fftLength);  % （TF每一列就是每一个时间维度的频率）
subplot(4,2,[1,3]);
% 这样画图是以(-fs/2,fs/2)为纵轴, 对应(0 ,length(TF))的初始纵轴
mesh(T,F,abs(TF)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' ); 
view(0,90); % 设置初始视角为俯视角
% mesh(angle(TF));
title('抑制前');

%% 一次谐波时频谱
component1 = [fm-fm/2,fm+fm/2];
component2 = [2*fm-fm/2,2*fm+fm/2];
component3 = [3*fm-fm/2,3*fm+fm/2];
% TF1(2*windowLength+1:3*windowLength,:) = TF(1*windowLength+1:2*windowLength,:);  % 取出一次谐波
TF1 = takeHarmonicComponent2(TF,fs,fm,component3(1),component3(2));
TF_curb =TF1;
% 进行时频抑制
TF_curb = TF_inhibit1(TF1,V);
subplot(4,2,[2,4]);
mesh(T,F,abs(TF_curb)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' );
view(0,90);
title('一次谐波时频谱');


figure(1);
subplot(7,1,4);
p = (istft(TF_curb, fs, 'Window', window, 'OverlapLength', overlapLength, 'FFTLength', 5*windowLength))';
p_init1 = p;  % 用p_init1存储未归一化的信号
p_init1 = p_init1(padding+1:end-padding);  % ✔去掉之前的padding
p1 = imag(p_init1);
% p1 = 2 * (p1 - -abs(hilbert(-p1)))./(abs(hilbert(p1)) - -abs(hilbert(-p1))) - 1;
plot(p1);
title('一次谐波时域信号');

%% 二次谐波时频谱
figure(2);
TF2 = takeHarmonicComponent2(TF,fs,fm,component2(1),component2(2));
TF_curb = TF2;
% 进行时频抑制
TF_curb = TF_inhibit1(TF2,V);

subplot(4,2,[6,8]);
mesh(abs(TF_curb)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' );
view(0,90);
title('二次谐波时频谱');

figure(1);
subplot(7,1,5);
p = (istft(TF_curb, fs, 'Window', window, 'OverlapLength', overlapLength, 'FFTLength', 5*windowLength))';
p_init2 = p;  % 用p_init2存储未归一化的信号
p_init2 = p_init2(padding+1:end-padding);  % ✔去掉之前的padding
p2 = real(p_init2);
% p2 = 2 * (p2 - -abs(hilbert(-p2)))./(abs(hilbert(p2)) - -abs(hilbert(-p2))) - 1;
plot(p2);
title('二次谐波时域信号');


%% （去直流）
subplot(7,1,6);
% [top_p, loc_p, top_v, loc_v, top_r, loc_r, direction] = SMI_API_FRINGE(p1,N);
% [p1] = SMI_API_ELIMINATE_DC1(p1,direction,"time");
% [p2] = SMI_API_ELIMINATE_DC2(p_init,direction,"time");
% [p11,p22] = SMI_API_evenlopEXTRACT_HT_PRO(p1,N);
% [p33,p44] = SMI_API_evenlopEXTRACT_HT_PRO(p2,N);
% p1 = p11;
% p2 = p33;
% % plot(p2);
% title("去直流后的自混合信号");

%% 时频抑制
% windowLength = 128; % 窗长
% V = 0.65; % 抑制因子
% [T,F,TF,TF_curb,p1] = SMI_API_TFPM(p1,N,fs,windowLength,V);
% [T,F,TF,TF_curb,p2] = SMI_API_TFPM(p2,N,fs,windowLength,V);
% plot(p1);
% subplot(7,1,7);
% plot(p2);

%% 重构
figure(3);
subplot(5,1,1);
p1 = p1 ./ -besselj(1,2*h); 
p2 = p2 ./ besselj(2,2*h);
phiF_wrapped = atan(p1./p2);
% 方向问题
phiF_wrapped = dir * phiF_wrapped;
plot(phiF_wrapped);

%% 解包裹重构
subplot(5,1,2);
for i = 2:N
    if (phiF_wrapped(i) - phiF_wrapped(i-1) > pi/2)
        phiF_wrapped(i:end) = phiF_wrapped(i:end) -  pi;
    elseif (phiF_wrapped(i) - phiF_wrapped(i-1) < -pi/2)  % 找峰值点，碰到峰值则加pi
        phiF_wrapped(i:end) = phiF_wrapped(i:end) +  pi;            
    end
end

phiF_reconstruct = phiF_wrapped/1;
subplot(5,1,2);
plot(phiF_reconstruct);
title("phiF-reconstruct")

subplot(5,1,3);
phi0_reconstrut = phiF_reconstruct + 0*sin(phiF_reconstruct+atan(alpha));
Lt_reconstruct = phi0_reconstrut * lambda / (4 * pi);
Lt_reconstruct = sgolayfilt(Lt_reconstruct,1,11);
Lt_reconstruct = sgolayfilt(Lt_reconstruct,2,21);
Lt_reconstruct = sgolayfilt(Lt_reconstruct,3,31);

plot(Lt,'k');
hold on;
% Lt_reconstruct = Lt_reconstruct - mean(Lt_reconstruct); % 简写振动和余弦调制振动，重构后加上幅值A
% Lt_reconstruct = Lt_reconstruct + 1.5 * lambda;  % 重构后的随机振动信号要加上幅值1.5的波长，这是为啥我页不知道
plot(Lt_reconstruct,'r')
title(['解包裹重构后的信号，C-reconstruct=', num2str(C)]);

%% 误差分析
subplot(5,1,4);
plot(Lt-Lt_reconstruct);
RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
title(['绝对误差，RMSE=', num2str(RMSE)])











































%% subfunction
% 时频抑制算法1，嘎嘎好用
function TF_curb = TF_inhibit1(TF,V)
    weight1 = abs(TF)./max(abs(TF));
    weight2 = weight1;
    weight2( find(weight2 < V) ) = 0;  % 算法1：好用
    weight2 = weight2 ./ max(weight2);  % 抑制后按列归一化，用作权值矩阵
    TF_curb = TF .* weight2;  % 获得抑制后的信号
end

%% 时域取出谐波信号函数
function pp_ = takeHarmonicComponent(p_,X,Y)  % p_为输入频谱。其中XY分别为输入频谱的起始点和终点。N为采样点数
     N = length(p_);
     pp_ = p_;
     pp_(1:X-1) = 0;
     pp_(Y+1:end) = 0;
     pp_(N/2-(Y-X)/2:N/2+(Y-X)/2) = pp_(X:Y);   
     pp_(X:Y) = 0;
end

%% 根据指定范围从时频域取出谐波信号函数，根据[-fs/2,fs/2]对应[0,F]
function TF1 = takeHarmonicComponent2(TF,fs,fm,x1,x2)  % TF为作STFT后的信号。F为作时频谱后的F，即length(TF)
     max = fs/2;
     min = -fs/2;
     F = length(TF);
     % 把[-fs/2,fs/2]的数据归一化到[0,F]之间
     y1 = 0 + (x1 - min)/(max-min)*(F-0);
     y2 = 0 + (x2 - min)/(max-min)*(F-0);
     
     z1 = 0 + (-fm/2 - min)/(max-min)*(F-0);
     z2 = 0 + (fm/2 - min)/(max-min)*(F-0); 
     
     TF1 = zeros(size(TF));
     TF1(z1:z2,:) = TF(y1:y2,:);
end

%% 偶次幂函数！！！ 写了一个小时我真是醉了
function p_x = evenPower(x,p) % x为i
    if(x==0)
        p_x = p.^4-p.^2;  % i0
    else
        coff = 2^(x+1) -2;
        p_x = evenPower(x-1,p).^2 + 1 /2^coff *  evenPower(x-1,p);
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