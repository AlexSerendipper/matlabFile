clc;
clear all;
close all;
fs = 10000; 
N = 10000;
fv = 10; 
alpha = 5;
n = 3 ;  % 翻倍次数，艹啊。2和3一样
% arr = [834,2501,4168,5834,7501,9168];  % 去点
arr = [1251,3751,6251,8751];
arr1 = [1251,3751,6251,8751];  % 方向
arr1 = arr;
%% 产生自混合信号
figure(1);
subplot(9,1,1);
t = (0:N-1)/fs;  % 采样时间，设N=10, fs=200，即采样了0.05s，t为[0...0.045]
lambda = 650e-9;  % 波长
% A = 15*97.5*10^-9;  % 幅值（✔）默认幅值为0.3*lambda
A = 0.3*lambda;
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
    phiF(i) = solve_phiF(C, phi0(i), alpha);
    p(i) = beta * cos(solve_phiF(C, phi0(i), alpha));  % 遍历所有的phi0
end

p_init = p;

[p,h] = SMI_API_MODULATE(beta,phiF,300,2000,0,t);  % 调制深度/调制频率/调制信号初始相位
% [phi0,h] = SMI_API_MODULATE2(phi0,300,80000,0,t);
% for i = 1:N 
%     C = c(i);
%     p(i) = beta * cos(solve_phiF(C, phi0(i), alpha));  % 遍历所有的phi0
% end
% p = awgn(p,40);  % 10db，加高斯噪声
% p = p .* (1+0.2*cos(2*pi*75*t));  % 给自混合信号加包络，加了一个幅值为0.2，频率为75的包络
plot(p);
hold on;
title("自混合信号");

%% 傅里叶变换看频谱
figure(1);
subplot(9,1,2);
% w = hamming(N);
f = fs / N * (0 : 1 : N-1);  % Fs/N就是这个频谱中的最小频率间隔！！！！！所以N越大，分辨率会越高
%----------------------------
p_ = fft(p);
% 平移频域信号
fshift = (-N/2:N/2-1)*(fs/N);  % 平移后信号的频域范围
p_ = fftshift(p_);  % fftshift将零频分量移动到数组中心，重新排列
% amp1 = abs(p_) * 2 / N ;
amp1 = abs(p_);
plot(amp1);
title("平移后频域信号（未更改频域范围）");
subplot(9,1,3);
plot(fshift,amp1);
title("平移后频域信号（更改频域范围）");
%----------------------------
% f = fs / N * (0 : 1 : N-1);
% p_ = fft(p);
% amp1 = abs(p_);
% plot(f,amp1);
% title("平移后频域信号（未更改频域范围）");
%% 取出谐波成分，将平移后的频谱置零取出(丽萍学姐论文)
pp_ = takeHarmonicComponent(p_,fs,6000,8000);  % 这个范围还是不能太小嗷，6500-7500误差就变大了
% pp_ = takeHarmonicComponent(p_,fs,39960,119960);
subplot(9,1,4);
amp2 = abs(pp_) * 2 / N ;
plot(fshift,amp2);
title("取出的一次谐波成分（更改了频域范围）");

subplot(9,1,5);
p1 = ifft(ifftshift(pp_));
p1 = imag(p1);  % 呃。。。文章里一次谐波取的是虚部，这里为啥一次谐波取实部啊，所以这是sin
plot(p1);
title("一次谐波时域信号")

pp_ = takeHarmonicComponent(p_,fs,8000,10000);
% pp_ = takeHarmonicComponent(p_,fs,119960,199960);
subplot(9,1,6);
p2 = ifft(ifftshift(pp_));
p2 = real(p2);
plot(p2);
title("二次谐波时域信号")

%% 条纹翻倍
figure();
if besselj(2,2*h)==0 && besselj(1,2*h)==0
    error('贝塞尔函数0值处，无意义')
end
% 一类besselj、二类bessely、三类besselh
p1 = p1 ./ -besselj(1,2*h);  % 这里要删除掉对应的贝塞尔函数值吗
p2 = p2 ./ besselj(2,2*h);
subplot(7,1,1);
n = n;  % 翻倍次数
times = 2^n;
[phiF_wrapped,Pn,Pnn] = SMI_API_TANDOUBLE_specific(n,p1,p2);
phiF_wrapped_init = phiF_wrapped; 
plot(phiF_wrapped);
hold on;
plot(p_init);

% figure(3);
% subplot(2,1,1);
% plot(Pn);
% subplot(2,1,2);
% plot(Pnn);

%% 解包裹重构
for i = 2:N
    if (phiF_wrapped(i) - phiF_wrapped(i-1) > pi/2)
        phiF_wrapped(i:end) = phiF_wrapped(i:end) -  pi;
    elseif (phiF_wrapped(i) - phiF_wrapped(i-1) < -pi/2)  % 找峰值点，碰到峰值则加pi
        phiF_wrapped(i:end) = phiF_wrapped(i:end) +  pi;            
    end
end
phiF_reconstruct = phiF_wrapped/times;
subplot(7,1,2);
plot(phiF_reconstruct);
title("phiF-reconstruct")

subplot(7,1,3);
phi0_reconstrut = phiF_reconstruct + C*sin(phiF_reconstruct+atan(alpha));
Lt_reconstruct = phi0_reconstrut * lambda / (4 * pi);
% Lt_reconstruct(1:padding)=[];
% Lt_reconstruct(end-padding+1:end)=[];
plot(Lt,'k');
hold on;
Lt_reconstruct = Lt_reconstruct - mean(Lt_reconstruct); % 简写振动和余弦调制振动，重构后加上幅值A
% Lt_reconstruct = Lt_reconstruct + 1.5 * lambda;  % 重构后的随机振动信号要加上幅值1.5的波长，这是为啥我页不知道
plot(Lt_reconstruct,'r')
title(['解包裹重构后的信号，C-reconstruct=', num2str(C)]);

%% 误差分析
subplot(7,1,4);
plot(Lt-Lt_reconstruct);
RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
title(['绝对误差，RMSE=', num2str(RMSE)])

%% 数条纹重构
% 找到所有的峰谷值（理论上包含跳变点）
subplot(7,1,5);
[top_p,loc_p] = findpeaks(phiF_wrapped_init);  % 找到所有大于0.1的峰值,这是因为驼峰区乱七八糟（当然这只针对该实验信号）
[top_v,loc_v] = findpeaks(-phiF_wrapped_init);  % 默认。 'minpeakheight',0.1, 'minpeakdistance',50
top_v = -top_v;
plot(phiF_wrapped_init);
hold on;
% scatter(loc_p,top_p);
% scatter(loc_v,top_v);

%% 当findpeaks无法精确定位实验信号的所有峰谷值时（峰谷值中包含翻转点，亦或是有一些异常点）,针对这些异常点的特殊处理
[loc_p,loc_v,top_p,top_v] = ditchReversePoing(arr,loc_p,loc_v,top_p,top_v);
scatter(loc_p,top_p);
scatter(loc_v,top_v);
hold on;

%% 确定方向3：手动确定方向（最方便）
dir = solve_dir(arr1,N);  % 输入翻转点arr数组(顺序输入) （方向有误乘-1即可）
direction = dir;
plot(dir)
%% 数条纹重构 
subplot(7,1,6);
step = lambda/4/times;  % 单根条纹对应的波长
% 用谷值点，我感觉适用于cos振动
% 用峰值点，我感觉适用于sin振动
trigger = loc_p;
Lt_reconstruct_beforeinternp = zeros(1,length(trigger));

for i = 2:length(Lt_reconstruct_beforeinternp)
    if direction(trigger(i))==direction(trigger(i-1))
        Lt_reconstruct_beforeinternp(i) = Lt_reconstruct_beforeinternp(i-1) + step*direction(loc_v(i));
    elseif direction(trigger(i))~=direction(trigger(i-1))
        Lt_reconstruct_beforeinternp(i) = Lt_reconstruct_beforeinternp(i-1);
    end
end
figure(1);
subplot(7,1,6);
plot(Lt_reconstruct_beforeinternp);
title("Lt-before-internp")

figure(2);
Lt_reconstruct = interp1(trigger, Lt_reconstruct_beforeinternp, 1:N, "spline");  % 这个插值，在开始后结束位置很容易出现大误差
Lt_reconstruct = Lt_reconstruct - mean(Lt_reconstruct);  % 不开这个误差差这么多！
plot(Lt);
hold on;
plot(Lt_reconstruct);
RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
title(['Lt(count fringe)&均方误差:RMSE=', num2str(RMSE)]);

%% 

subplot(7,1,7);
padding = N/100;
% plot(Lt(padding:end-padding+1)-Lt_reconstruct(padding:end-padding+1));
plot(Lt-Lt_reconstruct);
title(['均方误差:RMSE=', num2str(RMSE)]);



























%% subfunction ditch reversepoint，输入不想要的点（翻转点）的数组，在峰谷值中去掉翻转点
function [loc_p,loc_v,top_p,top_v] = ditchReversePoing(arr,loc_p,loc_v,top_p,top_v)
    for i=1:length(arr)
        aa = find(loc_p==arr(i));
        bb = find(loc_v==arr(i));
        flag1 = isempty(aa);
        flag2 = isempty(bb);
        if flag1==0
            loc_p(aa)=[];
            top_p(aa)=[];
        end
        if flag2==0
            loc_v(bb)=[];
            top_v(bb)=[];
        end
    end
end





%% subfunction direction，输入翻转点arr数组(顺序输入) （方向有误乘-1即可）
function direction = solve_dir(arr,N)
    arr1 = [1, sort([arr, arr+1]), N];
    temp = zeros(1,N);
    flag = 1;
    
    for i = 1:2:length(arr1)-1      
        temp(arr1(i):arr1(i+1)) = flag;
        flag = flag * -1;
    end
    direction = sign(temp);
end
%% 取出谐波信号函数,使用fshift坐标，所以用对应关系进行转换
function pp_ = takeHarmonicComponent(p_,fs,X,Y)  % p_为输入频谱。其中XY分别为输入fshift频谱中的起始点和终点。N为采样点数
     N = length(p_);
%      X = (X-1)*fs/N - 20*N;  % 坐标转换
%      Y = (Y-1)*fs/N - 20*N;  % 坐标转换
%      X = (X+20*N)*N/fs + 1;  % 坐标转换
%      Y = (Y+20*N)*N/fs + 1;  % 坐标转换    
     pp_ = p_;
     pp_(1:X-1) = 0;
     pp_(Y+1:end) = 0;
     pp_(N/2-(Y-X)/2:N/2+(Y-X)/2) = pp_(X:Y);   
     pp_(X:Y) = 0;
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