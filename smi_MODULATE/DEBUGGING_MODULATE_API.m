%% 全局变量
clc;
clear all;
close all;
fs = 100000;  % 采样率
N = 8000;  
fv = 30;  % 震动频

alpha = 5;
dir = 1; % 方向
C = [0.2]; 
h = 300; % 调制深度
fm = 10000;  % 调制频率
gamma = 0;  % 调制初相位
beta = 1;
n = 0 ;  % 翻倍次数，
arr = [14,1244,1252,1260,3745,3754,3759,6245,6254,6261];
arr1 = [1252,3754,6254];  % 方向

%% 产生自混合信号
figure(1);
subplot(7,1,1);
[t, lambda, L0, Lt, phi0, phiF, p, c] = MOVE_API_HARMONIC(fs, N, fv, C, alpha);
% [t, lambda, L0, Lt, phi0, p, c, phiF] = MOVE_API_ALEATORY_LOAD(fs, N, C, alpha);
p_init = p;
[p,h] = SMI_API_MODULATE(beta,phiF,h,fm,gamma,t);  % 调制深度/调制频率/调制信号初始相位

p = awgn(p,10);  % 10db，加高斯噪声
% p = p .* (1+0.5*cos(2*pi*75*t));  % 给自混合信号加包络，加了一个幅值为0.2，频率为75的包络

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
fshift = (-N/2:N/2-1)*(fs/N);  % 平移后信号的频域范围
p_ = fftshift(p_);  % fftshift将零频分量移动到数组中心，重新排列
% amp1 = abs(p_) * 2 / N ;
amp1 = abs(p_);
plot(amp1);
title("平移后频域信号（未更改频域范围）");
subplot(7,1,3);
plot(fshift,amp1);
title("平移后频域信号（更改频域范围）");
%----------------------------
% f = fs / N * (0 : 1 : N-1);
% p_ = fft(p);
% amp1 = abs(p_);
% plot(f,amp1);
% title("平移后频域信号（未更改频域范围）");
%% 取出谐波成分，将平移后的频谱置零取出(丽萍学姐论文)
pp_ = takeHarmonicComponent2(p_,fm,fshift,1);  % 这个范围还是不能太小嗷，6500-7500误差就变大了
% pp_ = takeHarmonicComponent(p_,fs,39960,119960);
subplot(7,1,4);
amp2 = abs(pp_) * 2 / N ;
plot(fshift,amp2);
title("取出的一次谐波成分（更改了频域范围）");

subplot(7,1,5);
p1 = ifft(ifftshift(pp_));
p1 = imag(p1);  % 呃。。。文章里一次谐波取的是虚部，这里为啥一次谐波取实部啊，所以这是sin
plot(p1);
title("一次谐波时域信号")

pp_ = takeHarmonicComponent2(p_,fm,fshift,2);
% pp_ = takeHarmonicComponent(p_,fs,119960,199960);
subplot(7,1,6);
p2 = ifft(ifftshift(pp_));
p2 = real(p2);
plot(p2);
title("二次谐波时域信号")

%% 条纹翻倍
figure(2);
subplot(7,1,1);
if besselj(2,2*h)==0 && besselj(1,2*h)==0
    error('贝塞尔函数0值处，无意义')
end
% 一类besselj、二类bessely、三类besselh
p1 = p1 ./ -besselj(1,2*h);  % 这里要删除掉对应的贝塞尔函数值吗
p2 = p2 ./ besselj(2,2*h);
n = n;  % 翻倍次数
times = 2^n;
[phiF_wrapped,Pn,Pnn] = SMI_API_TANDOUBLE_specific(n,p1,p2);
phiF_wrapped = dir * phiF_wrapped;
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
% Lt_reconstruct = Lt_reconstruct - mean(Lt_reconstruct); % 简写振动和余弦调制振动，重构后加上幅值A
% Lt_reconstruct = Lt_reconstruct - 1.5 * lambda;  % 重构后的随机振动信号要加上幅值1.5的波长，这是为啥我页不知道
plot(Lt_reconstruct,'r')
title(['解包裹重构后的信号，C-reconstruct=', num2str(C)]);

%% 误差分析
subplot(7,1,4);
plot(Lt-Lt_reconstruct);
RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
title(['解包裹重构RMSE=', num2str(RMSE)])

%% 数条纹重构
% 找到所有的峰谷值（理论上包含跳变点）
subplot(7,1,5);
[top_p,loc_p] = findpeaks(phiF_wrapped_init);  % 找到所有大于0.1的峰值,这是因为驼峰区乱七八糟（当然这只针对该实验信号）
[top_v,loc_v] = findpeaks(-phiF_wrapped_init);  % 默认。 'minpeakheight',0.1, 'minpeakdistance',50
top_v = -top_v;
plot(phiF_wrapped_init);
hold on;
title("数条纹重构的phiF包裹相位")
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

plot(Lt_reconstruct_beforeinternp);
title("Lt-before-internp")
hold on;
Lt_reconstruct = interp1(trigger, Lt_reconstruct_beforeinternp, 1:N, "spline");  % 这个插值，在开始后结束位置很容易出现大误差
% Lt_reconstruct = Lt_reconstruct - mean(Lt_reconstruct);  % 不开这个误差差这么多！

subplot(7,1,7);
plot(Lt);
hold on;
plot(Lt_reconstruct);
RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
title(['数条纹重构RMSE=', num2str(RMSE)]);






























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

%% 取出谐波信号函数,使用fshift前的坐标（即仅fft后的坐标），所以需要用到对应关系进行转换
function pp_ = takeHarmonicComponent(p_,fs,X,Y)  % p_为输入频谱。其中XY分别为输入fft频谱中的起始点和终点。
     N = length(p_);
%      X = (X-1)*fs/N - 20*N;  % 坐标转换
%      Y = (Y-1)*fs/N - 20*N;  % 坐标转换
%      X = (X+20*N)*N/fs + 1;  % 坐标转换
%      Y = (Y+20*N)*N/fs + 1;  % 坐标转换    
     pp_ = p_;
     pp_(1:X-1) = 0;
     pp_(Y+1:end) = 0;
     pp_(N/2-(Y-X)/2:N/2+(Y-X)/2) = pp_(X:Y);  % 这里就是移动到pp_的中间位置 
     pp_(X:Y) = 0;
end

%% 取出谐波信号函数,使用fshift的坐标
function pp_ = takeHarmonicComponent2(p_,fm,fshift,n)  % p_为fftshift后的频谱。n为需要的谐波的阶数
     N = length(p_); 
     pp_ = p_;
     X = find(fshift==n*fm)-1;  % fshift中移动n个fs，对应p_的位置
     step = (find(fshift==fm)-1-N/2)/2;
     
     pp_(1:X-step-1) = 0;
     pp_(X+step+1:end) = 0;
     
     pp_(N/2-step:N/2+step) = pp_(X-step:X+step);   
     pp_(X-step:X+step) = 0;
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