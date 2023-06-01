%% 该程序用于调试，其实也是非常好用且重要的
%% 啊艹，数条纹实验信号比仿真信号还好做，我日了你大爸的

%% 全局变量
clc;
clear all;
close all;
% fs = 100000;  % 采样率，即fs(s)采一个点。采样率越大，调制频率就要越大，否则会失真，默认采样率10000==>调制频率2000，二者成倍数关系
% N = 10000;  % 调制默认值为10000
fv = 40;  % 震动频，调制默认值为10
alpha = 5;
n = 3;
% arr = [1251,3751,6251,8751];  % 去点
% arr1 = [1251,3751,6251,8751];  % 方向

%% 实验信号
figure(1);
subplot(9,1,1);

% path =  'D:\matlab save\smi_实验信号\EOM\9v 10hz 2k 200k.csv'; 
% path =  'D:\matlab save\smi_实验信号\EOM\11v 10hz 2k 200k.csv';  % 150000-50000,比较垃圾啊，不太行这个信号 
% path =  'D:\matlab save\smi_实验信号\EOM\11v 40hz 4k 200k.csv';  % 312100-50000，这个还可以，看看有没有更好的
path =  'D:\matlab save\smi_实验信号\EOM\9v 40hz 4k 200k.csv';  % 390000-50000，这是最好的目前
M = 390000; N = 50000;  [t, p, fs] = MOVE_API_EXPERIMENT(M, N, path);  % 1 加载.csv文件，从M点处开始取N个点
% load('xxx.mat');  % 2. 加载.mat文件
lambda = 650e-9;  % 波长
% p = sgolayfilt(p,1,11);
% p = sgolayfilt(p,2,21);
% p = sgolayfilt(p,3,31);

% SMI_API_CORR_FILTER(p,smoothingfactor,threshold) 使用自相关能够有效去噪，但是在后续处理方向上反而不太好，可能是信号带散斑的原因
% p = -1 + (p - min(p))/(max(p) - min(p)) * 2;  % 还是需要归一化一下，否则无法求acosp
plot(p);
hold on;

%% 傅里叶变换看频谱
figure(1);
subplot(9,1,2);
% w = hamming(N);
f = fs / N * (0 : 1 : N-1);  % Fs/N就是这个频谱中的最小频率间隔！！！！！所以N越大，分辨率会越高
p_ = fft(p);
% 平移频域信号
fshift = (-N/2:N/2-1)*(fs/N);  % 平移后信号的频域范围
p_ = fftshift(p_);  % fftshift将零频分量移动到数组中心，重新排列
% amp1 = abs(p_) * 2 / N ;
amp1 = abs(p_);
plot(amp1);
title("平移后频域信号（未更改频域范围）");
% plot(f(1:N/2), amp1(1:N/2));  % fft算法默认是双边谱,通常我们只取一半

%% 取出谐波成分，将平移后的频谱置零取出(丽萍学姐论文)
% pp_ = takeHarmonicComponent(p_,6000,8000);  % 这个范围还是不能太小嗷，6500-7500误差就变大了
pp_ = takeHarmonicComponent(p_,25500,26500);
% pp_ = takeHarmonicComponent(p_,25470,26480);
subplot(9,1,3);
amp2 = abs(pp_) * 2 / N ;
plot(fshift,amp2);
title("取出的一次谐波成分（更改了频域范围）");

subplot(9,1,4);
p1 = ifft(ifftshift(pp_));
p1 = imag(p1);  % 呃。。。文章里一次谐波取的是虚部，这里为啥一次谐波取实部啊，所以这是sin
plot(p1);
title("一次谐波时域信号")

% pp_ = takeHarmonicComponent(p_,8000,10000);
pp_ = takeHarmonicComponent(p_,26500,27500);
subplot(9,1,5);
p2 = ifft(ifftshift(pp_));
p2 = real(p2);
plot(p2);
title("二次谐波时域信号")

%% 条纹翻倍
figure();
N = 13610;
p1 = p1(1:N);  % 太长了，取短一点来用
p2 = p2(1:N);

% if besselj(2,2*h)==0 && besselj(1,2*h)==0
%     error('贝塞尔函数0值处，无意义')
% end
% 一类besselj、二类bessely、三类besselh
h = 50;
p1 = p1 ./ -besselj(1,2*h);  % 这里要删除掉对应的贝塞尔函数值吗
p2 = p2 ./ besselj(2,2*h);
subplot(7,1,1);
n = n;  % 翻倍次数
times = 2^n;
[phiF_wrapped,Pn,Pnn] = SMI_API_TANDOUBLE_specific(n,p1,p2);
phiF_wrapped_init = phiF_wrapped; 
% padding = N/100;
% phiF_wrapped = wextend('1','sym',phiF_wrapped,padding);  % 因为数条纹要插值，插值在信号开始和结尾处误差有问题
plot(phiF_wrapped);
% hold on;
% plot(p_init);


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

[C_reconstruct, alpha_reconstruct] = SMI_API_ESTIMATE_C(phiF_reconstruct);  
subplot(7,1,3);
phi0_reconstrut = phiF_reconstruct + C_reconstruct*sin(phiF_reconstruct+atan(alpha_reconstruct));
Lt_reconstruct = phi0_reconstrut * lambda / (4 * pi);
% Lt_reconstruct(1:padding)=[];
% Lt_reconstruct(end-padding+1:end)=[];
% plot(Lt,'k');
% hold on;
% Lt_reconstruct = Lt_reconstruct - mean(Lt_reconstruct); % 简写振动和余弦调制振动，重构后加上幅值A
% Lt_reconstruct = Lt_reconstruct + 1.5 * lambda;  % 重构后的随机振动信号要加上幅值1.5的波长，这是为啥我页不知道
Lt_reconstruct = sgolayfilt(Lt_reconstruct,1,11);
Lt_reconstruct = sgolayfilt(Lt_reconstruct,2,21);
Lt_reconstruct = sgolayfilt(Lt_reconstruct,3,31);
plot(Lt_reconstruct,'r')
Lt_reconstruct1 = Lt_reconstruct;
title(['解包裹重构后的信号，C-reconstruct=', num2str(C_reconstruct)]);

%% 误差分析
% subplot(7,1,4);
% plot(Lt-Lt_reconstruct);
% RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
% title(['绝对误差，RMSE=', num2str(RMSE)])

%% 数条纹重构
% 找到所有的峰谷值（理论上包含跳变点）
% phiF_wrapped_init = phiF_wrapped_init(1:20000);
% N = 20000;
subplot(7,1,5);
% phiF_wrapped_init = sgolayfilt(phiF_wrapped_init,1,11);  % 第二个参数是阶数，通常123就够啦。第三个参数是长度，长度越短，滤除的就是越高频的分量，因为噪声是高频分量，所以一般不要太大！！！
% phiF_wrapped_init = sgolayfilt(phiF_wrapped_init,2,21);
% phiF_wrapped_init = sgolayfilt(phiF_wrapped_init,3,31);
[top_p,loc_p] = findpeaks(phiF_wrapped_init,'minpeakheight',0.8);  % 找到所有大于0.1的峰值,这是因为驼峰区乱七八糟（当然这只针对该实验信号）
[top_v,loc_v] = findpeaks(-phiF_wrapped_init,'minpeakheight',0.8);  % 默认。 'minpeakheight',0.1, 'minpeakdistance',50
top_v = -top_v;
plot(phiF_wrapped_init);
hold on;
% scatter(loc_p,top_p);
% scatter(loc_v,top_v);

%% 当findpeaks无法精确定位实验信号的所有峰谷值时（峰谷值中包含翻转点，亦或是有一些异常点）,针对这些异常点的特殊处理
arr = [9583,12110,12187,9452,4529,4654,4758,2090,2049,2191,4560,4723,9428,9525,9641,14549,14693,19588,19677];
arr1 = [2111,4676,7093,9577,12144];
[loc_p,loc_v,top_p,top_v] = ditchReversePoing(arr,loc_p,loc_v,top_p,top_v);
scatter(loc_p,top_p);
scatter(loc_v,top_v);
hold on;
%% 确定方向3：手动确定方向（最方便）
dir = solve_dir(arr1,N);  % 输入翻转点arr数组(顺序输入) （方向有误乘-1即可）
direction = -dir;
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
        Lt_reconstruct_beforeinternp(i) = Lt_reconstruct_beforeinternp(i-1) + step*direction(trigger(i));
    elseif direction(trigger(i))~=direction(trigger(i-1))
        Lt_reconstruct_beforeinternp(i) = Lt_reconstruct_beforeinternp(i-1);
    end
end
% plot(Lt_reconstruct_beforeinternp);
% title("Lt-before-internp")

Lt_reconstruct = interp1(trigger, Lt_reconstruct_beforeinternp, 1:N, "spline");  % 这个插值，在开始后结束位置很容易出现大误差
% Lt_reconstruct = Lt_reconstruct - mean(Lt_reconstruct);
% plot(Lt);
% hold on;
plot(Lt_reconstruct);
hold on;
plot(Lt_reconstruct1);
RMSE = sqrt(mean((Lt_reconstruct1-Lt_reconstruct).^2));
title(['Lt(count fringe)&均方误差:RMSE=', num2str(RMSE)]);

%% 
subplot(7,1,7);
padding = N/450;
plot(Lt_reconstruct1(100:13400)-Lt_reconstruct(100:13400));
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
%% 取出谐波信号函数
function pp_ = takeHarmonicComponent(p_,X,Y)  % p_为输入频谱。其中XY分别为输入频谱的起始点和终点。N为采样点数
     N = length(p_);
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