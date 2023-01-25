%% 分段去直流，因为当反馈强度逐渐增大，且目标物远离激光器运动时会在信号中引入较大的直流成分，在时域上体现为自混合信号条纹向上偏移。

%% 全局变量
clc;
clear all;
close all;

fs = 200000;  % 采样率，即fs(s)采一个点。
N = 4000;  
fv = 200;  % 震动频
C = [0.5,0.5];
alpha = 4.6;

%% 全局变量
windowLength = 512; % 窗长
overlapLength = floor(windowLength * 0.9);  % OverlapLength后为指定的重叠长度
% overlapLength = floor(50);  % OverlapLength后为指定的重叠长度
window = hamming(windowLength, "periodic");  % 使用汉明窗作为滑动的窗口
% window = kaiser(windowLength,5);
fftLength = 5*windowLength;
% fftLength = 512;
M = C(1)/10; % 抑制因子
V = 0.9; % 抑制因子2

% %% 产生自混合信号
% figure(1);
% subplot(7,2,1);
% t = (0:N-1)/fs;  % 采样时间，设N=10, fs=200，即采样了0.05s，t为[0...0.045]
% lambda = 650e-9;  % 波长
% A = 2 * lambda;  % 幅值（✔）
% % A = 3e-6;  % 幅值（✔）
% L0 = 20 * lambda;  % 外腔距离（✔） 
% Lt = A.* sin(2*pi*fv*t);  % L0为标称位置，外腔长度
% beta = 1;  % the amplitude of selfmixing signal
% phi0 = 4*pi*(L0+Lt)/lambda;
% 
% p = zeros(1,N);
% 
% 
%  % C的变化是一个正弦曲线，不能随机数！
% C_upper = C(1);
% C_top = C(2);
% % 这个乘和加保证了c的上下限
% x = linspace(0, 3*pi, N);
% c = (C(2)-C(1))/2 * cos(x) + (C_top - (C(2)-C(1))/2);
% plot(1:N,c);ylabel('Amp.(a.u)','FontSize',7,'FontWeight','bold'); xlabel('sample(num)','FontSize',7,'FontWeight','bold');
% for i = 1:N 
%     CC = c(i);
%     p(i) = beta * cos(solve_phiF(CC, phi0(i), alpha));  % 遍历所有的phi0
% end
% % p = p + randn(size(p));
% % p = awgn(p,10); % 10db，加噪声在抑制是可以的，问题是加了噪声没法从原信号判断方向了。
% plot(p);
% p_init = p;
% title(['自混合信号,C=',num2str(C)]);

%% 自混合实验信号
figure(1);
subplot(7,2,1);
path =  'D:\matlab save\self-mixing\smi_实验信号\gexiong_f(100)_A(3um).csv';  % 5 文件路径
M = 112000; N = 16000;  [t, p, fs] = MOVE_API_EXPERIMENT(M, N, path);  % 5 从M点处取N个点
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
figure(1);
subplot(7,2,3);
plot(p);
title('未归一化');

% subplot(7,2,5);
% p = 2 * (p - min(p))./(max(p) - min(p)) - 1;
% p = p(padding+1:end-padding);  % 去掉之前的padding
% plot(p);
% title('归一化1');

subplot(7,2,7);
p = 2 * (p - -abs(hilbert(-p)))./(abs(hilbert(p)) - -abs(hilbert(-p))) - 1;
% p = p./abs(p);  % ×这个有问题，我自己乱推的
p = p(padding+1:end-padding);  % 去掉之前的padding
p = sgolayfilt(p,4,31);  % ✔平滑滤波器，这里的参数也需要细调！！！第一个参数好像是重复的次数，第二个是强度，重复次数必须小于强度，强度必须为奇数          c小的时候可以设置4，31   c大的时候可以设置2，21
plot(p);
title('归一化2');



%% 找到所有的峰谷值（包含跳变点）
hold on;
[top_p,loc_p] = findpeaks(p);
% loc_p_convert = N2T(loc_p,N,T);
[top_v,loc_v] = findpeaks(-p);
top_v = -top_v;
% loc_v_convert = N2T(loc_v,N,T);
scatter(loc_p,top_p);
scatter(loc_v,top_v);

%% 包络确定方向！
% subplot(7,2,9);
diffp = diff(p);  % diff是相邻两个数的差分，对第一个位置补0
diffp = [0,diffp];
diff_acosp = diff(acos(p));
diff_acosp = [0,diff_acosp];
plot(diffp);
hold on;

[top_diffp_p,loc_diffp_p] = findpeaks(diffp);  % 拿到极值和索引值
[top_diffp_v,loc_diffp_v] = findpeaks(-diffp);
% scatter(loc_diffp_p,top_diffp_p);
% scatter(loc_diffp_v,-top_diffp_v);
en_top = interp1(loc_diffp_p,top_diffp_p,1:N,'spline');  % 三次样条插值，曲线更平滑
en_bottom = interp1(loc_diffp_v,-top_diffp_v,1:N,'spline');

en_median = (en_top - (en_top - en_bottom)/2);
dir = -sign(en_median);  % 让dir暂时指明方向

% plot(en_top);
% plot(en_bottom);
% plot(dir);
% title("diffp及其包络确定方向")

%% 利用平均峰值，求出翻转点，因为翻转点都是M型的
Avg_p = sum(diff(loc_p)) / (length(loc_p)) + 75;  % 平均峰值距离，这个length(loc_p)可以调...这里实验信号算出来平均峰值太小了，A= 1.3刚好，3um我直接补了个70...我醉了
Avg_v = sum(diff(loc_v)) / (length(loc_p)) + 150;  % 平均峰值距离，这个length(loc_p)可以调...这里实验信号算出来平均峰值太小了，A= 1.3刚好，3um我直接补了个70...我醉了

loc_r = [];
top_r = [];
for i = length(loc_p)-1:-1:2
    % 求峰值的翻转点
    if abs(loc_p(i) - loc_p(i-1)) > Avg_p && abs(loc_p(i) - loc_p(i+1)) > Avg_p % 每个点都和左右边的点比，距离都大于平均峰值就是翻转点
        loc_r = [loc_r,loc_p(i)];
        top_r = [top_r, top_p(i)];     
    end
    % 求谷值的翻转点
    if abs(loc_v(i) - loc_v(i-1)) > Avg_v && abs(loc_v(i) - loc_v(i+1)) > Avg_v % 每个点都和左右边的点比，距离都大于平均峰值就是翻转点
        loc_r = [loc_r,loc_v(i)];
        top_r = [top_r, top_v(i)];     
    end
    
end

subplot(7,2,9);
plot(p);
hold on;
scatter(loc_r,top_r,"g");


%% 挖去峰值中的跳变点，返回，峰值、谷值、跳变点！
for i = 1:length(loc_p)
    for j = 1:length(loc_r)
        if loc_p(i) == loc_r(j)
            loc_p(i) = nan;
            top_p(i) = nan;
        end
    end
end

% 挖去谷值中的跳变点
for i = 1:length(loc_v)
    for j = 1:length(loc_r)
        if loc_v(i) == loc_r(j)
            loc_v(i) = nan;
            top_v(i) = nan;
        end
    end
end

% 删除数组中的nan
loc_p(isnan(loc_p))=[];
loc_v(isnan(loc_v))=[];
top_p(isnan(top_p))=[];
top_v(isnan(top_v))=[];


%% 根据求出的翻转点修正一下dir信息！
subplot(7,2,11);
direction = zeros(1,N);
direction(1) = dir(1);
k = 1;
for i = 1:N
    direction(i) = k;
    if ismember(i, loc_r)
        k = k * -1;
    end
end

plot(p);
hold on;
scatter(loc_p,top_p);
scatter(loc_v,top_v);
scatter(loc_r,top_r);
plot(direction);
title("得到纯净的峰谷值和完全正确的方向")
 
%% just for test 




























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