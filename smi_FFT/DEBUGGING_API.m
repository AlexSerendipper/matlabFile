%% 该程序用于调试，其实也是非常好用且重要的
% zabit2021 IEEEsensorjournal，频域位移重构仿真

%% 全局变量
clc;
clear all;
close all;
fs = 200000;  % 采样率，即fs(s)采一个点。
N = 4000;  
fv = 200;  % 震动频
alpha = 6;

%% 产生自混合信号
figure(1);
subplot(7,1,1);
t = (0:N-1)/fs;  % 采样时间，设N=10, fs=200，即采样了0.05s，t为[0...0.045]
lambda = 650e-9;  % 波长
A = 2 * lambda;  % 幅值（✔）
L0 = 20 * lambda;  % 外腔距离（✔） 
Lt = A.* sin(2*pi*fv*t);  % L0为标称位置，外腔长度
beta = 1;  % the amplitude of selfmixing signal
phi0 = 4*pi*(L0+Lt)/lambda;

p = zeros(1,N);
C = [2, 2];
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

%% 找到所有的峰谷值（包含跳变点）
[top_p,loc_p] = findpeaks(p);
% loc_p_convert = N2T(loc_p,N,T);
[top_v,loc_v] = findpeaks(-p);
top_v = -top_v;
% loc_v_convert = N2T(loc_v,N,T);
scatter(loc_p,top_p);
scatter(loc_v,top_v);

%% 包络确定方向！
subplot(7,1,2);
diffp = diff(p);  % diff是相邻两个数的差分，对第一个位置补0
diffp = [0,diffp];
diff_acosp = diff(acos(p));
diff_acosp = [0,diff_acosp];

plot(diffp);
hold on;
[top_diffp_p,loc_diffp_p] = findpeaks(diffp);  % 拿到极值和索引值
[top_diffp_v,loc_diffp_v] = findpeaks(-diffp);
scatter(loc_diffp_p,top_diffp_p);
scatter(loc_diffp_v,-top_diffp_v);
en_top = interp1(loc_diffp_p,top_diffp_p,1:N,'spline');  % 三次样条插值，曲线更平滑
en_bottom = interp1(loc_diffp_v,-top_diffp_v,1:N,'spline');

en_median = (en_top - (en_top - en_bottom)/2);
dir = -sign(en_median);  % 让dir暂时指明方向

plot(en_top);
plot(en_bottom);
plot(dir);
title("diffp及其包络确定方向")


%% 利用方向信息，求出最合适的找翻转点的范围，大大增加了鲁棒性
subplot(7,1,3);
direction_seg1 = [];  % 方向发生变化的点(ˇ∀ˇ)
for i = 1:length(dir)-1
    if dir(i) ~=  dir(i+1)
        direction_seg1 = [direction_seg1, i];
    end
end

direction_seg2 = direction_seg1 + 1;
direction_seg = [1, sort([direction_seg1,direction_seg2]), length(dir)];  % 其中的1-2，3-4 为方向恒定的区间!!!!!!(ˇ∀ˇ)
top_diffp_seg = [];
loc_diffp_seg = [];  % 存储翻转点的区间

for i = 1 : 2 : length(direction_seg)
    if dir(direction_seg(i):direction_seg(i+1)) < 0  % 在方向小于0的时候求diffp【两端】的谷值！！！！这样找翻转点就不用缩小范围了！！！！！！！！！！！！！！！！！！！！
        [tem1, tem2] = findpeaks(-diffp(direction_seg(i):direction_seg(i+1)));
        top_diffp_seg = [top_diffp_seg, -tem1(1), -tem1(end)];
        loc_diffp_seg = [loc_diffp_seg, tem2(1)+ direction_seg(i) - 1 , tem2(end) + direction_seg(i) - 1];   % 这里返回的索引是个难点！由于我是分段进行峰值寻找，索引值应当加上初始值（假设在[x,y]中找极值，返回索引是2，第2是极值点,实际上对应的索引值为x+2-1）
        
    else
        [tem3, tem4] = findpeaks(diffp(direction_seg(i):direction_seg(i+1)));  % 在方向大于0的时候求diffp【两端】的峰值！！！！这样找翻转点就不用缩小范围了！！！！！！！！！！！！！！！！！！！！
        top_diffp_seg = [top_diffp_seg, tem3(1), tem3(end)];
        loc_diffp_seg = [loc_diffp_seg, tem4(1) + direction_seg(i) - 1, tem4(end) + direction_seg(i) - 1];
    end
end
% plot(diffp);
top_diffp_seg([1,end]) = [];
loc_diffp_seg([1,end]) = [];  % 最适合找翻转点的区间
plot(p);
hold on;
scatter(direction_seg,0,"r");
title("根据方向dir，找出用来求diffp极值的区间（红点）（恒定正负1）")

subplot(7,1,4);
plot(diffp);
hold on;
scatter(loc_diffp_seg,0,"g");
title("diffp，在direction>0的时候求极小值！！！，direction<0的时候求极大值！！！这样求出的区间(绿点)内不再包含条纹~~~")

%% 取出翻转点，与上述方法相同！
subplot(7,1,5);
top_r = [];
loc_r = [];
for i = 1:2:length(loc_diffp_seg)
    [temp1,temp2] = findpeaks(p(loc_diffp_seg(i):loc_diffp_seg(i+1)));
    [temp3,temp4] = findpeaks(-p(loc_diffp_seg(i):loc_diffp_seg(i+1)));
    top_r = [top_r, temp1, -temp3];
    loc_r = [loc_r, temp2 + loc_diffp_seg(i) - 1, temp4 + loc_diffp_seg(i) - 1 ];
end
plot(p);
hold on;
scatter(loc_r,top_r);
title("在上述区间内求出翻转点")

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
subplot(7,1,6);
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


%% 重构
figure(2);
subplot(7, 1, 1);
acosp = acos(p);
plot(acosp);
title("acosp");
hold on;

% 设定初值,固定
if sign(acosp(2) - acosp(1)) == direction(1)
    init = 1;
else
    init = -1;
end

% 遇到峰谷值乘-1
mul_op = init * ones(1,N);
for i = 1:N
    if ismember(i, [loc_v,loc_p]) == 1  % 当遇到翻转点，变一个方向（折叠）
        mul_op(i:end) = -mul_op(i:end);
    end
end

acosp_op1 =  acosp .* mul_op;

plot(acosp_op1);
hold on;
title("翻转点×-1");


add_op = zeros(1,N);  % 累加阶梯

% for i = 2:N
%     if ((acosp_op1(i)-acosp_op1(i-1)) > 1) && (direction(i)==-1)  % 在下降区碰到谷值减2pi
%         add_op(i:end) = add_op(i:end) -  2 * pi;
%     elseif ((acosp_op1(i)-acosp_op1(i-1)) < -1) && (direction(i)==1)  % 该上升区碰到谷值加2pi
%         add_op(i:end) = add_op(i:end) +  2 * pi;
%     end
% end

for i = 2:N  
    if ismember(i, loc_v) && (direction(i)==-1)
        add_op(i:end) = add_op(i:end) -  2 * pi;
    elseif ismember(i, loc_v) && (direction(i)==1)
        add_op(i:end) = add_op(i:end) +  2 * pi;
    end
end

% add_op = add_op - (max(add_op) + min(add_op))/2;
phiF_reconstruct = acosp_op1 + add_op;
subplot(7, 1, 2);
plot(phiF_reconstruct);


%% 频域处理
subplot(7, 1, 3);
% f = fs / N * (0 : 1 : N-1);  % Fs/N就是这个频谱中的最小频率间隔！！！！！所以N越大，分辨率会越高
% phiF_reconstruct_ = fft(phiF_reconstruct, N);
% amp1 = abs(phiF_reconstruct_) * 2 / N;  % 
% plot(amp1);  % fft算法默认是双边谱,通常我们只取一半f(1:N/2), 
% title('phiF-reconstruct-amp');

%%
phiF_reconstruct_ = fft(phiF_reconstruct);
phiF_reconstruct_sin = fft(sin(phiF_reconstruct));
phiF_reconstruct_cos = fft(cos(phiF_reconstruct));

%%
fshift = (-N/2:N/2-1)*(fs/N);  % 平移后信号的频域范围
phiF_reconstruct_ = fftshift(phiF_reconstruct_);  % fftshift将零频分量移动到数组中心，重新排列
amp1 = abs(phiF_reconstruct_) * 2 / N ;
plot(amp1);
title('phiF-reconstruct-amp');

%%
phiF_reconstruct_sin_ = fftshift(phiF_reconstruct_sin);  % fftshift将零频分量移动到数组中心，重新排列
amp2 = abs(phiF_reconstruct_sin_) * 2 / N ;
phiF_reconstruct_cos_ = fftshift(phiF_reconstruct_cos);  % fftshift将零频分量移动到数组中心，重新排列
amp3 = abs(phiF_reconstruct_cos_) * 2 / N ;
i2 = amp2;
i3 = amp3;

%% try
% amp1(1997) = 0;
% amp1(2005) = 0;
% std1 = std(amp1);

% weight1 = amp1./max(amp1);
% weight1(find(weight1 == 1)) = weight1(find(weight1 == 1)) - 1e-11;  % 算法1：是找到weight1中等于1的位置（理论上每列都有），并都减去num极小值，因为这个后边要用做分母，所以减去一个无穷小量，不要让他产生无穷大
% weight2 = (weight1.^128)./(1 - weight1.^128);  % 算法1：构建一个函数，让原来频率分量大的更大，小的更小，windowLenght此处为抑制因子
% weight2 = weight2 ./ max(weight2);  % 抑制后按列归一化，用作权值矩阵
% amp2 = amp1 .* weight2;  % 获得抑制后的信号
% subplot(7, 1, 4);
% plot(amp2);
%% 
omega = sinh(amp1) ./ (max(amp1));
std1 = std(omega);
wb = find(omega>std1);
% ww = max(wb);
% wb(1) = 4001-ww;
% wb(2) = ww;
wb(1) = wb(1)-50;
wb(2) = wb(2)+50;
% 2000 - (ww-2000)

i1 = amp1;
i1(1:wb(1))=0;  % 我总觉得应该留下主频~~~
i1(wb(2):end)=0;
subplot(7, 1, 4);
% i1 = ifftshift(i1);
plot(i1);


i2(1:wb(1))=0;
i2(wb(2):end)=0;
% i2(2001) = 0;
subplot(7, 1, 5);
% plot(amp2);
% i2 = ifftshift(i2);
plot(i2)
 
i3(1:wb(1))=0;
i3(wb(2):end)=0;
% i3(2001)=0;
subplot(7, 1, 6);
% i3 = ifftshift(i3);
plot(i3)

f = @(x) norm(i1 + x(1).*i2 + x(2).*i3,1); 
[xmin,fmin] = fminsearch(f,[0,0]);




















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