%% 该程序用于调试，其实也是非常好用且重要的
%% 这里的想法是拿掉驼峰区，再给他补回去

%% 全局变量
clc;
clear all;
close all;
fs = 200000;  % 采样率，即fs(s)采一个点。
N = 4000;  
fv = 300;  % 震动频
alpha = 5;

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
C = [0.5,0.5];
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
p = p .* (1+0.5*cos(2*pi*100*t));  % 给自混合信号加包络，加了一个幅值为0.2，频率为75的包络
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
    if dir(direction_seg(i):direction_seg(i+1)) > 0  % 在方向小于0的时候求diffp【两端】的谷值！！！！这样找翻转点就不用缩小范围了！！！！！！！！！！！！！！！！！！！！
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












%% （去直流）
figure(2);
set(gcf,'position',[15, 30, 800, 800]);
subplot(7,1,7);
[p] = SMI_API_ELIMINATE_DC(p,dir);
p = sgolayfilt(p,1,21); 
p = sgolayfilt(p,2,21); 
p = sgolayfilt(p,3,21); 
plot(p); ylabel('p(t)(nm)','FontSize',7,'FontWeight','bold'); xlabel('sample(num)','FontSize',7,'FontWeight','bold');  % ✔这玩意赶紧学 
title(['自混合信号,C=',num2str(C), 'alpha=',num2str(alpha)]);

%% loc_diffp_seg就是我要设置的区间（驼峰区） —— 时域上去除驼峰区
% p_ = p;
% for i = 1:2:length(loc_diffp_seg)-1
%     p(loc_diffp_seg(i):loc_diffp_seg(i+1)) = 0;
% end
% subplot(7,1,7);
% plot(p);
%% 全局变量
windowLength = 128; % 窗长
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
p = [zeros(1,padding) p zeros(1,padding)];

%% 时频分析（抑制）
figure(2);
[TF,F,T] = stft(p, fs, 'Window', window, 'OverlapLength', overlapLength, 'FFTLength', fftLength);  % （TF每一列就是每一个时间维度的频率）
weight1 = abs(TF)./max(abs(TF));  % 要把这个矩阵当作权值使用，所以按列先归一化

% weight1(find(weight1 == 1)) = weight1(find(weight1 == 1)) - 1e-11;  % 算法1：是找到weight1中等于1的位置（理论上每列都有），并都减去num极小值，因为这个后边要用做分母，所以减去一个无穷小量，不要让他产生无穷大
% weight2 = (weight1.^64)./(1 - weight1.^64); % 算法1：构建一个函数，让原来频率分量大的更大，小的更小，windowLenght此处为抑制因子

weight2 = weight1;
weight2( find(weight2 < V) ) = 0; % 算法2：嘎嘎好用
weight2 = weight2 ./ max(weight2);  % 抑制后按列归一化
TF_curb = TF .* weight2;  % 获得抑制后的信号
subplot(7,2,[2,4,6]);
mesh(abs(TF)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' ); 
% view(0,90); % 设置初始视角为俯视角
% mesh(angle(TF));
title('抑制前');

%% 
gg = max(abs(TF));
hh = findpeaks(gg,'MinPeakDistance',length(T)/(length(loc_diffp_seg)/2 + 1 + 5)); % 这个峰值距离应该设置为横轴宽度/峰值个数+5, 小一点我想的是
jj = [];  % 拿到那些峰值所在的列！！！！！！！

for i =1:length(hh)
    jj = [jj, find(gg==hh(i))];
end

% 从这些峰值中左右各拿掉4列的数据,,,但是这样弄有用么...好像从TF_curb中拿掉比较合理
TF_curb_save = zeros(size(TF_curb));
for i = 1:length(jj)
    TF_curb_save(:,jj(i)-4:jj(i)+4) = TF_curb(:,jj(i)-4:jj(i)+4);
    TF_curb(:,jj(i)-4:jj(i)+4) = nan;
end

% 把驼峰区也拉平
TF_curb_save(find(TF_curb_save~=0)) = TF_curb_save(find(TF_curb_save~=0))./ abs(TF_curb_save(find(TF_curb_save~=0))) * 55; % 找到其不为0处的索引,把这些处的值设置除其各自的模值，这样这些复数的模长都为1，只不过相角不同

%% 把之前除去的驼峰区补回来
TF_curb(find(TF_curb~=0)) = TF_curb(find(TF_curb~=0))./ abs(TF_curb(find(TF_curb~=0))) * 55; % 找到其不为0处的索引,把这些处的值设置除其各自的模值，这样这些复数的模长都为1，只不过相角不同

for i = 1:length(jj)
    TF_curb(:,jj(i)-4:jj(i)+4) = TF_curb_save(:,jj(i)-4:jj(i)+4);
end

subplot(7,2,[10,12,14]);
mesh(abs(TF_curb)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' );
% mesh(angle(TF_curb));
% view(0,90);
title('抑制后');

%% 
p = (istft(TF_curb, fs, 'Window', window, 'OverlapLength', overlapLength, 'FFTLength', 5*windowLength))';
p_init2 = p;  % 用p_init2存储未归一化的信号
p_init2 = p_init2(padding+1:end-padding);  % ✔去掉之前的padding
p_init2 = sgolayfilt(p_init2,4,31);  % ✔ 超参数，平滑滤波器, 这里的参数也需要细调！！！c小的时候可以设置4，31   c大的时候可以设置，2，21 

figure(2);
subplot(7,2,3);
plot(p_init2);
title('未归一化');




















 
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