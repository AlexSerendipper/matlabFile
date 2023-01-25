%% 该程序用于调试，其实也是非常好用且重要的
%% 这里不讨论去直流，换一种重构方式PUM看看

%% 全局变量
clc;
clear all;
close all;
fs = 200000;  % 采样率，即fs(s)采一个点。
N = 4000;  
fv = 100;  % 震动频
alpha = 4;

%% 产生自混合信号
figure(1);
subplot(8,1,1);
t = (0:N-1)/fs;  % 采样时间，设N=10, fs=200，即采样了0.05s，t为[0...0.045]
lambda = 650e-9;  % 波长
A = 2 * lambda;  % 幅值（✔）
L0 = 20 * lambda;  % 外腔距离（✔） 
Lt = A.* sin(2*pi*fv*t);  % L0为标称位置，外腔长度
beta = 1;  % the amplitude of selfmixing signal
phi0 = 4*pi*(L0+Lt)/lambda;

p = zeros(1,N);
C = [2.5,2.5];
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
% p = p .* (1+0.3*cos(2*pi*100*t));  % 给自混合信号加包络，加了一个幅值为0.2，频率为75的包络
plot(p);
p_init = p;
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
subplot(8,1,2);
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
subplot(8,1,3);
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

subplot(8,1,4);
plot(diffp);
hold on;
scatter(loc_diffp_seg,0,"g");
title("diffp，在direction>0的时候求极小值！！！，direction<0的时候求极大值！！！这样求出的区间(绿点)内不再包含条纹~~~")

%% 取出翻转点，与上述方法相同！
subplot(8,1,5);
top_r_temp = [];
loc_r_temp = [];
for i = 1:2:length(loc_diffp_seg)
    [temp1,temp2] = findpeaks(p(loc_diffp_seg(i):loc_diffp_seg(i+1)));
    [temp3,temp4] = findpeaks(-p(loc_diffp_seg(i):loc_diffp_seg(i+1)));
    top_r_temp = [top_r_temp, temp1, -temp3];
    loc_r_temp = [loc_r_temp, temp2 + loc_diffp_seg(i) - 1, temp4 + loc_diffp_seg(i) - 1 ];
end
plot(p);
hold on;
scatter(loc_r_temp,top_r_temp);
title("在上述区间内求出翻转点")

%% 挖去峰值中的跳变点，返回，峰值、谷值、跳变点！
for i = 1:length(loc_p)
    for j = 1:length(loc_r_temp)
        if loc_p(i) == loc_r_temp(j)
            loc_p(i) = nan;
            top_p(i) = nan;
        end
    end
end

% 挖去谷值中的跳变点
for i = 1:length(loc_v)
    for j = 1:length(loc_r_temp)
        if loc_v(i) == loc_r_temp(j)
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
subplot(8,1,6);
direction = zeros(1,N);
direction(1) = dir(1);
k = 1;
for i = 1:N
    direction(i) = k;
    if ismember(i, loc_r_temp)
        k = k * -1;
    end
end

plot(p);
hold on;
scatter(loc_p,top_p);
scatter(loc_v,top_v);
scatter(loc_r_temp,top_r_temp);
plot(direction);
title("得到纯净的峰谷值和完全正确的方向")









%% 去直流分量
figure(1);
subplot(8,1,7);
plot(p);
hold on;
% 0原始做法，现在用原始做法配合时域去直流试试
% direction_seg1 = [];
% for i = 1:length(direction)-1
%     if direction(i) ~=  direction(i+1)
%         direction_seg1 = [direction_seg1, i];
%     end
% end
% direction_seg2 = direction_seg1 + 1;
% direction_seg = [1, sort([direction_seg1,direction_seg2]), length(direction)];  % 其中的1-2，3-4 为方向恒定的区间!!!!!!(ˇ∀ˇ)

% 1我要所有段，全部去直流(不行，有问题)
% loc_diffp_seg2 = loc_diffp_seg + 1;
% direction_seg = [1, sort([loc_diffp_seg,loc_diffp_seg2]), length(direction)];  % 其中的1-2，3-4 为方向恒定的区间!!!!!!(ˇ∀ˇ)

% 2我要所有段（除驼峰区），全部去直流
direction_seg = [1, loc_diffp_seg, length(direction)];  % 其中的1-2，3-4 为方向恒定的区间!!!!!!(ˇ∀ˇ)

% 3以导数信号作为条纹判断的依据，每个条纹区间都减去其平均值
% [top_diffp_p,loc_diffp_p] = findpeaks(diffp,'minpeakheight',0.2);  % 以导数信号作为条纹判断的依据，每个条纹区间都减去其平均值
% [top_diffp_v,loc_diffp_v] = findpeaks(-diffp,'minpeakheight',0.2);
% scatter(loc_diffp_p,1,'r');
% scatter(loc_diffp_v,-1,'b');
% direction_seg =[1, sort([loc_diffp_p,loc_diffp_p+1,loc_diffp_v,loc_diffp_v+1]), length(direction)];

% 4上条纹加驼峰区去直流
% loc_diffp_seg_1 = [];
% for i = 2:4:length(loc_diffp_seg)-1
%     loc_diffp_seg_1 = [loc_diffp_seg_1, loc_diffp_seg(i), loc_diffp_seg(i+1)];
% end
% direction_seg = [1, sort([loc_diffp_seg_1,loc_diffp_seg_1+1]), length(direction)];  % 其中的1-2，3-4 为我想要的恒定的区间!!!!!!(ˇ∀ˇ)

% 5





for i=1:2:length(direction_seg)-1  % 这里确实是1：2
    ave = mean(p(direction_seg(i):direction_seg(i+1)));
    p((direction_seg(i):direction_seg(i+1))) = p((direction_seg(i):direction_seg(i+1))) - ave;
end

plot(p);
hold on;

% 这样去直流后会有不连续点，现在试图通过插值来解决
% for i = 1:length(p)
%     for j = 1:length(loc_diffp_seg)
%         if (i == loc_diffp_seg(j))
%             p(i) = nan;
%         end
%     end
% end
% p = interp1(p,1:N,'spline');
% subplot(8,1,8);
% plot(p)








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
p = wextend('1D','sym',p,padding);  % 数据延拓
% p = [zeros(1,padding) p zeros(1,padding)];

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


%% 拉平时频谱去除包络
% TF_curb(find(TF_curb~=0)) = TF_curb(find(TF_curb~=0))./ abs(TF_curb(find(TF_curb~=0))) * 55; % 找到其不为0处的索引,把这些处的值设置除其各自的模值，这样这些复数的模长都为1，只不过相角不同
subplot(7,2,[10,12,14]);
mesh(abs(TF_curb)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' );
% mesh(angle(TF_curb));
% view(0,90);
title('抑制后');


%% 
p = (istft(TF_curb, fs, 'Window', window, 'OverlapLength', overlapLength, 'FFTLength', 5*windowLength))';
p_init2 = p;  % 用p_init2存储未归一化的信号
p_init2 = p_init2(padding+1:end-padding);  % ✔去掉之前的padding
p_init2 = sgolayfilt(p_init2,2,31);  % ✔ 超参数，平滑滤波器, 这里的参数也需要细调！！！c小的时候可以设置4，31   c大的时候可以设置，2，21 

figure(2);
subplot(7,2,1);
plot(p_init2);
title('未归一化');
%% 方向信息
direction = dir;  % 目前使用dir的方向信息

%% PUM变换重构
figure(2);
subplot(7, 2, 3);
p = p_init2;
plot(p);
hold on;
[top_p,loc_p] = findpeaks(p);
% loc_p_convert = N2T(loc_p,N,T);
[top_v,loc_v] = findpeaks(-p);
top_v = -top_v;
% loc_v_convert = N2T(loc_v,N,T);
scatter(loc_p,top_p);
scatter(loc_v,top_v);
plot(dir);

%%
subplot(7, 2, 5);
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

top_diffp_seg([1,end]) = [];
loc_diffp_seg([1,end]) = [];  % 最适合找翻转点的区间
plot(p);
hold on;
scatter(direction_seg,0,"r");
title("根据方向dir，找出用来求diffp极值的区间（红点）（恒定正负1）")

subplot(7,2,7);
plot(diffp);
hold on;
scatter(loc_diffp_seg,0,"g");
title("diffp，在direction>0的时候求极大值！！！，direction<0的时候求极小值！！！这样求出的区间(绿点)内不再包含条纹~~~")

%% 取出翻转点，与上述方法相同！
subplot(7,2,9);
top_r_temp = [];
loc_r_temp = [];
for i = 1:2:length(loc_diffp_seg)
    [temp1,temp2] = findpeaks(p(loc_diffp_seg(i):loc_diffp_seg(i+1)));
    [temp3,temp4] = findpeaks(-p(loc_diffp_seg(i):loc_diffp_seg(i+1)));
    top_r_temp = [top_r_temp, temp1, -temp3];
    loc_r_temp = [loc_r_temp, temp2 + loc_diffp_seg(i) - 1, temp4 + loc_diffp_seg(i) - 1 ];
end
plot(p);
hold on;
loc_r = [];  % 这里先根据规律取出来看看效果 !!!!!!!!!!!!!!!!!!!!
top_r = [];
loc_r(1) = loc_r_temp(1);
loc_r(2) = loc_r_temp(6);
loc_r(3) = loc_r_temp(7);
loc_r(4) = loc_r_temp(12);

top_r(1) = top_r_temp(1);
top_r(2) = top_r_temp(6);
top_r(3) = top_r_temp(7);
top_r(4) = top_r_temp(12);
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










%% 重构
subplot(7, 2, 13);
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

%% 参数估算
[C_reconstruct, alpha_reconstruct] = SMI_API_ESTIMATE_C(phiF_reconstruct);  
% [C_reconstruct, alpha_reconstruct] = SMI_API_ESTIMATE_ALPHA(direction,loc_v,loc_p,top_v,top_p);  


%% 重构
figure(3);
subplot(5,1,1);
phi0_reconstrut = phiF_reconstruct + C_reconstruct*sin(phiF_reconstruct+atan(alpha_reconstruct));  % 这里的alpha如果用估算的，就会引入蛮大的误差
Lt_reconstruct = phi0_reconstrut * lambda / (4 * pi);


Lt_reconstruct = Lt_reconstruct - mean(Lt_reconstruct); % 简写振动和余弦调制振动，重构后加上幅值A
% Lt_reconstruct = Lt_reconstruct + 1.5 * lambda;  % 重构后的随机振动信号要加上幅值1.5的波长，这是为啥我页不知道
plot(Lt,'k')
hold on;

plot(Lt_reconstruct,'r')
title("重构后的信号");

%% 误差分析
subplot(5,1,2);
plot(Lt-Lt_reconstruct)
RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
title(['绝对误差，RMSE=', num2str(RMSE)])


















 
%% just for test 
%% 可视化区域
figure(4);
subplot(5,1,1);
plot(Lt);
title('Lt')
subplot(5,1,2);
plot(p);
title('p')
subplot(5,1,3);
plot(phiF_reconstruct);
title('phiF-reconstruct')
subplot(5,1,4);
plot(phi0_reconstrut);
title('phi0-reconstrut')
subplot(5,1,5);
% plot(Lt_reconstruct);
% title('Lt-reconstruct')
plot(c)



























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


%% subfunction （eliminate_dc,输入一段信号p，及其长度N ，输出去除直流分量后的信号p_）
function g_ = eliminate_dc(p_)
    p_= fft(p_);
    p_(1) = 0;
    g_ = ifft(p_);
end