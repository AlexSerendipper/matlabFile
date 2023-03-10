%% 该程序用于调试（实验信号有两种思路：）
% 一是求出完美的峰谷值（不包括翻转点），利用包络求出的方向（所以结果并不完全准确）直接进行重构
% 二是求出完美的峰谷值（包括翻转点），求翻转点可以利用平均峰值距离，但是并不是所有的信号都能求出翻转点
% 总的来说，实验信号并不具备普适性，所以需要针对不同的实验信号进行不同的处理，所以直接用方法1，配合自定义方向和自定义翻转点是个不错的方案
%% 这里主要就是对偶次幂是否能提高精度查看验证
%% 全局变量
clc;
clear all;
close all;
C = 0.2;
fs = 200000;  % 采样率，即fs(s)采一个点。
N = 4000;  
fv = 100;  % 震动频
alpha = 5;

%% 实验信号
figure(1);
subplot(7,1,1);

% path =  'D:\matlab save\smi_实验信号\L2_weak_sanban_2240481_2255480';  % 1 文件路径
% M = 2240481; N = 15000;  [t, p, fs] = MOVE_API_EXPERIMENT(M, N, path);  % 1 加载.csv文件，从M点处开始取N个点
load('D:\matlab save\self-mixing\smi_MODULATE\temp_experi._p2.mat');  % 2. 加载.mat文件
p = p2;
lambda = 650e-9;  % 波长
p = sgolayfilt(p,1,11);
p = sgolayfilt(p,2,21);
p = sgolayfilt(p,3,31);

% SMI_API_CORR_FILTER(p,smoothingfactor,threshold) 使用自相关能够有效去噪，但是在后续处理方向上反而不太好，可能是信号带散斑的原因
N = 13610;
p = p(1:N);

p = -1 + (p - min(p))/(max(p) - min(p)) * 2;  % 还是需要归一化一下，否则无法求acosp
plot(p);
hold on;

%% 找到所有的峰谷值（理论上包含跳变点）
% [top_p,loc_p] = findpeaks(p,'minpeakheight',0.1, 'minpeakdistance',50);  % 找到所有大于0.1的峰值,这是因为驼峰区乱七八糟（当然这只针对该实验信号）
[top_p,loc_p] = findpeaks(p, 'minpeakdistance',500);  % 找到所有大于0.1的峰值,这是因为驼峰区乱七八糟（当然这只针对该实验信号）
[top_v,loc_v] = findpeaks(-p, 'minpeakdistance',500);  % 默认。 'minpeakheight',0.1, 'minpeakdistance',50
top_v = -top_v;
% scatter(loc_p,top_p);
% scatter(loc_v,top_v);


%% 当findpeaks无法精确定位实验信号的所有峰谷值时（峰谷值中包含翻转点，亦或是有一些异常点）,针对这些异常点的特殊处理
% arr = [754,3266,5731,5848,8177,8339,10808,10881,13273,13192,13398,15741,15664,18241,18291];
% arr = [2071,2201,4649,4719,7078,9572,12140,12178,14639,14685,17159,19675];
arr = [2191,4655,7091,9524,11943,12222];
[loc_p,loc_v,top_p,top_v] = ditchReversePoing(arr,loc_p,loc_v,top_p,top_v);

% temp_posi = find(loc_p==1235);
% loc_p(temp_posi)=[];
% top_p(temp_posi)=[];  % 他们对应的位置相同
% temp_posi = find(loc_v==6242);
% loc_v(temp_posi)=[];
% top_v(temp_posi)=[];  % 他们对应的位置相同

scatter(loc_p,top_p);
scatter(loc_v,top_v);

%% 确定方向1：手动确定方向（最方便）yyds
subplot(7,1,2);
% arr = [754,3266,5731,8177,10808,13273,15741,18291];
% arr = [2143,4649,7078,9572,12140,14639,17159,19675];
arr = [2191,4655,7091,9524,11943];
dir = solve_dir(arr,N);  % 输入翻转点arr数组(顺序输入) （方向有误乘-1即可）
plot(p);
hold on;
plot(dir)
direction = dir;


%% part2: 重构
subplot(7, 1, 3);
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

for i = 2:N  
    if ismember(i, loc_v) && (direction(i)==-1)
        add_op(i:end) = add_op(i:end) -  2 * pi;
    elseif ismember(i, loc_v) && (direction(i)==1)
        add_op(i:end) = add_op(i:end) +  2 * pi;
    end
end

% add_op = add_op - (max(add_op) + min(add_op))/2;
phiF_reconstruct = acosp_op1 + add_op;

%% 参数估算，这里可以不开，实验信号就必须开
[C_reconstruct, alpha_reconstruct] = SMI_API_ESTIMATE_C(phiF_reconstruct);  
% [C_reconstruct, alpha_reconstruct] = SMI_API_ESTIMATE_ALPHA(direction,loc_v,loc_p,top_v,top_p);  


%% 重构
subplot(7,1,4);
phi0_reconstrut = phiF_reconstruct + C_reconstruct*sin(phiF_reconstruct+atan(alpha_reconstruct));  % 这里的alpha如果用估算的，就会引入蛮大的误差
Lt_reconstruct = phi0_reconstrut * lambda / (4 * pi);


Lt_reconstruct = Lt_reconstruct - mean(Lt_reconstruct); % 简写振动和余弦调制振动，重构后加上幅值A
% Lt_reconstruct = Lt_reconstruct + 1.5 * lambda;  % 重构后的随机振动信号要加上幅值1.5的波长，这是为啥我页不知道
plot(Lt_reconstruct,'r')
hold on;
title("重构后的信号");

%% 实验信号通常需要对其平滑化处理
subplot(7,1,5);
Lt_reconstruct = sgolayfilt(Lt_reconstruct,1,11);  % 第二个参数是阶数，通常123就够啦。第三个参数是长度，长度越短，滤除的就是越高频的分量，因为噪声是高频分量，所以一般不要太大！！！
Lt_reconstruct = sgolayfilt(Lt_reconstruct,2,21);
Lt_reconstruct = sgolayfilt(Lt_reconstruct,3,31);
plot(Lt_reconstruct,'r')
hold on;
title("重构后的信号-平滑处理后");

%% 误差分析(仿真信号才有)
% subplot(7,1,5);
% plot(Lt);
% hold on;
% plot(Lt_reconstruct);
% title("LtandLt_reconstruct")
% subplot(7,1,6);
% plot(Lt-Lt_reconstruct)
% RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
% title(['绝对误差，RMSE=', num2str(RMSE)])


















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