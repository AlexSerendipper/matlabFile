%% 该程序用于调试（实验信号有两种思路：）
% 一是求出完美的峰谷值（不包括翻转点），利用包络求出的方向（所以结果并不完全准确）直接进行重构
% 二是求出完美的峰谷值（包括翻转点），求翻转点可以利用平均峰值距离，但是并不是所有的信号都能求出翻转点
% 总的来说，实验信号并不具备普适性，所以需要针对不同的实验信号进行不同的处理，所以直接用方法1，配合自定义方向和自定义翻转点是个不错的方案

%% 全局变量
clc;
clear all;
close all;

%% 自混合信号
C = 2;
fs = 200000;  % 采样率，即fs(s)采一个点。
N = 4000;  
fv = 100;  % 震动频
alpha = 5;
lambda = 650e-9;  % 波长
% [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_HARMONIC(fs, N, fv, C, alpha);
% subplot(7,1,1);
% plot(p);
% hold on;

%% 实验信号
figure(1);
subplot(7,1,1);
% path =  'D:\matlab save\smi_实验信号\experSignal\1_w_smooth_365424_15000';  % 1 文件路径
path =  'D:\matlab save\smi_实验信号\experSignal_standard\4_w_110000_10000_f100_A3um.csv';  % 1 文件路径
M = 110000; N = 15000;  [t, p, fs] = MOVE_API_EXPERIMENT(M, N, path);  % 1 加载.csv文件，从M点处开始取N个点

% load('xxx.mat');  % 2. 加载.mat文件
p = sgolayfilt(p,1,11);
p = sgolayfilt(p,2,21);
p = sgolayfilt(p,3,31);
plot(p);
hold on;
arr = [1685,4265,6805,9243,11792,14314];
direction = SMI_API_DIR(arr,N);
plot(direction);
% SMI_API_CORR_FILTER(p,smoothingfactor,threshold) 使用自相关能够有效去噪，但是在后续处理方向上反而不太好，可能是信号带散斑的原因
% p = -1 + (p - min(p))/(max(p) - min(p)) * 2;  % 还是需要归一化一下，否则无法求acosp

%% 峰谷值1：找到所有的峰谷值（理论上包含跳变点）
% [top_p,loc_p] = findpeaks(p,'minpeakheight',0.1, 'minpeakdistance',50);  % 找到所有大于0.1的峰值,这是因为驼峰区乱七八糟（当然这只针对该实验信号）
[top_p,loc_p] = findpeaks(p);  % 找到所有大于0.1的峰值,这是因为驼峰区乱七八糟（当然这只针对该实验信号）
[top_v,loc_v] = findpeaks(-p);  % 默认。 'minpeakheight',0.1, 'minpeakdistance',50
top_v = -top_v;
scatter(loc_p,top_p);
scatter(loc_v,top_v);

%% 峰谷值2：当findpeaks无法精确定位实验信号的所有峰谷值时（峰谷值中包含翻转点，亦或是有一些异常点）,针对这些异常点的特殊处理
% arr = [501,1501,2501,3501];
% [loc_p,loc_v,top_p,top_v] = ditchReversePoing(arr,loc_p,loc_v,top_p,top_v);
% scatter(loc_p,top_p);
% scatter(loc_v,top_v);

%% 确定方向1：手动确定方向（最方便）yyds
% subplot(7,1,2);
% dir = solve_dir(arr,N);  % 输入翻转点arr数组(顺序输入) （方向有误乘-1即可）
% direction = dir;
% plot(p);
% hold on;
% plot(dir)

%% 确定方向2：包络确定方向！
subplot(7,1,2);
diffp = diff(p);  % diff是相邻两个数的差分，对第一个位置补0
diffp = [0,diffp];
% diff_acosp = diff(acos(p));
% diff_acosp = [0,diff_acosp];
% 非实验信号无需平滑化
% diffp = sgolayfilt(diffp,1,11);
% diffp = sgolayfilt(diffp,2,21);
% diffp = sgolayfilt(diffp,3,31);
plot(diffp);
hold on;
[top_diffp_p,loc_diffp_p] = findpeaks(diffp,'minpeakheight',0, 'minpeakdistance',150);  % 拿到极值和索引值
[top_diffp_v,loc_diffp_v] = findpeaks(-diffp,'minpeakheight',0, 'minpeakdistance',150);
top_diffp_v = -top_diffp_v;

scatter(loc_diffp_p,top_diffp_p);
scatter(loc_diffp_v,top_diffp_v);
en_top = interp1(loc_diffp_p,top_diffp_p,1:N,'spline');  % 三次样条插值，曲线更平滑
en_bottom = interp1(loc_diffp_v,top_diffp_v,1:N,'spline');

en_median = (en_top - (en_top - en_bottom)/2);
dir = -sign(en_median);  % 让dir暂时指明方向

plot(en_top);
plot(en_bottom);
plot(dir);
direction = dir;
title("diffp及其包络确定方向")


%% 重构：PUM transform
[phiF_reconstruct,Lt_reconstruct] = SMI_API_RECON_PUM(p,loc_p,loc_v,direction,N,lambda,c,alpha);
% [phiF_wrapped,phiF_reconstruct,Lt_reconstruct] = SMI_API_RECON_HT(p,direction,N,lambda,1,1);

%% 
% subplot(5, 1, 1);
% plot(phiF_wrapped);
% hold on;
% title("phiF_wrapped")

subplot(7,1,4);
plot(phiF_reconstruct)
title("重构后的phiF");

subplot(7,1,5);
plot(Lt);
hold on;
plot(Lt_reconstruct,'r')
title("重构后的信号");

%% 误差分析
subplot(7,1,6);
plot(Lt-Lt_reconstruct)
RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
title(['绝对误差，RMSE=', num2str(RMSE)])


















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