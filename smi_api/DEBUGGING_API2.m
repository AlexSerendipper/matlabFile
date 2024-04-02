%% 该程序用于调试，其实也是非常好用且重要的

%% 产生自混合信号
% addpath(genpath('..\smi_api'));  % 添加文件夹到搜索目录
% rmpath(genpath('..\smi_api'));  % 删除文件夹到搜索目录
clc;
clear all;
close all;
subplot(7,1,1);
fs = 200000; % 采样率
N = 4000;  % 采样点
fv = 100;  % 震动频率
C = [0.7];
c=C(1);
alpha = 5;
[t, lambda, L0, Lt, phi0, phiF, p, c] = MOVE_API_HARMONIC(fs, N, fv, C, alpha);  % 1 简谐振动的自混合信号
% [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_SIN_MODULATE_FM(fs, N, C, alpha);   % 2 余弦调制信号的自混合信号
cut = 500;  % cut降采样，输入一个能被N整除的数，将N分为N/cut段  % 3
% [t, lambda, L0, Lt, phi0, p] = MOVE_API_ALEATORY(fs, N, cut, C, alpha);  % 3 产生随机振动时，方向×负
% [t, lambda, L0, Lt, phi0, p] = MOVE_API_ALEATORY_LOAD(fs, N, C, alpha);  % 4 加载储存好的随机振动，方向×负

% [p,PS] = mapminmax(p,-1,1);  % 归一化函数
plot(Lt);
subplot(7,1,2);
plot(p);
% p = awgn(p,40);  % 10db，加高斯白噪声
% p = p .* (1+0.2*cos(2*pi*75*t));  % 给自混合信号加包络，加了一个幅值为0.2，频率为75的包络
hold on;
title(['自混合信号，C=', num2str(c)]);

%% 条纹检测
[top_p, loc_p, top_v, loc_v, top_r, loc_r, direction,direc] = SMI_API_FRINGE(p,N);
subplot(7,1,3);
plot(p);
hold on;
scatter(loc_p,top_p);
scatter(loc_v,top_v);
scatter(loc_r,top_r);
plot(direction);
title("得到纯净的峰谷值和完全正确的方向")
 
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