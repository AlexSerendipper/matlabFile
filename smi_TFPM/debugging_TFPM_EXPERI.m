%% 对于实验信号，需要先进行时频谱抑制才能判断方向~ 所以步骤与仿真实验不同
clc;
clear all;
close all;
set(0,'defaultfigurecolor','w');  % 设置画布背景色为白色

%% 全局变量
windowLength = 128; % 窗长
overlapLength = floor(windowLength * 0.9);  % OverlapLength后为指定的重叠长度
window = hamming(windowLength, "periodic");  % 使用汉明窗作为滑动的窗口
fftLength = 5*windowLength;
V = 0.65; % 抑制因子2

%% 产生自混合信号
fs = 200000;  % 采样率，即1/fs(s)采一个点。
N = 4000;  
fv = 200;  % 震动频率
C = [1.5];  % C设置一个从a到b变化的值
alpha = 4;

% [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_HARMONIC(fs, N, fv, C, alpha);  % 1 简谐振动的自混合信号

% [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_SIN_MODULATE(fs, N, C, alpha);  % 2 正弦调制信号的自混合信号

% cut = 200; [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_ALEATORY(fs, N, cut, C, alpha);  % 3 产生随机振动时，方向×负，输入一个能被N整除的数，将N分为N/cut段

% [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_ALEATORY_LOAD(fs, N, C, alpha);  % 4 加载储存好的随机振动，方向×负

% path =  'D:\matlab save\self-mixing\smi_脉搏波\scope_1.csv';  % 5 文件路径, M/N/win/w = 250515/16000/128/150
path =  'D:\matlab save\self-mixing\smi_TFPM\temp\data\音叉\TEK00008.csv'; 
M = 10000; N = 5000; lambda = 650e-9; [t, p, fs] = MOVE_API_EXPERIMENT2(M, N, path);  % 5 从M点处取N个点

p_init = p;
% p = awgn(p,10); % 10db，加噪声在抑制是可以的，问题是加了噪声没法从原信号判断方向了
figure(1);
set(gcf,'position',[15, 30, 800, 800]);
subplot(7,1,1);
plot(p); ylabel('p(t)(nm)','FontSize',7,'FontWeight','bold'); xlabel('sample(num)','FontSize',7,'FontWeight','bold');  % ✔这玩意赶紧学 
title(['自混合信号,C=',num2str(C), 'alpha=',num2str(alpha)]);
% p = p - mean(p);


%% 时频分析
figure(2);
subplot(4,2,[1,2,3.4]);
% load("p.mat");
% N = length(p);
% fs = 50000;
% windowLength = 128;
% p(1:N) = p(1:N)+0.2;
% % f1 = @(x) -2/N^2 * x.^2 + 2/N * x;
% f2 = @(x) 1/600000 * x;
% p = p + f2(1:N);
% p = imag(hilbert(p));
[T,F,TF,TF_curb,p] = SMI_API_TFPM(p,N,fs,windowLength,V);
mesh(abs(TF)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' ); 
colormap jet;
% colormap Cold; 
view(0,90); % 设置初始视角为俯视角
title('抑制前');

subplot(4,2,[5,6,7,8]);
mesh(abs(TF_curb)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' );
% colormap Hot; 
view(0,90); % 设置初始视角为俯视角
title('抑制后');


%% 归一化(p = 2 * (p - -abs(hilbert(-p)))./(abs(hilbert(p)) - -abs(hilbert(-p))) - 1;)
figure(1);
subplot(7,1,3);
plot(p);

p = sgolayfilt(p,1,11);
p = sgolayfilt(p,2,21);
p = sgolayfilt(p,3,31);
hold on;
[top_p,loc_p] = findpeaks(p);  % 找到所有大于0.1的峰值,这是因为驼峰区乱七八糟（当然这只针对该实验信号）
[top_v,loc_v] = findpeaks(-p);  % 默认。 'minpeakheight',0.1, 'minpeakdistance',50
top_v = -top_v;
% scatter(loc_p,top_p);
% scatter(loc_v,top_v);
title('时频谱处理后的时域信号');

%% 得到重构所需的方向信息
arr = [214,842,1471,2101,2719,3356,3977,4591];
direction = SMI_API_DIR(arr,N);
% [top_p, loc_p, top_v, loc_v, top_r, loc_r, direction] = SMI_API_FRINGE(p_init,N);  % 法1求方向，求的是初始信号的方向，所以都能用
% W = 150; direction = SMI_API_TFPM_DIRECTION(p,N,W);  % 法2求方向，利用平均峰值距离（要求均为w形的驼峰），且必须用归一化2的信号
% direction = -direction;  % 如果初始震动用的cos，那方向是反的，一定要注意
plot(direction);

%% 计算一下抑制后，C的值
% [C_reconstruct, alpha_reconstruct] = SMI_API_ESTIMATE_C(phiF_reconstruct); 


%% 希尔伯特变换重构
subplot(7,1,4);
[phiF_wrapped,phiF_reconstruct,Lt_reconstruct] = SMI_API_RECON_HT(p,direction,N,lambda,0,0);
plot(phiF_wrapped);
title("phiF_wrapped")

subplot(7,1,5);
plot(phiF_reconstruct);
title("重构后的phiF");

subplot(7,1,6);
% plot(Lt);
% hold on;
% plot(Lt_reconstruct,'r')
% title("重构后的信号");

%% 误差分析
subplot(7,1,7);
% plot((Lt-Lt_reconstruct)); 
% RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
% title(['绝对误差，RMSE=', num2str(RMSE)]);

%% 这部分就是我个人的特殊用途了
% a1 = phiF_reconstruct + 0*sin(phiF_reconstruct+atan(alpha));
% b1 = a1 * lambda / (4 * pi);
% b1 = b1 - mean(b1);
% c1 =  sqrt(mean((Lt-b1).^2));
% 
% a2 = phiF_reconstruct + 0.01*sin(phiF_reconstruct+atan(alpha));
% b2 = a2 * lambda / (4 * pi);
% b2 = b2 - mean(b2);
% c2 =  sqrt(mean((Lt-b2).^2));
% 
% a3 = phiF_reconstruct + 0.02*sin(phiF_reconstruct+atan(alpha));
% b3 = a3 * lambda / (4 * pi);
% b3 = b3 - mean(b3);
% c3 =  sqrt(mean((Lt-b3).^2));
