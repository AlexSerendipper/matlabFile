clc;
clear all;
close all;
%% 产生自混合信号
figure(1);
subplot(7,1,1);
fs = 200000;  % 采样率，即1/fs(s)采一个点。
N = 4000;  
fv = 100;  % 震动频率
C = [0.3,0.3];  % C设置一个从a到b变化的值
c=C(1);
alpha = 4;
[t, lambda, L0, Lt, phi0, p, c] = MOVE_API_HARMONIC(fs, N, fv, C, alpha);  % 1 简谐振动的自混合信号
p_init = p;
plot(t,p_init);
%% 得到重构所需的方向信息
[top_p, loc_p, top_v, loc_v, top_r, loc_r, direction] = SMI_API_FRINGE(p_init,N);  % 法1求方向，求的是初始信号的方向，所以都能用

%% （去直流）
subplot(7,1,2);
% dc1有效
[p2] = SMI_API_ELIMINATE_DC1(p_init,direction,"time");
% [p2] = SMI_API_ELIMINATE_DC2(p_init,direction,"time");
% [p1,p2] = SMI_API_evenlopEXTRACT_HT_PRO(p,N);

p = p2;
plot(p2);
title("去直流后的自混合信号");

%% 全局变量
windowLength = 512; % 窗长
V = 0.65; % 抑制因子

%% 时频分析
figure(2);
subplot(4,2,[1,2,3.4]);
[T,F,TF,TF_curb,p] = SMI_API_TFPM(p,N,fs,windowLength,V);

mesh(T,F,abs(TF)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' ); 
view(0,90); % 设置初始视角为俯视角
title('抑制前');

subplot(4,2,[5,6,7,8]);
mesh(T,F,abs(TF_curb)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' );
view(0,90); % 设置初始视角为俯视角
title('抑制后');

%% 某时刻频谱分析
% plot(1:length(F), abs(TF(:,20)));
% plot(1:length(F), abs(TF_curb(:,20)))

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
scatter(loc_p,top_p);
scatter(loc_v,top_v);
title('时频谱处理后的时域信号');

%% 希尔伯特变换重构
subplot(7,1,4);
[phiF_wrapped,phiF_reconstruct,Lt_reconstruct] = SMI_API_RECON_HT(p,direction,N,lambda,0,0);
plot(phiF_wrapped);
title("phiF_wrapped")

subplot(7,1,5);
plot(phiF_reconstruct);
title("重构后的phiF");

subplot(7,1,6);
plot(Lt);
hold on;
plot(Lt_reconstruct,'r')
title("重构后的信号");

%% 误差分析
subplot(7,1,7);
plot((Lt-Lt_reconstruct)); 
RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
title(['绝对误差，RMSE=', num2str(RMSE)]);

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
