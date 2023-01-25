clc;
clear all;
close all;
%% 产生自混合信号
fs = 200000; % 采样率
N = 4000;  % 采样点
fv = 100;  % 震动频率
C = [1.6, 1.6];  
alpha = 4;
[t, lambda, L0, Lt, phi0, p] = MOVE_API_HARMONIC(fs, N, fv, C, alpha);  % 1 简谐振动的自混合信号
% [t, lambda, L0, Lt, phi0, p] = MOVE_API_COS(fs, N, C, alpha);  % 2 余弦调制信号的自混合信号

% cut = 200;  % cut降采样，输入一个能被N整除的数，将N分为N/cut段  % 3
% [t, lambda, L0, Lt, phi0, p] = MOVE_API_ALEATORY(fs, N, cut, C, alpha);  % 3 产生随机振动时，方向×负

% [t, lambda, L0, Lt, phi0, p] = MOVE_API_ALEATORY_LOAD(fs, N, C, alpha);  % 4 加载储存好的随机振动，方向×负

%% 得到重构所需的相关信息
[top_p, loc_p, top_v, loc_v, top_r, loc_r, direction] = SMI_API_FRINGE(p,N);
% direction = direction;  % 如果初始震动用的cos，那方向是反的，一定要注意

%% early
figure(1);
subplot(6, 1, 1);
plot(Lt);
title(['外部简谐振动,C=',num2str(C)]);
subplot(6, 1, 2);
title("自混合信号");
% p = awgn(p,10);  % 10db，加高斯白噪声
% p = p .* (1+0.2*cos(2*pi*75*t));  % 给自混合信号加包络，加了一个幅值为0.2，频率为75的包络
plot(p);
hold on;
scatter(loc_p,top_p);
scatter(loc_v,top_v);
scatter(loc_r,top_r);
plot(direction);
title("pks, vls, rev")

%% hilbert transform
subplot(6, 1, 3);
inverse_hb = (imag(hilbert(p))) .* direction;  % 得到sin 
plot(inverse_hb,"b");
hold on;
plot(p,'r');
title("blue,inverse-HT and orginal-HT","red,self-mixing signal");

subplot(6, 1, 4);
phiF_wrapped = atan(inverse_hb./p);
plot(phiF_wrapped);
hold on;
title("phiF_wrapped")

%% arctan相位解包裹
for i = 2:N
    if (phiF_wrapped(i) - phiF_wrapped(i-1) > pi/2)
        phiF_wrapped(i:end) = phiF_wrapped(i:end) -  pi;
    elseif (phiF_wrapped(i) - phiF_wrapped(i-1) < -pi/2)  % 找峰值点，碰到峰值则加pi
        phiF_wrapped(i:end) = phiF_wrapped(i:end) +  pi;            
    end
end

phiF_reconstruct = phiF_wrapped;
% subplot(6, 1, 5);
% plot(phiF_reconstruct)
% title("重构后的phiF")

%% 参数估算
[C_reconstruct, alpha_reconstruct] = SMI_API_ESTIMATE_C(phiF_reconstruct);  

%% 重构
subplot(6, 1, 5);
phi0_reconstrut = phiF_reconstruct + C_reconstruct*sin(phiF_reconstruct+atan(alpha));
Lt_reconstruct = phi0_reconstrut * lambda / (4 * pi);
plot(Lt,'k')
hold on;

Lt_reconstruct = Lt_reconstruct - mean(Lt_reconstruct); % 简写振动和余弦调制振动，重构后加上幅值A
% Lt_reconstruct = Lt_reconstruct + 1.5 * lambda;  % 重构后的随机振动信号要加上幅值1.5的波长，这是为啥我页不知道

plot(Lt_reconstruct,'r')
title("重构后的信号");

%% 误差分析
subplot(5,1,5);
plot(Lt-Lt_reconstruct)
RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
title(['绝对误差，RMSE=', num2str(RMSE)])
