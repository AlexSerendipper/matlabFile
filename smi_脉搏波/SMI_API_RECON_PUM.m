clc;
clear all;
close all;
%% 产生自混合信号
subplot(5, 1, 1);
fs = 200000;  % 采样率，即fs(s)采一个点。
N = 4000;
fv = 100;  % 震动频
C = [0.5];  % C设置一个从a到b变化的值
alpha = 5;
% [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_HARMONIC(fs, N, fv, C, alpha);  % 1 简谐振动的自混合信号

% [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_SIN_MODULATE(fs, N, C, alpha);  % 2 正弦调制信号的自混合信号

% cut = 200; [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_ALEATORY(fs, N, cut, C, alpha);  % 3 产生随机振动时，方向×负，输入一个能被N整除的数，将N分为N/cut段

% [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_ALEATORY_LOAD(fs, N, C, alpha);  % 4 加载储存好的随机振动，方向×负

% path =  'D:\matlab save\self-mixing\smi_实验信号\gexiong_f(100)_A(3um).csv';  % 5 文件路径
% M = 20000; N = 16000;  [t, p, fs] = MOVE_API_EXPERIMENT(M, N, path);  % 5 从M点处取N个点

[t, lambda, L0, Lt, phi0, p, c] = MOVE_API_PPG_LOAD(fs, N, C, alpha);  % 加载脉搏波

plot(Lt);
title(['外部简谐振动,C=',num2str(C)]);
subplot(5, 1, 2);
% p = awgn(p,40);  % 10db，加高斯白噪声
% p = p .* (1+0.2*cos(2*pi*75*t));  % 给自混合信号加包络，加了一个幅值为0.2，频率为75的包络
plot(p(1:1000));
title("自混合信号");

%% 得到重构所需的相关信息
[top_p, loc_p, top_v, loc_v, top_r, loc_r, direction] = SMI_API_FRINGE(p,N);
direction = -direction;  % 如果初始震动用的cos，或采样随机振动，方向×负，一定要注意!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
plot(p);
hold on;
scatter(loc_p,top_p);
scatter(loc_v,top_v);
scatter(loc_r,top_r);
plot(direction);
title("pks, vls, rev")


%% 重构
subplot(5, 1, 3);
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
figure(1);
subplot(5,1,4);
phi0_reconstrut = phiF_reconstruct + C_reconstruct*sin(phiF_reconstruct+atan(alpha_reconstruct));  % 这里的alpha如果用估算的，就会引入蛮大的误差
Lt_reconstruct = phi0_reconstrut * lambda / (4 * pi);


Lt_reconstruct = Lt_reconstruct - mean(Lt_reconstruct); % 简写振动和余弦调制振动，重构后加上幅值A
% Lt_reconstruct = Lt_reconstruct + 1.5 * lambda;  % 重构后的随机振动信号要加上幅值1.5的波长，这是为啥我页不知道
plot(Lt,'k')
hold on;
plot(Lt_reconstruct,'r')
title("重构后的信号");

%% 误差分析
subplot(5,1,5);
plot(Lt-Lt_reconstruct)
RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
title(['绝对误差，RMSE=', num2str(RMSE)])


%% 可视化区域
figure(2);
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

