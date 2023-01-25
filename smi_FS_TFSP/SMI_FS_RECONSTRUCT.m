%% FS/TFSP就是对phiF做傅里叶变换，抑制后(得到主频)，逆傅里叶变换出来的值近似为phi0咯。 
% 而TFPM就是对p做傅里叶变换，抑制后逆傅里叶变换出来的值近似于C=0！！！
clc;
clear all;
close all;

%% 产生自混合信号
subplot(5, 1, 1);
fs = 128000;  % 采样率，即fs(s)采一个点
N = 4000;  
fv = 128;  % 震动频
C = 2;
alpha = 4.6;
[t, lambda, L0, Lt, phi0, p] = SMI_API_HARMONIC(fs, N, fv, C, alpha);  % 简谐振动的自混合信号
% [t, lambda, L0, Lt, phi0, p] = SMI_API_SIN(fs, N, C, alpha);  % 正弦调制信号的自混合信号

% cut = 500;  % cut降采样，输入一个能被N整除的数，将N分为N/cut段
% [t, lambda, L0, Lt, phi0, p] = SMI_API_ALEATORY(fs, N, cut, C, alpha);  % 使用随机振动时，方向×负

plot(Lt);
title(['外部简谐振动,C=',num2str(C)]);
subplot(5, 1, 2);
plot(p);
title("自混合信号");

%% 得到重构所需的相关信息
[top_p, loc_p, top_v, loc_v, top_r, loc_r, direction] = SMI_FRINGE(p,N);
% direction = -direction;  % 如果初始震动用的cos，或采样随机振动，方向×负，一定要注意!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
plot(p);
hold on;
scatter(loc_p,top_p);
scatter(loc_v,top_v);
scatter(loc_r,top_r);
plot(direction);
title("pks, vls, rev")

%% 重构
% subplot(5, 1, 3);
acosp = acos(p);
% plot(acosp);
% title("acosp");
% hold on;

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

% plot(acosp_op1);
hold on;
% title("翻转点×-1");


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
% subplot(5,1,4);
% plot(phiF_reconstruct);
% title(phiF)

%% 幅度谱
subplot(5, 1, 3);
phiF_reconstruct = phiF_reconstruct - mean(phiF_reconstruct);  % 消除直流分量
% w = hamming(N);
% p_ = fft(w'.* p, N) * 2;  % 加窗傅里叶变换，具体为何这样写 我也不懂
f = fs / N * (0 : 1 : N-1);  % Fs/N就是这个频谱中的最小频率间隔！！！！！所以N越大，分辨率会越高
phiF_reconstruct_ = fft(phiF_reconstruct, N);
Fv = phiF_reconstruct_;  % 用来储存傅里叶变换后的值
Av = abs(phiF_reconstruct_) * 2 / N;
plot(f(1:N/2), Av(1:N/2));  % fft算法默认是双边谱,通常我们只取一半
title('自混合信号的幅度谱');


%% （FS是frequent sampling）频率采样，利用该原理可以绕过C和alpha的估计实现重构（求出Xa近似于phi0）
maxAv = max(Av);
NAv = (Av/maxAv) - 10^-10;  % 10^-10，保证后边不会求出无穷小
Fl = NAv.^1 ./ (1 - NAv.^1);  % 权重函数，频谱中原来大的更大，小的更小
std_Fl = std(Fl);
Fl_filtered = [];
Fv_filtered = [];

% 选出Fl中大于其标准差的成分，因为是双边谱，所以有两个一样的值
for i = 1:N
    if Fl(i) > std_Fl
        Fl_filtered = [Fl_filtered, Fl(i)];
        Fv_filtered = [Fv_filtered, Fv(i)];
    else
        Fv_filtered = [Fv_filtered, 0];
    end  
end

plot(f(1:N/2), -Fv_filtered(1:N/2))
phia = ifft(Fv_filtered);

subplot(5,1,4);
phi0_reconstrut = phia;
Lt_reconstruct = phi0_reconstrut * lambda / (4 * pi);
plot(Lt_reconstruct);
title("重构后的信号");
hold on;
plot(Lt);

subplot(5,1,5);
plot(Lt-Lt_reconstruct);
RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
title(['绝对误差，RMSE=', num2str(RMSE)])









