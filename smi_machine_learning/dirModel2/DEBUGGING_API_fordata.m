%% 制作普通自混合信号的数据集
%% 重采样为1000个点太多了！！大概20~50之间就够！
clc;
clear all;
close all;
%% 产生自混合信号
subplot(5, 1, 1);
fs = 200000;  % 采样率，即fs(s)采一个点。
N = 8000;  
fv = 50;  % 震动频
C = [0.2];  % C设置一个从a到b变化的值
alpha = 5;
[t, lambda, L0, Lt, phi0, p, c] = MOVE_API_HARMONIC(fs, N, fv, C, alpha);  % 1 简谐振动的自混合信号
plot(Lt);
title(['外部简谐振动,C=',num2str(C)]);
subplot(5, 1, 2);
% p = awgn(p,40);  % 10db，加高斯白噪声
% p = p .* (1+0.2*cos(2*pi*75*t));  % 给自混合信号加包络，加了一个幅值为0.2，频率为75的包络
plot(p);
title("自混合信号");

%% 得到重构所需的相关信息
[top_ov,loc_ov,top_p, loc_p, top_v, loc_v, top_r, loc_r, direction] = SMI_API_FRINGE(p,N);
% direction = -direction;  % 如果初始震动用的cos，或采样随机振动，方向×负，一定要注意!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
plot(p);
hold on;
scatter(loc_ov,top_ov);
% scatter(loc_p,top_p);
% scatter(loc_v,top_v);
% scatter(loc_r,top_r);
plot(direction);
title("pks, vls, rev")

%% 拿信号
fringeData = [];
for i=2:length(loc_ov)
    % 得到谷值区间点数（条纹）
    N = loc_ov(i) - loc_ov(i-1);
    % 拿到谷值区间内所有的点
    dd = loc_ov(i-1)+1:loc_ov(i)-1;
    % 使用众数记录当前区间方向
    m = mode(direction(loc_ov(i-1):loc_ov(i)));
    % 判断谷值区间 内！，是否含有翻转点，如果有则 if=false
    if isempty(intersect(dd, loc_r))
        % 不包含，记录当前方向
        dir = m;
    else
        % 包含,则方向置零
        dir = 0;
    end
    loc = p(loc_ov(i-1):loc_ov(i));
    int_ = 30;
    loc_ = [SMI_API_RESAMPLE(loc,int_) dir];
    fringeData = [fringeData;loc_];
end

subplot(5, 1, 3);
plot(fringeData(3,:))
subplot(5, 1, 4);
plot(fringeData(4,:))
subplot(5, 1, 5);
plot(fringeData(5,:))