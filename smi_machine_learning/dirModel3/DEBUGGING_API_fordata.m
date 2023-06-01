%% 获得降C后的自混合信号
%% 若要根据降C后的数据作训练集，思路是根据方向得到方向变换的点，501、1501...
%  若方向变换的点附近(70)有谷值，将那点替换成谷值，此时方向变换的点当作翻转点来用就好了

clc;
clear all;
close all;
%% 产生自混合信号
C=1.8;
fs = 200000;  % 采样率，即1/fs(s)采一个点。
N = 4000;  
fv = 100;  % 震动频率
alpha = 5;
[t, lambda, L0, Lt, phi0, p, c] = MOVE_API_HARMONIC(fs, N, fv, C, alpha);  % 1 简谐振动的自混合信号
p_init = p;

%% 得到重构所需的方向信息
[top_ov,loc_ov,top_p, loc_p, top_v, loc_v, top_r, loc_r, direction] = SMI_API_FRINGE(p,N);
direction_seg1 = [];  % 根据方向信息，求出方向变化的点
for i = 1:length(direction)-1
    if direction(i) ~=  direction(i+1)
        direction_seg1 = [direction_seg1, i];
    end
end

%% （去直流）
% dc1有效
[p] = SMI_API_ELIMINATE_DC1(p_init,direction,"time");
%     [p] = SMI_API_ELIMINATE_DC2(p_init,direction,"time");
%     [p1,p2] = SMI_API_evenlopEXTRACT_HT_PRO(p,N);
%     p = p2;

%% 全局变量
windowLength = 512; % 窗长
V = 0.65; % 抑制因子

%% 时频分析
[T,F,TF,TF_curb,p] = SMI_API_TFPM(p,N,fs,windowLength,V);

%% 归一化(p = 2 * (p - -abs(hilbert(-p)))./(abs(hilbert(p)) - -abs(hilbert(-p))) - 1;)
p = sgolayfilt(p,1,11);
p = sgolayfilt(p,2,21);
p = sgolayfilt(p,3,31);


[top_p,loc_p] = findpeaks(p);  % 找到所有大于0.1的峰值,这是因为驼峰区乱七八糟（当然这只针对该实验信号）
[top_v,loc_v] = findpeaks(-p);  % 默认。 'minpeakheight',0.1, 'minpeakdistance',50
top_v = -top_v;


subplot(2,1,1)
plot(p);
hold on;
plot(direction);
hold on;
scatter(loc_p,top_p);
scatter(loc_v,top_v);

%% 方向变换的点附近有谷值
for i = 1:length(direction_seg1)
    for j = 1:length(loc_v)
        if direction_seg1(i)>loc_v(j)-70  && direction_seg1(i)<loc_v(j)+70  % 附近有谷值
             direction_seg1(i) = loc_v(j);
        end
    end
end
% 此时direction_seg1，即作为翻转点即可~~~
loc_r = direction_seg1;

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
    int_ = 1000;
    loc_ = [SMI_API_RESAMPLE(loc,int_) dir];  % 方向作为最后一个点
    fringeData = [fringeData;loc_];
end


subplot(5, 1, 3);
plot(fringeData(6,:))