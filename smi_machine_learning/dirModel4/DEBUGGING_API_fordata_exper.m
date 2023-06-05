clc;
clear all;
close all;
%% 产生自混合信号
% fs = 200000;  % 采样率，即fs(s)采一个点。
% N = 8000;  
% fv = 50;  % 震动频
% C = [0.5];  % C设置一个从a到b变化的值
% alpha = 5;
% [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_HARMONIC(fs, N, fv, C, alpha);  % 1 简谐振动的自混合信号
% subplot(6, 1, 1);
% % p = p .* (1+0.2*cos(2*pi*75*t));  % 给自混合信号加包络，加了一个幅值为0.2，频率为75的包络
% p = -1 + (p-min(p))/(max(p)-min(p))*2;
% p = awgn(p,20); 
% plot(p);
% hold on;
% [top_v,loc_v] = findpeaks(-p);  % 'minpeakdistance'
% top_v = -top_v;
% scatter(loc_v,top_v);
% title("自混合信号");


%% 实验信号
subplot(6, 1, 1);
path =  'D:\matlab save\smi_实验信号\f(100)_A(3um).csv';  % 5 文件路径, M/N/win/w = 250515/16000/128/150
M = 40000; N = 10000; lambda = 650e-9; [t, p, fs] = MOVE_API_EXPERIMENT(M, N, path);  % 5 从M点处取N个点
plot(p);

%% first step ,calculate the MPD
subplot(6, 1, 2);
[c, lags] = xcorr(p,'unbiased');
c = sgolayfilt(c,1,11);
c = sgolayfilt(c,2,21);
c = sgolayfilt(c,3,31);

[top_c,loc_c] = findpeaks(c);
MPD = mean(diff(loc_c));
plot(c);
hold on;
scatter(loc_c,top_c);

subplot(6, 1, 3);
plot(p);
hold on;
[top_v,loc_v] = findpeaks(-p, 'minpeakdistance',MPD/2);
top_v = -top_v;
scatter(loc_v,top_v);



%% second step,zero-crossing,所有大于0的都不要了
subplot(6, 1, 4);
loc_ov=[];
top_ov=[];
for i = 1:length(loc_v)
    if(top_v(i)<0)
        loc_ov = [loc_ov loc_v(i)];
        top_ov = [top_ov top_v(i)];
    end
end

plot(p);
hold on;
scatter(loc_ov,top_ov);

%% third step,average && variance
avg = mean(top_ov);  % 计算所有谷值点的平均值
mea = mean(abs(top_ov - avg));  % 计算所有点 距离 谷点平均值 的平均值
for i = 1:length(loc_ov)
    temp = top_ov(i) - avg;
    if (temp > mea)
        loc_ov(i) = nan;
        top_ov(i) = nan;
    end
end
loc_ov(isnan(loc_ov))=[];  
top_ov(isnan(top_ov))=[]; 
subplot(6, 1, 5);
plot(p);
hold on;
scatter(loc_ov,top_ov);

%% MDP judge,根据平均峰值判定，小于平均峰值的取更小点保留。。。
% 需要多加一个判断是否两点的差值>mea/2，因为当C比较小的时候，是不需要平均峰值这一步的
for i = 1:length(loc_ov)-1
    if(loc_ov(i+1)-loc_ov(i)<MPD/1.25)
        if(top_ov(i+1)>top_ov(i))
            top_ov(i+1) = nan;
            loc_ov(i+1) = nan;
            i = i+2;
        else
            top_ov(i) = nan;
            loc_ov(i) = nan;
            i = i+2;
        end
    end
end
subplot(6, 1, 6);
plot(p);
hold on;
scatter(loc_ov,top_ov);
