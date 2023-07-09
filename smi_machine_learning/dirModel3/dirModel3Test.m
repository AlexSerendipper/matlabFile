
%% 产生自混合信号
clc;
clear all;
close all;
C=2.5;
fs = 200000;  % 采样率，即1/fs(s)采一个点。
N = 4000;  
fv =  100;  % 震动频率
alpha = 5;
[t, lambda, L0, Lt, phi0, p, c] = MOVE_API_HARMONIC(fs, N, fv, C, alpha);  % 1 简谐振动的自混合信号
p_init = p;
    
% 得到重构所需的方向信息
[top_ov,loc_ov,top_p, loc_p, top_v, loc_v, top_r, loc_r, direction] = SMI_API_FRINGE(p,N);
direction_seg1 = [];  % 根据方向信息，求出方向变化的点
for i = 1:length(direction)-1
    if direction(i) ~=  direction(i+1)
        direction_seg1 = [direction_seg1, i];
    end
end

% （去直流）
% dc1有效
[p] = SMI_API_ELIMINATE_DC1(p_init,direction,"time");
%     [p] = SMI_API_ELIMINATE_DC2(p_init,direction,"time");
%     [p1,p2] = SMI_API_evenlopEXTRACT_HT_PRO(p,N);
%     p = p2;

% 全局变量
windowLength = 512; % 窗长
V = 0.65; % 抑制因子

% 时频分析
[T,F,TF,TF_curb,p] = SMI_API_TFPM(p,N,fs,windowLength,V);

% 归一化(p = 2 * (p - -abs(hilbert(-p)))./(abs(hilbert(p)) - -abs(hilbert(-p))) - 1;)
p = sgolayfilt(p,1,11);
p = sgolayfilt(p,2,21);
p = sgolayfilt(p,3,31);


[top_p,loc_p] = findpeaks(p);  % 找到所有大于0.1的峰值,这是因为驼峰区乱七八糟（当然这只针对该实验信号）
[top_v,loc_v] = findpeaks(-p);  % 默认。 'minpeakheight',0.1, 'minpeakdistance',50
top_v = -top_v;
loc_ov = loc_v;

subplot(5,1,1)
plot(p);
hold on;
plot(direction);
hold on;
scatter(loc_p,top_p);
scatter(loc_v,top_v);
title(['自混合信号，C=', num2str(C)]);
%% 基于神经网络的方向判断（注意要保证loc_v的连续性）
load("net.mat")
dir = zeros(1,N);
loc_r = [];
rect = [];
int_ = 30;  % 插值倍数(重采样放大倍数，默认为1000)
fs = int_;  % 因为原本采样率为1，放大后采样率就是int_
for i = 2:length(loc_ov)
    N = loc_ov(i)-loc_ov(i-1);
    % 将原信号重采样为原来的int_倍
    p_ = interp1(1:int_:int_* N,p(loc_ov(i-1):loc_ov(i)-1),1:int_* N,'spline');
    % 再将信号降采样为原来的N倍，即 N * 1000 / N， 故最后信号长度为int_
    p_ = p_(1:(loc_ov(i)-loc_ov(i-1)):length(p_));
    
    judge = sim(net,p_');  % 基于神经网络判断的方向
    [val,idx] = max(judge);
    if(idx==1)
        judge = -1;
    elseif(idx==2)   
        judge = 0;
    elseif(idx==3)
        judge = 1;
    end
    

    if(i==2)
        dir(1:loc_ov(i)) =  judge;
        lastJudge = judge;
    elseif(i==length(loc_ov))   
        dir(loc_ov(i-1):end) =  judge;
    else
        % 先获取下一段的方向
        N2 = loc_ov(i+1)-loc_ov(i);
        % 将原信号重采样为原来的int_倍
        p2_ = interp1(1:int_:int_* N2,p(loc_ov(i):loc_ov(i+1)-1),1:int_* N2,'spline');
        % 再将信号降采样为原来的N倍，即 N * 1000 / N， 故最后信号长度为int_
        p2_ = p2_(1:(loc_ov(i+1)-loc_ov(i)):length(p2_));
        
        nextJudge = sim(net,p2_');  % 基于神经网络判断的方向
        [val,idx] = max(nextJudge);
        if(idx==1)
            nextJudge = -1;
        elseif(idx==2)   
            nextJudge = 0;
        elseif(idx==3)
            nextJudge = 1;
        end

        %% 处理个别判断错误（只能处理 单个判断1为-1的错误，若连续出现则不行）
        if(judge == -lastJudge && judge == -nextJudge)  % 偶尔判断失误的处理, 即与上一段与下一段均均为要当前段方向相反，则使用相反的方向即可
           dir(loc_ov(i-1):loc_ov(i)) = -judge;  
        else
           dir(loc_ov(i-1):loc_ov(i)) =  judge;
           lastJudge = judge;
        end 
%          dir(loc_ov(i-1):loc_ov(i)) =  judge;
    end
end
subplot(5, 1, 2);
plot(p);
hold on;
plot(dir);
title("初始判断的方向");


%% 处理判断为驼峰区 
i = 2;
N = length(dir);
while i<N
    if(dir(i)==0)
        pre = dir(i-1);  % 记录前值
        % j = i + 1;
        for j = i+1:N
            if(dir(j)~=0)
                break
            end
        end
        n = floor((j-i)/2);
        
        
        % 将0的前半部分设为前值
        for k = i:i+n-1
            dir(k) = pre;
        end

        % 将0的后半部分设为前值取反
        for k = i+n : j-1
            dir(k) = -pre;
        end
        i = i + n -1;
    end
    i = i+1;
end
subplot(5, 1, 3);
plot(p);
hold on;
plot(dir);

















