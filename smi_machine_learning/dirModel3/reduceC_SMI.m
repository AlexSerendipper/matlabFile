%% 获得降C后的自混合信号
%% 若要根据降C后的数据作训练集，思路是根据方向得到方向变换的点，501、1501...
%  若方向变换的点附近(70)有谷值，将那点替换成谷值，此时方向变换的点当作翻转点来用就好了
function [p] = reduceC_SMI()
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
    
%     subplot(2,1,1)
%     plot(p);
%     hold on;
%     plot(direction);
%     hold on;
%     scatter(loc_p,top_p);
%     scatter(loc_v,top_v);
end

