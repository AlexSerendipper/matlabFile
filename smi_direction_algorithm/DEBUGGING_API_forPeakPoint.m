%% 纯通过算法的方式来判断降C后的自混合信号的方向：思路为，通过算法计算出带噪自混合信号的所有峰值点，
%   与降C后的自混合信号的峰值点进行比较即可~~驼峰区位置可能会出现方向判断错误的问题，但是由于驼峰区的方向重构没有影响，所以可以忽略错误
%   这个方法的劣势也十分明显，就是C值如果原本就是特别小的情况下，很难通过对比来得出方向，此外，如果耐噪范围也变小了，只能在25db以上的噪声下工作
%   并且这个方法对含有 过大的散斑 的信号鲁棒性不强
%   由于此方法要求保证驼峰区的正确性，所以建议随机振动
clc;
clear all;
close all;
%% 产生自混合信号
fs = 200000;  % 采样率，即fs(s)采一个点。
N = 4000;  
fv = 100;  % 震动频
C = [0.7,1.7];  % C设置一个从a到b变化的值
uptraParameter1 = 1.4; % MPD<MPD/uptralparameter1，用于设置降C前信号的平均峰值距离，相邻两点若小于该距离，取较大的点保留！
uptraParameter2 = 2; % MPD>MPD2*uptralparameter2，用于设置识别降C后信号的驼峰区的宽度，小于该宽度视为驼峰区

alpha = 5;
[t, lambda, L0, Lt, phi0, phiF, p, c] = MOVE_API_HARMONIC(fs, N, fv, C, alpha);  % 1 简谐振动的自混合信号
% [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_SIN_MODULATE(fs, N, C, alpha);
% [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_ALEATORY_LOAD(fs, N, C, alpha);

p = p .* (1+0.3*cos(2*pi*75*t));  % 给自混合信号加包络，加了一个幅值为0.2，频率为75的包络
% p = -1 + (p-min(p))/(max(p)-min(p))*2;
p = awgn(p,30);   % 噪声低于25，求出的峰值不够准确，并且在驼峰区位置会出现误差点

[top_p, loc_p, top_v, loc_v, top_r, loc_r, direction, direc, diffp, en_bottom, en_median, en_top] = SMI_API_FRINGE(p,N);
subplot(7, 1, 1);
plot(Lt);
subplot(7, 1, 2);
plot(p);
hold on;
[top_v,loc_v] = findpeaks(-p);  % 'minpeakdistance'
top_v = -top_v;
scatter(loc_v,top_v);loc_op0=loc_v;top_op0=top_v;
title("自混合信号");
p_init = p;

%% 全局变量
windowLength = 200; % 窗长
V = 0.65; % 抑制因子

%% first step ,calculate the MPD
subplot(7, 1, 3);
[c, lags] = xcorr(p,'unbiased');
c = sgolayfilt(c,1,11);
c = sgolayfilt(c,2,21);
c = sgolayfilt(c,3,31);

[top_c,loc_c] = findpeaks(c);
MPD = mean(diff(loc_c));
plot(c);
hold on;
scatter(loc_c,top_c);
title("根据自相关序列计算平均峰值距离");

subplot(7, 1, 4);
plot(p);
hold on;
[top_p,loc_p] = findpeaks(p, 'minpeakdistance',MPD/2);
% top_v = -top_v;
scatter(loc_p,top_p,'p','r','filled');loc_op1=loc_p;top_op1=top_p;
title("以平均峰值距离计算峰值点");


%% second step,zero-crossing,所有小于0的都不要了
subplot(7, 1, 5);
loc_op=[];
top_op=[];
for i = 1:length(loc_p)
    if(top_p(i)>0)
        loc_op = [loc_op loc_p(i)];
        top_op = [top_op top_p(i)];
    end
end

plot(p);
hold on;
scatter(loc_op,top_op);loc_op2=loc_op;top_op2=top_op;
title("过零检测，去除所有小于0的点")

%% third step,计算平均峰值，去除小于平均峰值的点
avg = mean(top_op);  % 计算所有 峰值 的平均值
for i = 1:length(loc_op)
    temp = top_op(i) - avg;
    if (top_op(i) < avg/1.5)  % 防止错漏，保证正确性
        loc_op(i) = nan;
        top_op(i) = nan;
    end
end
loc_op(isnan(loc_op))=[];  
top_op(isnan(top_op))=[]; 
subplot(7, 1, 6);
plot(p);
hold on;
scatter(loc_op,top_op);loc_op3=loc_op;top_op3=top_op;
title("去除所有小于平均峰值d的点");


%% MDP judge, 根据平均峰值距离判定，相邻的两个点若距离小于平均峰值距离。则取更大点保留。
% 需要多加一个判断是否两点的差值<MPD/1.25，因为当C比较小的时候，是不需要平均峰值这一步的
for i = 1:length(loc_op)-1
    if(loc_op(i+1)-loc_op(i)<MPD/uptraParameter1)  % uptraParameter1!!
        if(top_op(i+1)>top_op(i))
            top_op(i) = nan;
            loc_op(i) = nan;
        else
            top_op(i+1) = nan;
            loc_op(i+1) = nan;
        end
    end
end
% 删除数组中的nan
loc_op(isnan(loc_op))=[]; 
top_op(isnan(top_op))=[]; 

%% 
subplot(7, 1, 7);
plot(p);
hold on;
scatter(loc_op,top_op);loc_op4=loc_op;top_op4=top_op;
plot(direc);
title("相邻的两个点若距离小于平均峰值。则取更大点保留");


%% （去直流）
% dc1有效
% [p] = SMI_API_ELIMINATE_DC1(p_init,direction,"time");
%     [p] = SMI_API_ELIMINATE_DC2(p_init,direction,"time");
%     [p1,p2] = SMI_API_evenlopEXTRACT_HT_PRO(p,N);
%     p = p2;
% figure(2);
% subplot(6,1,1)
% plot(p);
% title("去直流后的自混合信号");

%% 时频分析
figure(3);
[T,F,TF,TF_curb,p] = SMI_API_TFPM(p,N,fs,windowLength,V);
mesh(T,F,abs(TF_curb)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' ); 
view(0,90); % 设置初始视角为俯视角

%% 归一化(p = 2 * (p - -abs(hilbert(-p)))./(abs(hilbert(p)) - -abs(hilbert(-p))) - 1;)
figure(2);
% p = 2 * (p - -abs(hilbert(-p)))./(abs(hilbert(p)) - -abs(hilbert(-p))) - 1;  % 归一化
p = sgolayfilt(p,1,11);
p = sgolayfilt(p,2,21);
% p = sgolayfilt(p,3,31);

[top_p,loc_p] = findpeaks(p);  % 找到所有大于0.1的峰值,这是因为驼峰区乱七八糟（当然这只针对该实验信号）
[top_v,loc_v] = findpeaks(-p);  % 默认。 'minpeakheight',0.1, 'minpeakdistance',50
top_v = -top_v;
loc_ov = loc_v;

subplot(6,1,2)
plot(p);
% hold on;
% plot(direc);
hold on;
scatter(loc_p,top_p);loc_op5=loc_p;top_op5=top_p;loc_ov0 = loc_v;top_ov0=top_v;
scatter(loc_v,top_v);
hold on;
title("降C后的自混合信号峰谷值 及 原始信号方向");

%% 对比降C后的信号与带噪信号峰值的位置
dir = zeros(1,N);
MPD2 = mean(diff(loc_v));  % 降C后信号的平均峰值距离，据此找驼峰区 

hump_seg = [];
for i = 1:length(loc_v)-1
    % loc_op是原始带噪信号的峰值位置，loc_v是降C后信号的谷值点，在此区间内找峰值！
    before_idx = find(loc_op>loc_v(i) & loc_op<loc_v(i+1));  % 找到原始带噪信号的峰值 位置
    after_idx = find(loc_p>loc_v(i) & loc_p<loc_v(i+1));  % 找到降C后的信号的峰值 位置
    
    
    if isempty(before_idx) || isempty(after_idx)  % ✔这是核心，因为降C前后的信号，驼峰区峰谷值都可能是乱的（多了或者少了），而其他区域一般都是对的，所以当判定是乱的，直接视为驼峰区！
                                                  % 当我们在降C后的信号的谷值区间找对应峰值时，只要找不到，直接置零
        dir(loc_v(i):loc_v(i+1)) = 0;
        hump_seg = [hump_seg,loc_v(i),loc_v(i+1)];
    elseif(loc_v(i+1)-loc_v(i)>MPD2*uptraParameter2)  % ✔ % uptraParameter2!!驼峰区置零，超参数！，通常设置为MPD2*1.2。
        dir(loc_v(i):loc_v(i+1)) = 0;
        hump_seg = [hump_seg,loc_v(i),loc_v(i+1)];  % 需要记录位置，同时这个地方存在一定的鲁棒性，即通常情况下，驼峰区通常是对称大于平均峰值距离的~~~通过方向矫正，都能得到正确的方向               
    elseif(loc_op(before_idx) > loc_p(after_idx))  % ✔比较降C前后的峰值点偏离方向
        dir(loc_v(i):loc_v(i+1)) = 1;
    else
        dir(loc_v(i):loc_v(i+1)) = -1;
    end
end
plot(dir);tempDir = dir;

%% 处理驼峰区 
N = length(dir);
for i=1:2:length(hump_seg)
    pre = dir(hump_seg(i)-1);  % 记录前值
    len = (hump_seg(i+1)-hump_seg(i))/2;
    n = floor(len);
    
    dir(hump_seg(i):hump_seg(i)+n) = pre;  % 将前半部分设为前值
    dir(hump_seg(i)+n:hump_seg(i+1)) = -pre;  % 后半部分设为前值取反
end

% i = 2;
% while i<N
%     if(dir(i)==0)
%         pre = dir(i-1);  % 记录前值
%         % j = i + 1;
%         for j = i+1:N
%             if(dir(j)~=0)
%                 break
%             end
%         end
%         n = floor((j-i)/2);
%         
%         % 将0的前半部分设为前值
%         for k = i:i+n-1
%             dir(k) = pre;
%         end
% 
%         % 将0的后半部分设为前值取反
%         for k = i+n : j-1
%             dir(k) = -pre;
%         end
%         i = i + n -1;
%     end
%     i = i+1;
% end

% 处理头和尾
dir(1:loc_v(1)) = dir(loc_v(1));
dir(loc_v(length(loc_v)):end) = dir(loc_v(length(loc_v))-1);


figure(2);
subplot(6,1,3);
plot(p);
hold on;
plot(dir);ddir = dir;

% dir = direc;
%% hilbert transform
% subplot(5, 1, 5);
inverse_hb = (imag(hilbert(p))) .* dir;  % 得到sin 
% plot(inverse_hb,"b");
% hold on;
% plot(p,'r');
% title("blue,inverse-HT and orginal-HT","red,self-mixing signal");

% subplot(5, 1, 5);
phiF_wrapped = atan(inverse_hb./p);
% plot(phiF_wrapped);
% hold on;
% title("phiF_wrapped")

%% arctan相位解包裹
for i = 2:N
    if (phiF_wrapped(i) - phiF_wrapped(i-1) > pi/2)
        phiF_wrapped(i:end) = phiF_wrapped(i:end) -  pi;
    elseif (phiF_wrapped(i) - phiF_wrapped(i-1) < -pi/2)  % 找峰值点，碰到峰值则加pi
        phiF_wrapped(i:end) = phiF_wrapped(i:end) +  pi;            
    end
end

phiF_reconstruct = phiF_wrapped;

%% 参数估算
[C_reconstruct, alpha_reconstruct] = SMI_API_ESTIMATE_C(phiF_reconstruct);  

%% 重构
subplot(6, 1, 4);
phi0_reconstrut = phiF_reconstruct + C_reconstruct*sin(phiF_reconstruct+atan(alpha));
Lt_reconstruct = phi0_reconstrut * lambda / (4 * pi);
plot(Lt,'k')
hold on;

Lt_reconstruct = Lt_reconstruct - mean(Lt_reconstruct); % 简写振动和余弦调制振动，重构后加上幅值A
% Lt_reconstruct = Lt_reconstruct + 2 * lambda;  % 重构后的随机振动信号要加上幅值1.5的波长，这是为啥我页不知道

plot(Lt_reconstruct,'r')
title("重构后的信号");







% fringeData = [];
% for i=2:length(loc_ov)
%     % 得到谷值区间点数（条纹）
%     N = loc_ov(i) - loc_ov(i-1);
%     % 拿到谷值区间内所有的点
%     dd = loc_ov(i-1)+1:loc_ov(i)-1;
%     % 使用众数记录当前区间方向
%     m = mode(direction(loc_ov(i-1):loc_ov(i)));
%     % 判断谷值区间 内！，是否含有翻转点，如果有则 if=false
%     if isempty(intersect(dd, loc_r))
%         % 不包含，记录当前方向
%         dir = m;
%     else
%         % 包含,则方向置零
%         dir = 0;
%     end
%     loc = p(loc_ov(i-1):loc_ov(i));
%     int_ = 1000;
%     loc_ = [SMI_API_RESAMPLE(loc,int_) dir];  % 方向作为最后一个点
%     fringeData = [fringeData;loc_];
% end

