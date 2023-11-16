clc;
clear all;
close all;
%% 产生自混合信号
% subplot(5, 1, 1);
% fs = 200000;  % 采样率，即fs(s)采一个点。
% N = 8000;  
% fv = 50;  % 震动频
% C = [0.05];  % C设置一个从a到b变化的值
alpha = 5;
% [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_HARMONIC(fs, N, fv, C, alpha);  % 1 简谐振动的自混合信号
% plot(Lt);
% title(['外部简谐振动,C=',num2str(C)]);
% subplot(5, 1, 2);
% % p = awgn(p,40);  % 10db，加高斯白噪声
% % p = p .* (1+0.2*cos(2*pi*75*t));  % 给自混合信号加包络，加了一个幅值为0.2，频率为75的包络
% plot(p);
% title("自混合信号");

%% 实验信号判断
subplot(5, 1, 1);
path =  'D:\matlab save\smi_实验信号\L3_moderate_38332_43331.csv';  % 5 文件路径, M/N/win/w = 250515/16000/128/150
M = 40000; N = 10000; lambda = 650e-9; [t, p, fs] = MOVE_API_EXPERIMENT(M, N, path);  % 5 从M点处取N个点
plot(p);

%% TFPM信号的方向判断
% subplot(5, 1, 1);
% load("TFPM.mat");
% N=4000;
% plot(p);


%% 得到重构所需的相关信息
subplot(5, 1, 2); 
p = sgolayfilt(p,1,11);  % 第二个参数是阶数，通常123就够啦。第三个参数是长度，长度越短，滤除的就是越高频的分量，因为噪声是高频分量，所以一般不要太大！！！
p = sgolayfilt(p,2,21);
p = sgolayfilt(p,3,31);
[top_v,loc_v] = findpeaks(-p);
top_v = -top_v;
plot(p);
hold on;

% arr = [8678,11149,11170,13680];
% arr = [3760,8756,13745];
arr = [1526,3528];
[loc_ov,top_ov] = ditchReversePoing2(arr,loc_v,top_v);
scatter(loc_ov,top_ov);

%% 基于神经网络的方向判断（注意要保证loc_v的连续性）
load("matlab2.mat")
dir = zeros(1,N);
loc_r = [];
rect = [];
int_ = 1000;  % 插值倍数(重采样放大倍数，默认为1000)
fs = int_;  % 因为原本采样率为1，放大后采样率就是int_
for i = 2:length(loc_ov)
    N = loc_ov(i)-loc_ov(i-1);
    % 将原信号重采样为原来的int_倍
    p_ = interp1(1:int_:int_* N,p(loc_ov(i-1):loc_ov(i)-1),1:int_* N,'spline');
    % 再将信号降采样为原来的N倍，即 N * 1000 / N， 故最后信号长度为int_
    p_ = p_(1:(loc_ov(i)-loc_ov(i-1)):length(p_));
    judge = DirModel2.predictFcn(p_);  % 基于神经网络判断的方向
    
    if(i==2)
        dir(1:loc_ov(i)) =  judge;
    elseif(i==length(loc_ov))   
        dir(loc_ov(i-1):end) =  judge;
    else
        dir(loc_ov(i-1):loc_ov(i)) =  judge;
    end
end
subplot(5, 1, 3);
plot(p);
hold on;
plot(dir);
title("初始判断的方向");

%%  
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
subplot(5, 1, 4);
plot(p);
hold on;
plot(dir);
title("初始判断的方向");


%% hilbert transform
% subplot(5, 1, 5);
inverse_hb = (imag(hilbert(p))) .* dir;  % 得到sin 
plot(inverse_hb,"b");
% hold on;
% plot(p,'r');
% title("blue,inverse-HT and orginal-HT","red,self-mixing signal");

subplot(5, 1, 5);
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
subplot(5, 1, 5);
phi0_reconstrut = phiF_reconstruct + C_reconstruct*sin(phiF_reconstruct+atan(alpha));
Lt_reconstruct = phi0_reconstrut * lambda / (4 * pi);
% plot(Lt,'k')
% hold on;

Lt_reconstruct = Lt_reconstruct - mean(Lt_reconstruct); % 简写振动和余弦调制振动，重构后加上幅值A
% Lt_reconstruct = Lt_reconstruct + 1.5 * lambda;  % 重构后的随机振动信号要加上幅值1.5的波长，这是为啥我页不知道

plot(Lt_reconstruct,'r')
title("重构后的信号");




















%% subfunction ditch reversepoint，输入不想要的点（翻转点）的数组，在峰谷值中去掉翻转点
function [loc_p,loc_v,top_p,top_v] = ditchReversePoing(arr,loc_p,loc_v,top_p,top_v)
    for i=1:length(arr)
        aa = find(loc_p==arr(i));
        bb = find(loc_v==arr(i));
        flag1 = isempty(aa);
        flag2 = isempty(bb);
        if flag1==0
            loc_p(aa)=[];
            top_p(aa)=[];
        end
        if flag2==0
            loc_v(bb)=[];
            top_v(bb)=[];
        end
    end
end

%% subfunction ditch reversepoint，输入不想要的点谷值点的数组，在谷值点中去掉翻转点
function [loc_v,top_v] = ditchReversePoing2(arr,loc_v,top_v)
    for i=1:length(arr)
        bb = find(loc_v==arr(i));
        flag2 = isempty(bb);
        if flag2==0
            loc_v(bb)=[];
            top_v(bb)=[];
        end
    end
end

%% subfunction direction，输入翻转点arr数组(顺序输入) （方向有误乘-1即可）
function direction = solve_dir(arr,N)
    arr1 = [1, sort([arr, arr+1]), N];
    temp = zeros(1,N);
    flag = 1;
    
    for i = 1:2:length(arr1)-1      
        temp(arr1(i):arr1(i+1)) = flag;
        flag = flag * -1;
    end
    direction = sign(temp);
end












