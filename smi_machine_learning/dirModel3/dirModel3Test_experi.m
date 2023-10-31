clc;
clear all;
close all;

%% 实验信号
path =  'D:\matlab save\smi_实验信号\experSignal\3_m_38332_5000.csv';  % 5 文件路径, M/N/win/w = 250515/16000/128/150
% path =  'D:\matlab save\smi_实验信号\L1_weak_365424_380423.csv';  % 5 文件路径, M/N/win/w = 250515/16000/128/150
% path =  'D:\matlab save\smi_实验信号\f(100)_A(3um).csv';
% path =  ' D:\matlab save\smi_实验信号\大论文数据_GSZ\TEK00012.CSV';
M = 35000; N = 5000; lambda = 650e-9; [t, p, fs] = MOVE_API_EXPERIMENT(M, N, path);  % 5 从M点处取N个点
subplot(5,1,1);
plot(p);
title("原始信号");
% direction = solve_dir([1519,4012,6541,9027],N);
% direction = -direction;
% （去直流）
%     [p] = SMI_API_ELIMINATE_DC1(p,direction,"time");
%     [p] = SMI_API_ELIMINATE_DC2(p_init,direction,"time");
%     [p1,p2] = SMI_API_evenlopEXTRACT_HT_PRO(p,N);
%     p = p2;

% 全局变量
windowLength = 512; % 窗长
V = 0.45; % 抑制因子

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

subplot(5,1,2)
plot(p);
hold on;
% plot(direction);
scatter(loc_p,top_p);
scatter(loc_v,top_v);
title(['自混合信号，C', num2str("≈0")]);


%% 得到重构所需的相关信息
% subplot(5, 1, 1); 
% p = sgolayfilt(p,1,11);  % 第二个参数是阶数，通常123就够啦。第三个参数是长度，长度越短，滤除的就是越高频的分量，因为噪声是高频分量，所以一般不要太大！！！
% p = sgolayfilt(p,2,21);
% p = sgolayfilt(p,3,31);
% [top_v,loc_v] = findpeaks(-p);
% top_v = -top_v;
% plot(p);
% hold on;
% 
% % arr = [8678,11149,11170,13680];
% % arr = [3760,8756,13745];
% arr = [1526,3528];
% [loc_ov,top_ov] = ditchReversePoing2(arr,loc_v,top_v);
% scatter(loc_ov,top_ov);

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
subplot(5, 1, 3);
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












