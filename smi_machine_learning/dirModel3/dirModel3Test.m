%% 产生自混合信号
subplot(5, 1, 1);
fs = 200000;  % 采样率，即fs(s)采一个点。
N = 8000;  
fv = 50;  % 震动频
C = [0.05];  % C设置一个从a到b变化的值
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
% scatter(loc_ov,top_ov);
scatter(loc_p,top_p);
scatter(loc_v,top_v);
scatter(loc_r,top_r);
plot(direction);
title("pks, vls, rev")

%% 基于神经网络的方向判断（注意要保证loc_v的连续性）
load("matlab2.mat")
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
















