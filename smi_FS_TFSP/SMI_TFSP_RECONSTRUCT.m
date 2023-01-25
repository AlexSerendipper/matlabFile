%% FS/TFSP就是对phiF做傅里叶变换，抑制后逆傅里叶变换出来的值近似为phi0。 
% 而TFPM就是对p做傅里叶变换，抑制后逆傅里叶变换出来的值近似于C=0！！！

clc;
clear all;
close all;

%% 产生自混合信号
subplot(6, 1, 1);
fs = 200000;  % 采样率，即fs(s)采一个点
N = 4000;  
fv = 300;  % 震动频
C = 2;
alpha = 5;
[t, lambda, L0, Lt, phi0, p] = SMI_API_HARMONIC(fs, N, fv, C, alpha);  % 简谐振动的自混合信号
% [t, lambda, L0, Lt, phi0, p] = SMI_API_SIN(fs, N, C, alpha);  % 正弦调制信号的自混合信号

% cut = 500;  % cut降采样，输入一个能被N整除的数，将N分为N/cut段
% [t, lambda, L0, Lt, phi0, p] = SMI_API_ALEATORY(fs, N, cut, C, alpha);  % 使用随机振动时，方向×负

plot(Lt);
title(['外部简谐振动,C=',num2str(C)]);
subplot(6, 1, 2);
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
subplot(6, 1, 3);
phiF_reconstruct = phiF_reconstruct - mean(phiF_reconstruct);  % 消除直流分量
f = fs / N * (0 : 1 : N-1);  % Fs/N就是这个频谱中的最小频率间隔！！！！！所以N越大，分辨率会越高
phiF_reconstruct_ = fft(phiF_reconstruct, N);
Fv = phiF_reconstruct_;  % 用来储存傅里叶变换后的值
Av = abs(phiF_reconstruct_) * 2 / N;
plot(f(1:N/2), Av(1:N/2));  % fft算法默认是双边谱,通常我们只取一半
title('PhiF的幅度谱');

%% TFSP重构，说实话我感觉这就是FS简化版，我都没看出来有用到阈值，那篇文章就是包装了一下把
subplot(6, 1, 4);
plot(phiF_reconstruct);
hold on;
peakst = find( Av == max(Av) );  % 找到最大值的索引，双边谱所以两个
% powerf(peakst) = min(powerf);  % 将最小值赋值给最大值处
rebuild = zeros(1,N);
rebuild(peakst) = Fv(peakst);  % 最大值索引处的Fv值，赋值给rebuild

rebuildtime = ifft(rebuild, N);

plot(rebuildtime);
title("标准的phiF and 近似的phi0")


subplot(6,1,5);
phi0_reconstrut = rebuildtime;
Lt_reconstruct = phi0_reconstrut * lambda / (4 * pi);

plot(Lt_reconstruct);
title("重构后的信号");


hold on;
plot(Lt);

subplot(6,1,6);
plot(Lt-Lt_reconstruct);
RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
title(['绝对误差，RMSE=', num2str(RMSE)])

%% SP(这里的SP是稀疏性提升)
% figure(3);
% subplot(3,1,1)
% plot(phi0_reconstrut);
% subplot(3,1,2)
% gg = translate(phi0_reconstrut,1000);
% plot(gg);
% subplot(3,1,3)
% ggg = translate(phi0_reconstrut,-250);
% plot(ggg);
argmin = 10000000;  % 找最小值，所以初值设置大一点
xf = phiF_reconstruct;
% j = 1.4056;
% j = 0.9;
tempp = [];
iindex = [];
for i = -20:1:20
    x0 = phi0_reconstrut;  % 主播
    for k = 0.5:0.01:1.5
        x0 = phi0_reconstrut;  % 主播
        for j = 1.2490:0.01:1.4601  % 这里是α的范围，这里已经包含了arctan(3)~arctan(6)，注意这里的步长必须得是0.01，有点夸张
            x0 = translate(x0, i);  % 每一次我平移,  
            x0 = k * x0;
            Ck = (x0 - xf) ./ (sin(xf + j));  % 这里我暂时遍历alpha
            Ck = median(Ck);
            x0_ = xf + Ck .* sin(xf + j);
            temp = norm(fft(x0_),1);
    %       temp = sum( diff(x0_).^4 );

            if temp < argmin
                argmin = temp;
                CC_index = i;
                CC = Ck;
                AA = tan(j);
                tempp = [tempp, temp];
                iindex = [iindex, i];
                KK = k;
            end
        end
    end
end

figure(3);
subplot(3,1,3);
plot(phi0_reconstrut,"k");
phi0_reconstrut = KK * translate(phi0_reconstrut,CC_index);

hold on;
plot(phi0_reconstrut);


% 重新赋值
% Ck = CC;
% argmin = 1000000000;
% for i = -40:1:10
%     x0 = phi0_reconstrut;
%     x0 = x0 + Ck  - i;
%     Uk = 1/2 .* ((x0 - xf) .* cos(xf) - sin(xf)) .* sqrt(abs(Ck.^2 - (x0 - xf).^2));
%     Uk = median(Uk);
%     a = Uk / sqrt(abs(1-Uk^2));
%     
%     x0_ = xf + Ck .* sin(xf + a);
% %   temp = sum( diff(x0_) .^2 );
%     temp = norm(fft(x0_),1);
% %     if temp < argmin
%         argmin = temp;
%         index = i;
%         AA = [AA, a];
% %     end
% end




%% 
figure(3);
subplot(3,1,1);
Lt_reconstruct = phi0_reconstrut * lambda / (4 * pi);
plot(Lt_reconstruct);
title("重构后的信号");
hold on;
plot(Lt);




subplot(3,1,2);
plot(Lt-Lt_reconstruct);
RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
title(['绝对误差，RMSE=', num2str(RMSE)])





%% translate 仅适用于周期信号平移操作,x<0向左平移，x>0向右平移
function pp = translate(p, x)  % p为平移的信号，x为平移的距离
    L = length(p);
    if(x < 0)  % 向左平移
        x = -x;
        p(L+1:L+x) = p(1:x);
        pp = p(x+1:L+x);
    elseif( x > 0 )  % 向右平移
        temp(x+1:L+x) = p(1:L);  % 整体右移
        temp(1:x) = p(L-x+1:L);
        pp = temp(1:L);
    else
        pp = p;
    end
    
end



%% scaling