clc;
clear all;
close all;
set(0,'defaultfigurecolor','w');  % 设置画布背景色为白色

%% 全局变量
windowLength = 256; % 窗长
overlapLength = floor(windowLength * 0.9);  % OverlapLength后为指定的重叠长度
window = hamming(windowLength, "periodic");  % 使用汉明窗作为滑动的窗口
fftLength = 5*windowLength;
V = 0.65; % 抑制因子2

%% 产生自混合信号
fs = 200000;  % 采样率，即1/fs(s)采一个点。
N = 4000;  
fv = 100;  % 震动频率
C = [2];  % C设置一个从a到b变化的值
alpha = 4;

[t, lambda, L0, Lt, phi0, p, c] = MOVE_API_HARMONIC(fs, N, fv, C, alpha);  % 1 简谐振动的自混合信号

% [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_SIN_MODULATE(fs, N, C, alpha);  % 2 正弦调制信号的自混合信号

% cut = 200; [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_ALEATORY(fs, N, cut, C, alpha);  % 3 产生随机振动时，方向×负，输入一个能被N整除的数，将N分为N/cut段

% [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_ALEATORY_LOAD(fs, N, C, alpha);  % 4 加载储存好的随机振动，方向×负

% path =  'D:\matlab save\self-mixing\smi_实验信号\gexiong_f(100)_A(3um).csv';  % 5 文件路径
% M = 62000; N = 8000; lambda = 650e-9; [t, p, fs] = MOVE_API_EXPERIMENT(M, N, path);  % 5 从M点处取N个点 M60000/62000   N8000  E69 E20!!!!!!!!!

p_init = p;
% p = awgn(p,40); % 10db，加噪声在抑制是可以的，问题是加了噪声没法从原信号判断方向了
% p = p .* (1+3*cos(2*pi*100*t));  % 给自混合信号加包络，加了一个幅值为0.2，频率为75的包络

%% 得到重构所需的方向信息（实验不需要）
[top_p, loc_p, top_v, loc_v, top_r, loc_r, direction] = SMI_API_FRINGE(p_init,N);  % 法1求方向，求的是初始信号的方向，所以都能用
% W = 80; direction = SMI_API_TFPM_DIRECTION(p,N,W);  % 法2求方向，利用平均峰值距离（要求均为w形的驼峰），且必须用归一化2的信号
% direction = -direction;  % 如果初始震动用的cos，那方向是反的，一定要注意

%% 实验信号需要
% [p, top_p,loc_p,top_v,loc_v,direction,        diffp, top_diffp_p,loc_diffp_p,top_diffp_v,loc_diffp_v,en_top,en_bottom,en_median ] = SMI_API_FRINGE_EXPERI(p,N,0.5,10,0,10);  % 实验信号所需1
% [Lt] = SMI_API_RECON_PUM_EXPERI(p,N,lambda,0.5,10,0,10);  % 实验信号所需2

%% （去直流）
[p] = SMI_API_ELIMINATE_DC(p,direction);
figure(1);
set(gcf,'position',[15, 30, 800, 800]);
subplot(7,2,1);
plot(p); ylabel('p(t)(nm)','FontSize',7,'FontWeight','bold'); xlabel('sample(num)','FontSize',7,'FontWeight','bold');  % ✔这玩意赶紧学 
title(['自混合信号,C=',num2str(C), 'alpha=',num2str(alpha)]);

%% padding,逆变换不可避免的会让信号变短，要padding
judge = true; 
padding = windowLength;  % 逆变换不可避免的会让信号变短，要padding
while judge
    if mod(((2 * padding + N - overlapLength) / (windowLength - overlapLength)), 1) == 0
        judge = false;
    else
        padding = padding - 1;
    end
end
p = wextend('1D','sym',p,padding);  % 数据延拓
% p = [zeros(1,padding) p zeros(1,padding)];

%% 时频分析（抑制）
[TF,F,T] = stft(p, fs, 'Window', window, 'OverlapLength', overlapLength, 'FFTLength', fftLength);  % （wf每一列就是每一个时间维度的频率）
weight1 = abs(TF)./max(abs(TF));  % 要把这个矩阵当作权值使用，所以按列先归一化

TF_curb = TF_inhibit1(TF,V);

subplot(7,2,[2,4,6]);
mesh(abs(TF)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' ); 
% view(0,90); % 设置初始视角为俯视角
% mesh(angle(TF));
title('抑制前');

subplot(7,2,[10,12,14]);
mesh(abs(TF_curb)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' );
% mesh(angle(TF_curb));
% view(0,90);
title('抑制后');

%% 某时刻频谱分析
figure(2);
set(gcf,'position',[30, 30, 243/1.25, 258/1.25]);  % 画布真的对应figure中的画布区域，设置值为测量值÷1.25
subplot(2,1,1);
plot(1:length(F), abs(TF(:,20)));
subplot(2,1,2);
plot(1:length(F), abs(TF_curb(:,20)))
% plot(1:size(1)/2, abs(TF_curb(1:size(1)/2 , 80)));

%% 归一化
p = (istft(TF_curb, fs, 'Window', window, 'OverlapLength', overlapLength, 'FFTLength', 5*windowLength))';
p_init2 = p;  % 用p_init2存储未归一化的信号
p_init2 = p_init2(padding+1:end-padding);  % ✔去掉之前的padding
% p_init2 = sgolayfilt(p_init2,4,31);  % ✔ 超参数，平滑滤波器, 这里的参数也需要细调！！！c小的时候可以设置4，31   c大的时候可以设置，2，21 

figure(1);
subplot(7,2,3);
plot(p_init2);
title('未归一化');

% subplot(7,2,5);
% plot(c);
% legend("Variation of C");

subplot(7,2,7);
p = 2 * (p - -abs(hilbert(-p)))./(abs(hilbert(p)) - -abs(hilbert(-p))) - 1;
p = p(padding+1:end-padding);  % ✔去掉之前的padding
p = sgolayfilt(p,2,21);  % ✔ 超参数，平滑滤波器, 这里的参数也需要细调！！！c小的时候可以设置4，31   c大的时候可以设置，2，21 
plot(p); ylabel('p(t)(nm)','FontSize',7,'FontWeight','bold'); xlabel('sample(num)','FontSize',7,'FontWeight','bold');
title('归一化2');
hold on;
% plot(direction);

%% 希尔伯特变换重构
subplot(7, 2, 9);
inverse_hb = (imag(hilbert(p_init2))) .* direction;  % 得到sin ，使用未归一化的信号，精度更高
phiF_wrapped = atan(inverse_hb./p_init2);  %  使用未归一化的信号，精度更高

% inverse_hb = (imag(hilbert(p))) .* direction;  % 使用归一化的信号
% phiF_wrapped = atan(inverse_hb./p);  %  使用归一化的信号

plot(phiF_wrapped);
hold on;
title("phiF-wrapped")

%% arctan相位解包裹
for i = 2:N
    if (phiF_wrapped(i) - phiF_wrapped(i-1) > pi/2)
        phiF_wrapped(i:end) = phiF_wrapped(i:end) -  pi;
    elseif (phiF_wrapped(i) - phiF_wrapped(i-1) < -pi/2)  % 找峰值点，加pi
        phiF_wrapped(i:end) = phiF_wrapped(i:end) +  pi;            
    end
end

phiF_reconstruct = phiF_wrapped;
% subplot(7, 1, 6);
% plot(phiF_reconstruct)
% title("重构后的phiF")

%% 计算一下抑制后，C的值
[C_reconstruct, alpha_reconstruct] = SMI_API_ESTIMATE_C(phiF_reconstruct); 

%% 重构
subplot(7, 2, 11);
phi0_reconstrut = phiF_reconstruct + C_reconstruct*sin(phiF_reconstruct+atan(alpha));
Lt_reconstruct = phi0_reconstrut * lambda / (4 * pi);


Lt_reconstruct = Lt_reconstruct - mean(Lt_reconstruct); % 简写振动和余弦调制振动，重构后加上幅值A
% Lt_reconstruct = Lt_reconstruct + 1.5 * lambda;  % 重构后的随机振动信号要加上幅值1.5的波长，这是为啥我页不知道
plot(Lt,'k')
hold on;
plot(Lt_reconstruct,'r'); ylabel('Disp(um)','FontSize',7,'FontWeight','bold'); xlabel('sample(num)','FontSize',7,'FontWeight','bold');
title(['重构后的信号，C_reconstruct=', num2str(C_reconstruct)]);

%% 误差分析
subplot(7,2,13);
plot((Lt-Lt_reconstruct)*10^9); ylabel('error(nm)','FontSize',7,'FontWeight','bold'); xlabel('sample(num)','FontSize',7,'FontWeight','bold');
RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
title(['绝对误差，RMSE=', num2str(RMSE)]);

%% 这部分就是我个人的特殊用途了
% a1 = phiF_reconstruct + 0*sin(phiF_reconstruct+atan(alpha));
% b1 = a1 * lambda / (4 * pi);
% b1 = b1 - mean(b1);
% c1 =  sqrt(mean((Lt-b1).^2));
% 
% a2 = phiF_reconstruct + 0.01*sin(phiF_reconstruct+atan(alpha));
% b2 = a2 * lambda / (4 * pi);
% b2 = b2 - mean(b2);
% c2 =  sqrt(mean((Lt-b2).^2));
% 
% a3 = phiF_reconstruct + 0.02*sin(phiF_reconstruct+atan(alpha));
% b3 = a3 * lambda / (4 * pi);
% b3 = b3 - mean(b3);
% c3 =  sqrt(mean((Lt-b3).^2));

%% subfunction
% 时频抑制算法1，嘎嘎好用
function TF_curb = TF_inhibit1(TF,V)
    weight1 = abs(TF)./max(abs(TF));
    weight2 = weight1;
    weight2( find(weight2 < V) ) = 0;  % 算法1：好用
    weight2 = weight2 ./ max(weight2);  % 抑制后按列归一化，用作权值矩阵
    TF_curb = TF .* weight2;  % 获得抑制后的信号
end

% 时频抑制算法2，一般
function TF_curb = TF_inhibit2(TF)
    weight1 = abs(TF)./max(abs(TF));
    weight1(find(weight1 == 1)) = weight1(find(weight1 == 1)) - 1e-11;  % 算法1：是找到weight1中等于1的位置（理论上每列都有），并都减去num极小值，因为这个后边要用做分母，所以减去一个无穷小量，不要让他产生无穷大
    weight2 = (weight1.^64)./(1 - weight1.^64);  % 算法1：构建一个函数，让原来频率分量大的更大，小的更小，windowLenght此处为抑制因子
    weight2 = weight2 ./ max(weight2);  % 抑制后按列归一化，用作权值矩阵
    TF_curb = TF .* weight2;  % 获得抑制后的信号
end
