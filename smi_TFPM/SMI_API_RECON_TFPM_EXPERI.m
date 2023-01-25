%% 对于实验信号，需要先进行时频谱抑制才能判断方向~ 所以步骤与仿真实验不同
clc;
clear all;
close all;
set(0,'defaultfigurecolor','w');  % 设置画布背景色为白色

%% 全局变量
windowLength = 128; % 窗长
overlapLength = floor(windowLength * 0.9);  % OverlapLength后为指定的重叠长度
window = hamming(windowLength, "periodic");  % 使用汉明窗作为滑动的窗口
fftLength = 5*windowLength;
V = 0.65; % 抑制因子2

%% 产生自混合信号
fs = 200000;  % 采样率，即1/fs(s)采一个点。
N = 4000;  
fv = 200;  % 震动频率
C = [1.5];  % C设置一个从a到b变化的值
alpha = 4;

% [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_HARMONIC(fs, N, fv, C, alpha);  % 1 简谐振动的自混合信号

% [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_SIN_MODULATE(fs, N, C, alpha);  % 2 正弦调制信号的自混合信号

% cut = 200; [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_ALEATORY(fs, N, cut, C, alpha);  % 3 产生随机振动时，方向×负，输入一个能被N整除的数，将N分为N/cut段

% [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_ALEATORY_LOAD(fs, N, C, alpha);  % 4 加载储存好的随机振动，方向×负

path =  'D:\matlab save\self-mixing\smi_脉搏波\scope_1.csv';  % 5 文件路径, M/N/win/w = 250515/16000/128/150
M = 30000; N = 30000; lambda = 650e-9; [t, p, fs] = MOVE_API_EXPERIMENT(M, N, path);  % 5 从M点处取N个点

p_init = p;
% p = awgn(p,10); % 10db，加噪声在抑制是可以的，问题是加了噪声没法从原信号判断方向了
figure(1);
set(gcf,'position',[15, 30, 800, 800]);
subplot(7,2,1);
plot(p); ylabel('p(t)(nm)','FontSize',7,'FontWeight','bold'); xlabel('sample(num)','FontSize',7,'FontWeight','bold');  % ✔这玩意赶紧学 
title(['自混合信号,C=',num2str(C), 'alpha=',num2str(alpha)]);
% p = p - mean(p);
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
p = [zeros(1,padding) p zeros(1,padding)];

%% 时频分析（抑制）
TF = stft(p, fs, 'Window', window, 'OverlapLength', overlapLength, 'FFTLength', fftLength);  % （wf每一列就是每一个时间维度的频率）
weight1 = abs(TF)./max(abs(TF));  % 要把这个矩阵当作权值使用，所以按列先归一化
weight2 = weight1;
weight2( find(weight2 < V) ) = 0; % 算法2：好用
weight2 = weight2 ./ max(weight2);  % 抑制后按列归一化
TF_curb = TF .* weight2;  % 获得抑制后的信号

size = size(TF);
[x,y] = meshgrid(1:size(2),1:size(1));
subplot(7,2,[2,4,6]);
mesh(x,y,abs(TF)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' ); 
% mesh(x,y,angle(TF));
view(0,90)
title('抑制前')

subplot(7,2,[10,12,14]);
mesh(x,y,abs(TF_curb)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' );
% mesh(x,y,angle(TF_curb));
view(0,90);
title('抑制后');


%%
figure(2);
subplot(7,2,[2,4,6]);
% [a,b] = find(TF_curb ~= 0);
% TF(a,:) = 0;
% TF_curb(find(TF ~= 0)) = 0;
TTF = TF-TF_curb;
mesh(x,y,abs(TTF));
view(0,90)
% mesh(x,y,abs(TF_curb - TF));




weight1 = abs(TTF)./max(abs(TTF));  % 要把这个矩阵当作权值使用，所以按列先归一化
weight2 = weight1;
weight2( find(weight2 < V) ) = 0; % 算法2：好用
weight2 = weight2 ./ max(weight2);  % 抑制后按列归一化
TF_curb = TTF .* weight2;  % 获得抑制后的信号
subplot(7,2,[10,12,14]);
mesh(x,y,abs(TTF-3*TF_curb)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' );

%% 某时刻频谱分析
% figure(2);
% set(gcf,'position',[30, 30, 243/1.25, 258/1.25]);  % 画布真的对应figure中的画布区域，设置值为测量值÷1.25
% subplot(2,1,1);
% plot(1:size(1), abs(TF(:,20)));
% subplot(2,1,2);
% plot(1:size(1), abs(TF_curb(:,20)))
% % plot(1:size(1)/2, abs(TF_curb(1:size(1)/2 , 80)));

%% 去包络


%% 归一化
p = (istft(TF_curb, fs, 'Window', window, 'OverlapLength', overlapLength, 'FFTLength', 5*windowLength))';
p_init2 = p;  % 用p_init2存储未归一化的信号
p_init2 = p_init2(padding+1:end-padding);  % ✔去掉之前的padding
p_init2 = sgolayfilt(p_init2,4,31);  % ✔ 超参数，平滑滤波器, 这里的参数也需要细调！！！c小的时候可以设置4，31   c大的时候可以设置，2，21 

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


%% 得到重构所需的方向信息
% [top_p, loc_p, top_v, loc_v, top_r, loc_r, direction] = SMI_API_FRINGE(p_init,N);  % 法1求方向，求的是初始信号的方向，所以都能用
W = 150; direction = SMI_API_TFPM_DIRECTION(p,N,W);  % 法2求方向，利用平均峰值距离（要求均为w形的驼峰），且必须用归一化2的信号
% direction = -direction;  % 如果初始震动用的cos，那方向是反的，一定要注意
plot(direction);

%% （去直流）
% [p] = SMI_API_ELIMINATE_DC(p,direction);
% figure(1);
% set(gcf,'position',[15, 30, 800, 800]);
% subplot(7,2,1);
% plot(p); ylabel('p(t)(nm)','FontSize',7,'FontWeight','bold'); xlabel('sample(num)','FontSize',7,'FontWeight','bold');  % ✔这玩意赶紧学 
% title(['自混合信号,C=',num2str(C), 'alpha=',num2str(alpha)]);

%% 希尔伯特变换重构
% inverse_hb = (imag(hilbert(p_init2))) .* direction;  % 得到sin ，使用未归一化的信号，精度更高
% subplot(7, 2, 9);
% phiF_wrapped = atan(inverse_hb./p_init2);  %  使用未归一化的信号，精度更高
% plot(phiF_wrapped);
% hold on;
% title("phiF-wrapped")

inverse_hb = (imag(hilbert(p))) .* direction;  % 得到sin ，使用未归一化的信号，精度更高
subplot(7, 2, 9);
phiF_wrapped = atan(inverse_hb./p);  %  使用未归一化的信号，精度更高
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
Lt = (3*10^-6) .* sin(2*pi*100*t);  % 葛兄采的数据是这样
plot(Lt,'k')
hold on;
plot(Lt_reconstruct,'r'); ylabel('Disp(um)','FontSize',7,'FontWeight','bold'); xlabel('sample(num)','FontSize',7,'FontWeight','bold');
title(['重构后的信号，C_reconstruct=', num2str(C_reconstruct)])

%% 误差分析
subplot(7,2,13);
plot((Lt-Lt_reconstruct)*10^9); ylabel('error(nm)','FontSize',7,'FontWeight','bold'); xlabel('sample(num)','FontSize',7,'FontWeight','bold');
RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
title(['绝对误差，RMSE=', num2str(RMSE)])

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


