%% 该程序用于调试，其实也是非常好用且重要的
%% 目前看来，调制最大的好处就是不用判断方向，这还是挺牛的，直接得到PhiF然后解包裹

%% 全局变量
clc;
clear all;
close all;
fs = 10000;  % 采样率，即fs(s)采一个点。
N = 10000;  
fv = 10;  % 震动频
alpha = 5;

%% 产生自混合信号
figure(1);
subplot(9,1,1);
t = (0:N-1)/fs;  % 采样时间，设N=10, fs=200，即采样了0.05s，t为[0...0.045]
lambda = 650e-9;  % 波长
A = 2 * lambda;  % 幅值（✔）
L0 = 20 * lambda;  % 外腔距离（✔） 
Lt = A.* sin(2*pi*fv*t);  % L0为标称位置，外腔长度
beta = 1;  % the amplitude of selfmixing signal
phi0 = 4*pi*(L0+Lt)/lambda;

h = 200;  % 调制深度200
fm = 2000;  % 调制频率2000
gamma = 0;  % 调制信号初始相位
phiM = h * sin(2*pi*fm*t + gamma);  % 调制信号

p = zeros(1,N);
C = [0.2, 0.2];
 % C的变化是一个正弦曲线，不能随机数！
C_lower = C(1);
C_upper = C(2);
% 这个乘和加保证了c的上下限，这里可以设置变换的周期！！但是这个变换周期需要长一点否则会报错！！
x = linspace(0, 3*pi, N);
c = (C_upper-C_lower)/2 * cos(x) + (C_upper - (C_upper-C_lower)/2);
% plot(x,c);

for i = 1:N 
    C = c(i);
    phiF(i) = solve_phiF(C, phi0(i), alpha);
    % p(i) = beta * cos(solve_phiF(C, phi0(i), alpha));  % 遍历所有的phi0
end

p = beta * cos(phiF + 2*phiM);  % 先后两次穿过置于外腔中的 EOM, 因此，实际的相位调制部分为2phiM    
% p = awgn(p,10);  % 10db，加高斯噪声
% p = p .* (1+0.2*cos(2*pi*75*t));  % 给自混合信号加包络，加了一个幅值为0.2，频率为75的包络
plot(p);
hold on;
title("自混合信号");

%% 全局变量
windowLength = 128; % 窗长
overlapLength = floor(windowLength * 0.9);  % OverlapLength后为指定的重叠长度
window = hamming(windowLength, "periodic");  % 使用汉明窗作为滑动的窗口
fftLength = 5*windowLength;  % 每个时刻傅里叶变换的长度
V = 0.65; % 抑制因子2

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

%% 先不进行抑制，看看时频谱
figure(2);
[TF,F,T] = stft(p, fs, 'Window', window, 'OverlapLength', overlapLength, 'FFTLength', fftLength);  % （TF每一列就是每一个时间维度的频率）
subplot(7,2,[1,3]);
mesh(abs(TF)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' ); 
view(0,90); % 设置初始视角为俯视角
% mesh(angle(TF));
title('抑制前');


%% 一次谐波时频谱
TF1 = zeros(size(TF));
% TF1(2*windowLength+1:3*windowLength,:) = TF(1:windowLength,:);
TF1(2*windowLength+1:3*windowLength,:) = TF(3*windowLength+1:4*windowLength,:);  % 取出一次谐波
subplot(7,2,[2,4]);


mesh(abs(TF1)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' );
view(0,90);
title('一次谐波时频谱');


%% 二次谐波时频谱
TF2 = zeros(size(TF));
% TF1(2*windowLength+1:3*windowLength,:) = TF(1:windowLength,:);
TF2(2*windowLength+1:3*windowLength,:) = TF(4*windowLength+1:5*windowLength,:);  % 取出二次谐波
subplot(7,2,[5,7]);

mesh(abs(TF2)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' );
view(0,90);
title('二次谐波时频谱');






%% 
p = (istft(TF2, fs, 'Window', window, 'OverlapLength', overlapLength, 'FFTLength', 5*windowLength))';
p = real(p);
subplot(7,2,[6,8]);
plot(p);








%% 时频分析（抑制）
hfigure(2);
[TF,F,T] = stft(p, fs, 'Window', window, 'OverlapLength', overlapLength, 'FFTLength', fftLength);  % （TF每一列就是每一个时间维度的频率）
weight1 = abs(TF)./max(abs(TF));  % 要把这个矩阵当作权值使用，所以按列先归一化

% weight1(find(weight1 == 1)) = weight1(find(weight1 == 1)) - 1e-11;  % 算法1：是找到weight1中等于1的位置（理论上每列都有），并都减去num极小值，因为这个后边要用做分母，所以减去一个无穷小量，不要让他产生无穷大
% weight2 = (weight1.^64)./(1 - weight1.^64); % 算法1：构建一个函数，让原来频率分量大的更大，小的更小，windowLenght此处为抑制因子

weight2 = weight1;
weight2( find(weight2 < V) ) = 0; % 算法2：嘎嘎好用
weight2 = weight2 ./ max(weight2);  % 抑制后按列归一化
TF_curb = TF .* weight2;  % 获得抑制后的信号
subplot(7,2,[2,4,6]);
mesh(abs(TF)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' ); 
% view(0,90); % 设置初始视角为俯视角
% mesh(angle(TF));
title('抑制前');


%% 拉平时频谱去除包络
% TF_curb(find(TF_curb~=0)) = TF_curb(find(TF_curb~=0))./ abs(TF_curb(find(TF_curb~=0))) * 55; % 找到其不为0处的索引,把这些处的值设置除其各自的模值，这样这些复数的模长都为1，只不过相角不同
subplot(7,2,[10,12,14]);
mesh(abs(TF_curb)); set(gca,'YTickLabel',[]); ylabel('Frq.', 'FontSize',7,'FontWeight','bold' );
% mesh(angle(TF_curb));
% view(0,90);
title('抑制后');


%% 
p = (istft(TF_curb, fs, 'Window', window, 'OverlapLength', overlapLength, 'FFTLength', 5*windowLength))';
p_init2 = p;  % 用p_init2存储未归一化的信号
% p_init2 = p_init2(padding+1:end-padding);  % ✔去掉之前的padding,因为如果padding是122，那就是1-122，123-4122，4123-4244三段，所以实际上x是
% p_init2 = p_init2(padding+20:end-padding+19);  % ✔去掉之前的padding!!!!!!!!!!!!!!!!!!!!!!!!!!!这里是最大的误差
p_init2 = p_init2(padding+22:end-padding+21);  % ✔去掉之前的padding!!!!!!!!!!!!!!!!!!!!!!!!!!!这里是最大的误差
% p_init2 = p_init2(padding-25:end-padding-26);  % ✔去掉之前的padding!!!!!!!!!!!!!!!!!!!!!!!!!!!这里是最大的误差
p_init2 = sgolayfilt(p_init2,4,31);  % ✔ 超参数，平滑滤波器, 这里的参数也需要细调！！！c小的时候可以设置4，31   c大的时候可以设置，2，21 

figure(2);
subplot(7,2,3);
plot(p_init2);
title('未归一化');

%% 希尔伯特变换重构
direction = dir;

figure(2);
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


%% 傅里叶变换看频谱
figure(1);
subplot(9,1,2);
w = hamming(N);
f = fs / N * (0 : 1 : N-1);  % Fs/N就是这个频谱中的最小频率间隔！！！！！所以N越大，分辨率会越高
p_ = fft(p);

% 平移频域信号
fshift = (-N/2:N/2-1)*(fs/N);  % 平移后信号的频域范围
p_ = fftshift(p_);  % fftshift将零频分量移动到数组中心，重新排列
amp1 = abs(p_) * 2 / N ;
plot(amp1);
title("平移后频域信号（未更改频域范围）");
% plot(f(1:N/2), amp1(1:N/2));  % fft算法默认是双边谱,通常我们只取一半

%% 取出谐波成分，将平移后的频谱置零取出(丽萍学姐论文)
pp_ = takeHarmonicComponent(p_,6000,8000);  % 这个范围还是不能太小嗷，6500-7500误差就变大了
subplot(9,1,3);
amp2 = abs(pp_) * 2 / N ;
plot(fshift,amp2);
title("取出的一次谐波成分（更改了频域范围）");

subplot(9,1,4);
p1 = ifft(ifftshift(pp_));
p1 = imag(p1);
plot(p1);
title("一次谐波时域信号")


pp_ = takeHarmonicComponent(p_,8000,10000);
subplot(9,1,5);
p2 = ifft(ifftshift(pp_));
p2 = real(p2);
plot(p2);
title("二次谐波时域信号")

subplot(9,1,6);
phiF_wrapped = atan(p1./p2);
plot(phiF_wrapped);

%% arctan相位解包裹
for i = 2:N
    if (phiF_wrapped(i) - phiF_wrapped(i-1) > pi/2)
        phiF_wrapped(i:end) = phiF_wrapped(i:end) -  pi;
    elseif (phiF_wrapped(i) - phiF_wrapped(i-1) < -pi/2)  % 找峰值点，碰到峰值则加pi
        phiF_wrapped(i:end) = phiF_wrapped(i:end) +  pi;            
    end
end

phiF_reconstruct = phiF_wrapped;
subplot(9,1,7);
plot(phiF_reconstruct);
title("phiF-reconstruct")
%% 参数估算
[C_reconstruct, alpha_reconstruct] = SMI_API_ESTIMATE_C(phiF_reconstruct);  

%% 重构
subplot(9,1,8);
phi0_reconstrut = phiF_reconstruct + C_reconstruct*sin(phiF_reconstruct+atan(alpha));
Lt_reconstruct = phi0_reconstrut * lambda / (4 * pi);
plot(Lt,'k')
hold on;

Lt_reconstruct = Lt_reconstruct - mean(Lt_reconstruct); % 简写振动和余弦调制振动，重构后加上幅值A
% Lt_reconstruct = Lt_reconstruct + 1.5 * lambda;  % 重构后的随机振动信号要加上幅值1.5的波长，这是为啥我页不知道

plot(Lt,'k');
plot(Lt_reconstruct,'r')
title(['重构后的信号，C_reconstruct=', num2str(C_reconstruct)]);








%% 误差分析
subplot(9,1,9);
plot(Lt-Lt_reconstruct)
RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
title(['绝对误差，RMSE=', num2str(RMSE)])



%% 取出谐波信号函数
function pp_ = takeHarmonicComponent(p_,X,Y)  % p_为输入频谱。其中XY分别为输入频谱的起始点和终点。N为采样点数
     N = length(p_);
     pp_ = p_;
     pp_(1:X-1) = 0;
     pp_(Y+1:end) = 0;
     pp_(N/2-(Y-X)/2:N/2+(Y-X)/2) = pp_(X:Y);   
     pp_(X:Y) = 0;
end

%% 偶次幂函数！！！ 写了一个小时我真是醉了
function p_x = evenPower(x,p) % x为i
    if(x==0)
        p_x = p.^4-p.^2;  % i0
    else
        coff = 2^(x+1) -2;
        p_x = evenPower(x-1,p).^2 + 1 /2^coff *  evenPower(x-1,p);
    end 
end





%% subfunction （selfmixing-power）
function phiF = solve_phiF(C,phi0,alpha)  % 求解出每一个phi0对应的phiF
    if C<=1  % 每个phi0的解区间
        [phiF_min, phiF_max] = bounds1(C,phi0);
    else
        [phiF_min, phiF_max] = bounds2(C,phi0,alpha);
    end
    
    excessphaze_equation = @(phiF)phiF-phi0+C*sin(phiF+atan(alpha));
    
    if (excessphaze_equation (phiF_min)>0)  % 文章的解释为,phiF_min值可能刚好比零点大一点点点，这时候取phiF_min为近似零点
        phiF = phiF_min;
    elseif (excessphaze_equation (phiF_max)<0)
    	phiF = phiF_max;
    else
    	phiF = fzero(excessphaze_equation,[phiF_min,phiF_max]);
    end  
end

%---C < 1时的解区间函数------------------------------------------
function [phiF_min, phiF_max] = bounds1(C,phi0) 
    phiF_min = phi0-C;
    phiF_max = phi0+C;
end

%---C > 1时的解区间函数------------------------------------------
function [phiF_min, phiF_max] = bounds2(C,phi0,alpha)
persistent m; 
if isempty (m); m = 0; end
mlower = ceil ((phi0 + atan (alpha) + acos (1/C)- sqrt (C*C- 1))/(2*pi)- 1.5);
mupper = floor ((phi0 + atan (alpha)- acos (1/C) + sqrt (C*C- 1))/(2*pi)- 0.5);
if (m < mlower); m = mlower; end
if (m > mupper); m = mupper; end
phiF_min = (2*m+1)*pi + acos (1/C)- atan (alpha); 
phiF_max = (2*m+3)* pi- acos (1/C)- atan (alpha); 
end


%% subfunction2(reconstruction_T_N间转换)
% 该函数实现，将N转换为对应的t，也就是从T=1，N=10这些简单的时候推出来N(i)与t的关系
% 其中N为需要转换的点，N为采样点数，T为模拟时间
function  temp = N2T(temp1,N,T)
temp = zeros(1,length(temp1));    
    for i = 1:length(temp1) 
        temp(i) = T*(temp1(i)-1)/N;
    end
end

function temp = T2N(temp1,N,T)
temp = zeros(1,length(temp1));    
    for i = 1:length(temp1) 
        temp(i) = N*temp1(i)/T + 1;
    end
end