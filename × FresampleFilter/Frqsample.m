function [xa] = Frqsample(x, l, num)
% 频谱采样滤波
% 输入：信号， 频谱抑制系数（l>0）
% 关于频谱抑制系数：1. l 值越大，非主频抑制越强
%                  2. 单主频信号一般设定 l > 1; 多主频信号设定 0 < l < 1, 主频数越少 l 越接近 1
%                  3. num = 1 为频谱抑制抽样; num = 2 为频谱能量抽样
N = length(x); fft_p = fft(x, N); 

if num == 1
    % 构建代数函数实现频谱抑制抽样
    % 1.设定中间函数，即归一化
    Na = zeros(1,N); absf = abs(fft_p);
    for i = 1:length(absf)
        if absf(i) ~= max(absf), Na(i) = absf(i) ./ max(absf); 
        elseif absf(i) == max(absf), Na(i) = absf(i) ./ max(absf) - 1e-5; 
        end
    end
    
    % 2.中间函数构建代数函数
    Pv = ((Na).^l)./(1 - (Na).^l); 
    
    % 3.求代数函数的标准差，并提取幅度大于标准差的分量所在的频率值
    thegma = std(Pv); 
    omiga = find( Pv > thegma ); 
    if mod(numel(omiga),2)==1, omiga(1) = []; end  % 保证是偶数
    rebuild = zeros(length(omiga)/2, N);
    
    for i = 1:length(omiga)/2
        rebuild(i, [omiga(i) omiga(end+1-i)]) = fft_p([omiga(i) omiga(end+1-i)]);
    end
    
    % 4.频谱进行 Inverse Fourier Transform 重构信号
    rebuildtime = ifft(rebuild, N, 2);
    re_sig = zeros(1,N);
    for i = 1:length(rebuildtime(:,1))
        re_sig = re_sig + rebuildtime(i,:);
    end

elseif num == 2
    % 计算能量进行频谱能量抽样
    % 1.得到频谱能量分布
    powerf = 10*log10(abs(fft_p));
    powerf(1) = min(powerf);

    % 2.取半排序提取出前两大峰值索引
    peakst = find( powerf == max(powerf) );
    powerf(peakst) = min(powerf);
    peaknd = find( powerf == max(powerf) );
    
    % 4.将两索引对应fft系数提出
    rebuild = zeros(2,N);
    rebuild(1,peakst) = fft_p(peakst);
    rebuild(2,peaknd) = fft_p(peaknd);

    % Inverse FFT 获取时域信号
    rebuildtime = ifft(rebuild, N, 2);
    re_sig = zeros(1,N);
    for i = 1:length(rebuildtime(:,1))
        re_sig = re_sig + rebuildtime(i,:);
    end

end

xa = re_sig;

end