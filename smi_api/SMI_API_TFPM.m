%% 时频抑制算法，最终返回抑制后逆变换回时域的信号,一般V为0.65即可
function  [T,F,TF,TF_curb,p] = SMI_API_TFPM(p,N,fs,windowLength,V)
    %% 全局变量
    overlapLength = floor(windowLength * 0.9);  % OverlapLength后为指定的重叠长度
    window = hamming(windowLength, "periodic");  % 使用汉明窗作为滑动的窗口
    fftLength = 5 * windowLength;

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

    %% 时频分析（抑制）
    [TF,F,T] = stft(p, fs, 'Window', window, 'OverlapLength', overlapLength, 'FFTLength', fftLength);  % （wf每一列就是每一个时间维度的频率）
    TF_curb = TF_inhibit1(TF,V);
%     TF_curb = TF_inhibit2(TF);
%     TF_curb = TF_inhibit3(TF,64);
    %% 归一化
    p = (istft(TF_curb, fs, 'Window', window, 'OverlapLength', overlapLength, 'FFTLength', 5*windowLength))';
    p = p(padding+1:end-padding);  % ✔去掉之前的padding
end


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
    weight1(find(weight1 == 1)) = weight1(find(weight1 == 1)) - 1e-11;  % 算法1：是找到weight1中等于1的位置（理论上每列都有，即主频位置），并都减去num极小值，因为这个后边要用做分母，所以减去一个无穷小量，不要让他产生无穷大
    weight2 = (weight1.^64)./(1 - weight1.^64);  % 算法1：构建一个函数，让原来频率分量大的更大，小的更小，64此处为抑制因子
    weight2 = weight2 ./ max(weight2);  % 抑制后按列归一化，用作权值矩阵
    TF_curb = TF .* weight2;  % 获得抑制后的信号
end

% 对TF_curb，抑制后的信号进行处理，因此此时主频大，其他成分极小，所以想实现放大C值的效果后用导数实现方向判断
function TF_curb = TF_inhibit3(TF,k)
    weight1 = abs(TF)./max(abs(TF));
    weight2 = (weight1.^k)./(weight1.^k+(1-weight1).^k);  % 构建的这个函数让靠近1的数无限趋近1，靠近0的数无限趋近0
    weight3 = weight1 .* weight2;
    TF_curb = TF .* weight3;  % 获得抑制后的信号
end

