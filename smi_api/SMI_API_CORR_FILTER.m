%% Adaptive Filter beasd on Sgolayfilter.
% 方法1，目前针对带噪信号的主要滤波手段为调用sgolayfilt函数，如p = sgolayfilt(p,1,11); p = sgolayfilt(p,2,21); p = sgolayfilt(p,3,31)。直到到达理想的效果
% 方法1 的优势是三次平滑化处理会让信号略微失真，在处理实验信号上反而有意想不到的效果
% 方法2：基于自相关的滤波方法相比与上述方法的优点在于参数可调性，通过调整smoothingfactor,threshold可以滤至我们想要的效果...参考参数（0.1，0.98）（0.2，0.9）

function x = SMI_API_CORR_FILTER(p,smoothingfactor,threshold)
    % 数据拓延，左右各拓百分之10
    extend_length = round(0.1*length(p));
    p_ = wextend("1D", "sym", p, extend_length);

    % 计算自混合信号的自相关系数及其滞后，plot(lags,c)后可以看见当lags等于0时，自相关性最强
    [c, lags] = xcorr(p_, p_);
    corr1 = c(lags == 0);

    % 初值
    flag = 1;
    x0 = p_;

    while flag
        x0 = smoothdata(x0,'sgolay','SmoothingFactor',smoothingfactor);  % 窗口因子大小，接近0的值会产生较小的移动窗口长度，从而导致较少的平滑处理。接近 1 的值会产生较大的移动窗口长度，从而促成较多的平滑处理。
        [c, lags] = xcorr(x0, p_);
        corr2 = c(lags == 0);
        % 用滤波后的自相关系数最大值与滤波前的自相关系数最大值做比值，小于阈值则停止循环
        J = corr2/corr1;
        if J <= threshold
            flag = 0;
        end  
    end

    % 去除延拓的数据
    x = x0(extend_length+1:end-extend_length);
end

