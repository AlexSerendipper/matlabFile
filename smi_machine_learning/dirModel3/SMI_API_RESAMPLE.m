%% 输入一段信号，重采样为N个点~
function loc_ = SMI_API_RESAMPLE(loc,int_)   % 输入信号为loc,重采样为int_个点
    fs = int_;  % 因为原本采样率为1，放大后采样率就是int_
    N = length(loc);
    % 将原信号重采样为原来的int_倍
    loc_ = interp1(1:int_:int_* N,loc,1:int_* N,'spline');
    % 再将信号降采样为原来的N倍，即 N * 1000 / N， 故最后信号长度为int_
    loc_ = loc_(1:N:length(loc_));
end










