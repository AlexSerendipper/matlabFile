%% 本方法是采用的是: 除驼峰区外，其他所有的上下条纹，每一段都进行去直流的操作
%% 通过平均峰值距离确认驼峰区，而后将剩余部分全部去直流。因此不需要外部方向信息

function [p] = SMI_API_ELIMINATE_DC3(p_init,method)
    p = sgolayfilt(p_init,1,11);
    p = sgolayfilt(p,2,21);
    p = sgolayfilt(p,3,31);
    N = length(p);
    %% 找到所有的峰谷值（包含跳变点）
    [top_p,loc_p] = findpeaks(p);
    [top_v,loc_v] = findpeaks(-p);
    top_v = -top_v;
    
    dir = zeros(1,N);
    MPD2 = mean(diff(loc_v))/0.85;
    for i = 1:length(loc_v)-1
         if(loc_v(i+1)-loc_v(i)>MPD2)  % 因为驼峰区存在很多异常情况，所以将驼峰区直接置零
                                      % 同时这个地方存在一定的鲁棒性，即通常情况下，驼峰区通常是对称大于平均峰值距离的~通过方向矫正，都能得到正确的方向
            dir(loc_v(i):loc_v(i+1)) = 0;
         end
    end
    
    direction_seg2 = [];
    for i=1:length(dir)
        if dir(i)==0
            direction_seg2 = [direction_seg2,i];
        end
    end
    direction_seg2 = [1,direction_seg2,N];
    
    
    %% 每一段都去除直流分量！
%     direction_seg2 = [];  % 用来存储方向恒大于0的区间。
%     for i = 1:length(direction_seg)
%         if direction(direction_seg(i)) > 0
%             direction_seg2 = [direction_seg2, direction_seg(i)];    % 其中的1-2，3-4 为方向恒大于0的区间!!!!!!(ˇ∀ˇ)
%         end
%     end

%     direction_seg2 = direction_seg;  %% 由于结果表明，每一段都去除直流分量效果更好，所以这里对direction_seg2的操作是无效的！ 
    
    %% 每一段在时域/频域上去直流
    for i = 1 : 2 : length(direction_seg2)
        p_ = p_init(direction_seg2(i):direction_seg2(i+1));
        if method=="time"
            g_ = eliminate_dc_time(p_);
        else
            g_ = eliminate_dc_freq(p_);
        end
        % g_ = -1 + (g_ - min(g_))/(max(g_)- min(g_)) * 2;  % 不能进行归一化，归一化过程中会引入很大的直流分量
        p_init(direction_seg2(i):direction_seg2(i+1)) = g_;
    end
    p = p_init;
end

%% subfunction （eliminate_dc,输入一段信号p，及其长度N ，输出去除直流分量后的信号p_）
function g_ = eliminate_dc_freq(p_)
    p_= fft(p_);
    p_(1) = 0;
    g_ = ifft(p_);
end

function g_ = eliminate_dc_time(p_)
    ave = mean(p_);
    g_ = p_-ave;
end
