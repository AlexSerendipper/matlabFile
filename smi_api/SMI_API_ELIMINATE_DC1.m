%% 需要提前知道方向信息，才能进行去直流
%% 本方法是采用的是以翻转点为分割，每一段在频域上进行去直流的操作，是论文中使用的最早的方法
%% 在尝试过
% （1）以翻转点为分割，所有段去直流，但是仅频域去直流有效✔
% （2）以翻转点为分割, 仅去上条纹【失效】
% （3）以翻转点为分割，去上条纹+驼峰区【失效】
% （4）以导数信号作为条纹判断的依据，每个条纹区间都减去其平均值【失效】
% （5）保留翻转点不动，上下条纹均去直流，时域频域去直流均有效✔
% （6）二次包络去直流，有效✔

function [p] = SMI_API_ELIMINATE_DC1(p_init,direction,method)
    % 重置方向发生变化的区间direction_seg1
    direction_seg1 = [];
    for i = 1:length(direction)-1
        if direction(i) ~=  direction(i+1)
            direction_seg1 = [direction_seg1, i];
        end
    end
    direction_seg2 = direction_seg1 + 1;
    direction_seg = [1, sort([direction_seg1,direction_seg2]), length(direction)];  % 其中的1-2，3-4 为方向恒定的区间!!!!!!(ˇ∀ˇ)

    %% 每一段都去除直流分量！
%     direction_seg2 = [];  % 用来存储方向恒大于0的区间。
%     for i = 1:length(direction_seg)
%         if direction(direction_seg(i)) > 0
%             direction_seg2 = [direction_seg2, direction_seg(i)];    % 其中的1-2，3-4 为方向恒大于0的区间!!!!!!(ˇ∀ˇ)
%         end
%     end
    direction_seg2 = direction_seg;  %% 由于结果表明，每一段都去除直流分量效果更好，所以这里对direction_seg2的操作是无效的！ 
    
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
