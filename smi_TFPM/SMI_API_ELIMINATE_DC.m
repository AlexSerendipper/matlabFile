%% 去除简谐振动中direction>0部分条纹的直流分量（当C越来越大，该部分引入越来越大的直流分量）
%% 所以实际上必须要知道方向信息，才能进行去直流（但结果表明，所有段去直流效果更佳）

function [p] = SMI_API_ELIMINATE_DC(p_init,direction)
    % 重置方向发生变化的区间direction_seg1
    direction_seg1 = [];
    for i = 1:length(direction)-1
        if direction(i) ~=  direction(i+1)
            direction_seg1 = [direction_seg1, i];
        end
    end
    direction_seg2 = direction_seg1 + 1;
    direction_seg = [1, sort([direction_seg1,direction_seg2]), length(direction)];  % 其中的1-2，3-4 为方向恒定的区间!!!!!!(ˇ∀ˇ)

    direction_seg2 = [];  % 用来存储方向恒大于0的区间。
    for i = 1:length(direction_seg)
        if direction(direction_seg(i)) > 0
            direction_seg2 = [direction_seg2, direction_seg(i)];    % 其中的1-2，3-4 为方向恒大于0的区间!!!!!!(ˇ∀ˇ)
        end
    end
    
    %% 每一段都去除直流分量！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！(ˇ∀ˇ)
    direction_seg2 = direction_seg; 
    
    %%
    % 对direction大于0的段进行去直流
    for i = 1 : 2 : length(direction_seg2)
        p_ = p_init(direction_seg2(i):direction_seg2(i+1));
        g_ = eliminate_dc(p_);
        % g_ = -1 + (g_ - min(g_))/(max(g_)- min(g_)) * 2;  % 不能进行归一化，归一化过程中会引入很大的直流分量
        p_init(direction_seg2(i):direction_seg2(i+1)) = g_;
    end
    p = p_init;
end





%% subfunction （eliminate_dc,输入一段信号p，及其长度N ，输出去除直流分量后的信号p_）
function g_ = eliminate_dc(p_)
    p_= fft(p_);
    p_(1) = 0;
    g_ = ifft(p_);
end