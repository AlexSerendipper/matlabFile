%% 需要提前知道方向信息，才能进行去直流
%% 本方法是采用的是: 除驼峰区外，其他所有的上下条纹，每一段都进行去直流的操作
%% 在尝试过
% （1）以翻转点为分割，所有段去直流，但是仅频域去直流有效✔
% （2）以翻转点为分割, 仅去上条纹【失效】
% （3）以翻转点为分割，去上条纹+驼峰区【失效】
% （4）以导数信号作为条纹判断的依据，每个条纹区间都减去其平均值【失效】
% （5）保留翻转点不动，上下条纹均去直流，时域频域去直流均有效✔
% （6）二次包络去直流，有效✔

function [p] = SMI_API_ELIMINATE_DC2(p_init,direction,method)
    % 需要传入初始方向
    dir = direction;
    diffp = diff(p_init);  % diff是相邻两个数的差分，对第一个位置补0
    diffp = [0,diffp];
    direction_seg1 = [];  % 方向发生变化的点(ˇ∀ˇ)
    for i = 1:length(dir)-1
        if dir(i) ~=  dir(i+1)
            direction_seg1 = [direction_seg1, i];
        end
    end

    direction_seg2 = direction_seg1 + 1;
    direction_seg = [1, sort([direction_seg1,direction_seg2]), length(dir)];  % 其中的1-2，3-4 为方向恒定的区间!!!!!!(ˇ∀ˇ)
    top_diffp_seg = [];
    loc_diffp_seg = [];  % 存储翻转点的区间

    for i = 1 : 2 : length(direction_seg)
        if dir(direction_seg(i):direction_seg(i+1)) > 0  % 在方向小于0的时候求diffp【两端】的谷值！！！！这样找翻转点就不用缩小范围了！！！！！！！！！！！！！！！！！！！！
            [tem1, tem2] = findpeaks(-diffp(direction_seg(i):direction_seg(i+1)));
            top_diffp_seg = [top_diffp_seg, -tem1(1), -tem1(end)];
            loc_diffp_seg = [loc_diffp_seg, tem2(1)+ direction_seg(i) - 1 , tem2(end) + direction_seg(i) - 1];   % 这里返回的索引是个难点！由于我是分段进行峰值寻找，索引值应当加上初始值（假设在[x,y]中找极值，返回索引是2，第2是极值点,实际上对应的索引值为x+2-1）
        else
            [tem3, tem4] = findpeaks(diffp(direction_seg(i):direction_seg(i+1)));  % 在方向大于0的时候求diffp【两端】的峰值！！！！这样找翻转点就不用缩小范围了！！！！！！！！！！！！！！！！！！！！
            top_diffp_seg = [top_diffp_seg, tem3(1), tem3(end)];
            loc_diffp_seg = [loc_diffp_seg, tem4(1) + direction_seg(i) - 1, tem4(end) + direction_seg(i) - 1];
        end
    end
    top_diffp_seg([1,end]) = [];
    loc_diffp_seg([1,end]) = [];  % 最适合找翻转点的区间
    
    direction_seg = [1, loc_diffp_seg, length(direction)]; 

    %% 每一段都去除直流分量！
%     direction_seg2 = [];  % 用来存储方向恒大于0的区间。
%     for i = 1:length(direction_seg)
%         if direction(direction_seg(i)) > 0
%             direction_seg2 = [direction_seg2, direction_seg(i)];    % 其中的1-2，3-4 为方向恒大于0的区间!!!!!!(ˇ∀ˇ)
%         end
%     end
    direction_seg2 = direction_seg;  %% 由于结果表明，每一段都去除直流分量效果更好，所以这里对direction_seg2的操作是无效的！ 
    
    %% 每一段在频域上去直流
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
