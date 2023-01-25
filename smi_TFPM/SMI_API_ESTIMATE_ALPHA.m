%% 输入p，实现C和alpha的估计（基于驼峰区的值），这个方法仅适用于适度反馈，在C大于2时表现良好
%  在C值的估计上会引入0.0几的误差，而传统算法不会引入误差。在C大于2时alpha引入的误差小于0.1而传统算法误差大于0.5。所以该算法优势在于在C值较大时对alpha的估计。

function [C_reconstruct, alpha_reconstruct] = SMI_API_ESTIMATE_ALPHA(direction,loc_v,loc_p,top_v,top_p)  
    direction_seg1 = [];
    for i = 1:length(direction)-1
        if direction(i) ~=  direction(i+1)
            direction_seg1 = [direction_seg1, i];
        end
    end

    direction_seg2 = direction_seg1 + 1;
    direction_seg = [1, sort([direction_seg1,direction_seg2]), length(direction)];  % 其中的1-2，3-4 为方向恒定的区间!!!!!!(ˇ∀ˇ)

    H_top_seg_p=[];  % 存储高条纹的峰值
    H_loc_seg_p=[];
    H_top_seg_v=[];  % 存储高条纹的谷值
    H_loc_seg_v=[];
    h_top_seg_p=[];  % 存储低条纹的峰值
    h_loc_seg_p=[];
    h_top_seg_v=[];  % 存储低条纹的谷值
    h_loc_seg_v=[];

    for i = 1 : 2 : length(direction_seg)
        if direction(direction_seg(i):direction_seg(i+1)) > 0  % 找高条纹
            temp1 = loc_p((loc_p > direction_seg(i)) & (loc_p < direction_seg(i+1)));  % 拿出高条纹区域的所有峰值
            H_loc_seg_p = [H_loc_seg_p, temp1(end)];  % 拿出高条纹区域的最后一个峰值
            m = find(loc_p==temp1(end));  % 找到最后一个估值对应位置的索引
            H_top_seg_p = [H_top_seg_p, top_p(m)];  

            temp2 = loc_v((loc_v > direction_seg(i)) & (loc_v < direction_seg(i+1)));  % 拿出高条纹区域的所有谷值
            H_loc_seg_v = [H_loc_seg_v, temp2(end)];  % 拿出高条纹区域的最后一个谷值
            m = find(loc_v==temp2(end));  % 找到最后一个谷值对应位置的索引
            H_top_seg_v = [H_top_seg_v, top_v(m)];
        else 
            temp3 = loc_p((loc_p > direction_seg(i)) & (loc_p < direction_seg(i+1)));  % 拿出低条纹区域的所有峰值
            h_loc_seg_p = [h_loc_seg_p, temp3(1)];  % 拿出低条纹区域的第一个峰值
            m = find(loc_p==temp3(1));  % 找到第一个峰值对应位置的索引
            h_top_seg_p = [h_top_seg_p, top_p(m)];  

            temp4 = loc_v((loc_v > direction_seg(i)) & (loc_v < direction_seg(i+1)));  % 拿出低条纹区域的所有谷值
            h_loc_seg_v = [h_loc_seg_v, temp4(1)];  % 拿出低条纹区域的第一个谷值
            m = find(loc_v==temp4(1));  % 找到第一个谷值对应位置的索引
            h_top_seg_v = [h_top_seg_v, top_v(m)];
        end        
    end
    
    %% 估计alpha值
    %  H = 1 + abs(top_seg(1));  % this is high fringe height!!!!
    %  h = 1 + abs(top_seg(2));  % this is low fringe height!!!!
    H = H_top_seg_p(1) - H_top_seg_v(1);  % 高条纹
    h = h_top_seg_p(1) - h_top_seg_v(1);  % 低条纹
    A = (h - H)^2;
    B = (h + H - 2)^2 + (h - H)^2 - 4;
    D = (h + H - 2)^2;
    alpha_reconstruct = sqrt( (-B + sqrt(B^2 - 4*A*D))/(2*A) );

    %% 估计C值
    gamma = acos((h-H)/(2*cos(atan(alpha_reconstruct))));
    f = @(x) acos(-1/x) + x*sin(acos(-1/x)) + gamma + x*sin(gamma) - 2*pi;
    C_reconstruct = fsolve(f,1);
end
    

    