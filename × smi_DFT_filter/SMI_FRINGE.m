%% 条纹识别计数API，返回峰谷值、翻转点、条纹方向以及相关必要信息,目前来说鲁棒性很强，C在0.1~4的范围内都能很好用
function  [top_p, loc_p, top_v, loc_v, top_r, loc_r, direction] = SMI_FRINGE(p,N)
    %% 找到所有的峰谷值（包含跳变点）
    [top_p,loc_p] = findpeaks(p);
    % loc_p_convert = N2T(loc_p,N,T); 
    [top_v,loc_v] = findpeaks(-p);
    top_v = -top_v;
    % loc_v_convert = N2T(loc_v,N,T); 
       
    %% 包络确定方向！
    diffp = diff(p);  % diff是相邻两个数的差分，对第一个位置补0
    diffp = [0,diffp];
    diff_acosp = diff(acos(p));
    diff_acosp = [0,diff_acosp];

    [top_diffp_p,loc_diffp_p] = findpeaks(diffp);  % 拿到极值和索引值
    [top_diffp_v,loc_diffp_v] = findpeaks(-diffp);  

    en_top = interp1(loc_diffp_p,top_diffp_p,1:N,'spline');  % 三次样条插值，曲线更平滑
    en_bottom = interp1(loc_diffp_v,-top_diffp_v,1:N,'spline');

    en_median = (en_top - (en_top - en_bottom)/2); 
    dir = -sign(en_median);  % 让direction指明方向
    
    %% 利用方向信息，求出最合适的找翻转点的范围，大大增加了鲁棒性
    direction_seg1 = [];
    for i = 1:length(dir)-1
        if dir(i) ~=  dir(i+1)
            direction_seg1 = [direction_seg1, i];
        end
    end
    
    direction_seg2 = direction_seg1 + 1;
    direction_seg = [1, sort([direction_seg1,direction_seg2]), length(dir)];  % 这就是恒定正负一的区间
    top_diffp_seg = []; 
    loc_diffp_seg = [];  % 存储寻找翻转点的区间  
    
    for i = 1 : 2 : length(direction_seg)
        if dir(direction_seg(i):direction_seg(i+1)) < 0   % 在方向大于0的时候求diffp两端的峰值！！！！这样找翻转点就不用缩小范围了！！！！！！！！！！！！！！！！！！！！
            [tem1, tem2] = findpeaks(-diffp(direction_seg(i):direction_seg(i+1)));
            top_diffp_seg = [top_diffp_seg, -tem1(1), -tem1(end)];  
            loc_diffp_seg = [loc_diffp_seg, tem2(1)+ direction_seg(i) - 1 , tem2(end) + direction_seg(i) - 1];   % 这里返回的索引是个难点！由于我是分段进行峰值寻找，索引值应当加上初始值（假设在[x,y]中找极值，返回索引是2，第2是极值点,实际上对应的索引值为x+2-1）

        else  
            [tem3, tem4] = findpeaks(diffp(direction_seg(i):direction_seg(i+1)));
            top_diffp_seg = [top_diffp_seg, tem3(1), tem3(end)];
            loc_diffp_seg = [loc_diffp_seg, tem4(1) + direction_seg(i) - 1, tem4(end) + direction_seg(i) - 1];
        end          
    end
    
    top_diffp_seg([1,end]) = [];
    loc_diffp_seg([1,end]) = [];  % 最适合找翻转点的区间
    
    %% 取出翻转点，与上述方法相同！
    top_r = [];
    loc_r = [];
    for i = 1:2:length(loc_diffp_seg)
        [temp1,temp2] = findpeaks(p(loc_diffp_seg(i):loc_diffp_seg(i+1))); 
        [temp3,temp4] = findpeaks(-p(loc_diffp_seg(i):loc_diffp_seg(i+1)));
        top_r = [top_r, temp1, -temp3];
        loc_r = [loc_r, temp2 + loc_diffp_seg(i) - 1, temp4 + loc_diffp_seg(i) - 1 ];
    end
    
    %% 挖去峰值中的跳变点，返回，峰值、谷值、跳变点！
    for i = 1:length(loc_p)
        for j = 1:length(loc_r)
            if loc_p(i) == loc_r(j)
                loc_p(i) = nan;
                top_p(i) = nan;
            end
        end
    end

    % 挖去谷值中的跳变点
    for i = 1:length(loc_v)
        for j = 1:length(loc_r)
            if loc_v(i) == loc_r(j)
                loc_v(i) = nan;
                top_v(i) = nan;
            end
        end
    end

    % 删除数组中的nan
    loc_p(isnan(loc_p))=[];  
    loc_v(isnan(loc_v))=[]; 
    top_p(isnan(top_p))=[]; 
    top_v(isnan(top_v))=[]; 
    
    %% 根据求出的翻转点修正一下dir信息！
    direction = zeros(1,N);
    direction(1) = dir(1);
    k = 1;
    for i = 1:N
        direction(i) = k;
        if ismember(i, loc_r)
            k = k * -1;
        end
    end
    
    
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
