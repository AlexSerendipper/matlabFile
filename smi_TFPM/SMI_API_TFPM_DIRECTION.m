%% 针对TFPM的方向识别，输入归一化的信号后，由于归一化的信号都是W型的SMI信号，利用平均峰值距离求出所有翻转点后，设置初始方向与包络求出的dir一致，碰到翻转点变号即可

function  direction = SMI_API_TFPM_DIRECTION(p, N, W)  % W用于增大平均峰值
    %% 找到所有的峰谷值（包含跳变点）
    [top_p,loc_p] = findpeaks(p);
    [top_v,loc_v] = findpeaks(-p);
    top_v = -top_v;

    %% 包络确定初步方向！
    diffp = diff(p);  % diff是相邻两个数的差分，对第一个位置补0
    diffp = [0,diffp];
    diff_acosp = diff(acos(p));
    diff_acosp = [0,diff_acosp];

    [top_diffp_p,loc_diffp_p] = findpeaks(diffp);  % 拿到极值和索引值
    [top_diffp_v,loc_diffp_v] = findpeaks(-diffp);
    % scatter(loc_diffp_p,top_diffp_p);
    % scatter(loc_diffp_v,-top_diffp_v);
    en_top = interp1(loc_diffp_p,top_diffp_p,1:N,'spline');  % 三次样条插值，曲线更平滑
    en_bottom = interp1(loc_diffp_v,-top_diffp_v,1:N,'spline');

    en_median = (en_top - (en_top - en_bottom)/2);
    dir = -sign(en_median);  % 让dir暂时指明方向


    %% 利用平均峰值，求出翻转点，因为翻转点都是M型的
    Avg_p = sum(diff(loc_p)) / length(loc_p) + W;  % 平均峰值距离，W可以增大平均峰谷值，因为发现条纹越密，平均峰值反而不够大
    
    loc_r = [];
    top_r = [];
    
    for i = length(loc_p)-1:-1:2
        % 求峰值的翻转点
        if abs(loc_p(i) - loc_p(i-1)) > Avg_p && abs(loc_p(i) - loc_p(i+1)) > Avg_p % 每个点都和左右边的点比，距离都大于平均峰值就是翻转点
            loc_r = [loc_r,loc_p(i)];
            top_r = [top_r, top_p(i)];     
        end
        % 求谷值的翻转点
        if abs(loc_v(i) - loc_v(i-1)) > Avg_p && abs(loc_v(i) - loc_v(i+1)) > Avg_p % 每个点都和左右边的点比，距离都大于平均峰值就是翻转点
            loc_r = [loc_r,loc_v(i)];
            top_r = [top_r, top_v(i)];     
        end  
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