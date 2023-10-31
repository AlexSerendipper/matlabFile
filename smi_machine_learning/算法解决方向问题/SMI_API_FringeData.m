%% 实现运行一次，拿到一个C值下的所有条纹数据(以谷点为间隔)，及其对应方向
%% 循环遍历该函数，即可实现得到数据集

function  [fringeData] = SMI_API_FringeData(C,SNR,HP)
    %% 产生自混合信号
    fs = 200000;  % 采样率，即fs(s)采一个点。
    N = 8000;  
    fv = 50;  % 震动频
%     C = [2.8];  % C设置一个从a到b变化的值
    alpha = 5;
    [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_HARMONIC(fs, N, fv, C, alpha);  % 1 简谐振动的自混合信号
    [top_ov,loc_ov,top_p, loc_p, top_v, loc_v, top_r, loc_r, direction] = SMI_API_FRINGE(p,N);
    
    p = -1 + (p-min(p))/(max(p)-min(p))*2;
    
    p = awgn(p,SNR); 
    
    [top_v,loc_v] = findpeaks(-p);  % 'minpeakdistance'
    top_v = -top_v;


    %% first step ,calculate the MPD
    [c, lags] = xcorr(p,'unbiased');
    c = sgolayfilt(c,1,11);
    c = sgolayfilt(c,2,21);
    c = sgolayfilt(c,3,31);

    [top_c,loc_c] = findpeaks(c);
    MPD = mean(diff(loc_c));
   
    [top_v,loc_v] = findpeaks(-p, 'minpeakdistance',MPD/2);
    top_v = -top_v;

    %% second step,zero-crossing,所有大于0的都不要了
    loc_ov=[];
    top_ov=[];
    for i = 1:length(loc_v)
        if(top_v(i)<0)
            loc_ov = [loc_ov loc_v(i)];
            top_ov = [top_ov top_v(i)];
        end
    end

    %% third step,average && variance
    avg = mean(top_ov);  % 计算所有谷值点的平均值
    mea = mean(abs(top_ov - avg));  % 计算所有点 距离 谷点平均值 的平均值
    for i = 1:length(loc_ov)
        temp = top_ov(i) - avg;
        if (temp > mea)
            loc_ov(i) = nan;
            top_ov(i) = nan;
        end
    end
    loc_ov(isnan(loc_ov))=[];  
    top_ov(isnan(top_ov))=[]; 
    
    %% MDP judge,根据平均峰值判定，小于平均峰值的取更小点保留。。。
    % 需要多加一个判断是否两点的差值>mea/1.25，因为当C比较小的时候，是不需要平均峰值这一步的
    % 目前看来，HP在C<1时，使用1.4合适。在C为1~2时，使用1.25为佳。在C为2.5~2.8时，使用1.3较好。
    for i = 1:length(loc_ov)-1
        if(loc_ov(i+1)-loc_ov(i)<MPD/HP)
            if(top_ov(i+1)>top_ov(i))
                top_ov(i+1) = nan;
                loc_ov(i+1) = nan;
                i = i+2;
            else
                top_ov(i) = nan;
                loc_ov(i) = nan;
                i = i+2;
            end
        end
    end
    
    % 删除数组中的nan
    loc_ov(isnan(loc_ov))=[]; 
    top_ov(isnan(top_ov))=[]; 
    %% 拿信号
    fringeData = [];
    for i=2:length(loc_ov)
        % 得到谷值区间点数（条纹）
        N = loc_ov(i) - loc_ov(i-1);
        % 拿到谷值区间内所有的点
        dd = loc_ov(i-1)+1:loc_ov(i)-1;
        % 使用众数记录当前区间方向
        m = mode(direction(loc_ov(i-1):loc_ov(i)));
        % 判断谷值区间 内！，是否含有翻转点，如果有则 if=false
        if isempty(intersect(dd, loc_r))
            % 不包含，记录当前方向
            dir = m;
        else
            % 包含,则方向置零
            dir = 0;
        end
        loc = p(loc_ov(i-1):loc_ov(i));
        int_ = 1000;
        loc_ = [SMI_API_RESAMPLE(loc,int_) dir];  % 方向作为最后一个点
        fringeData = [fringeData;loc_];
    end
end
