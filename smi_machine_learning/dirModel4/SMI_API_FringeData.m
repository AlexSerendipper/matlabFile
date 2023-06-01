%% 实现运行一次，拿到一个C值下的所有条纹数据(以谷点为间隔)，及其对应方向
%% 循环遍历该函数，即可实现得到数据集

function  [fringeData] = SMI_API_FringeData(C,SNR)
    %% 产生自混合信号
    fs = 200000;  % 采样率，即fs(s)采一个点。
    N = 8000;  
    fv = 50;  % 震动频
    alpha = 5;
    [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_HARMONIC(fs, N, fv, C, alpha);  % 1 简谐振动的自混合信号

    %% 得到重构所需的相关信息
    [top_ov,loc_ov,top_p, loc_p, top_v, loc_v, top_r, loc_r, direction] = SMI_API_FRINGE(p,N);
     p = awgn(p,SNR);  % 10db，加高斯白噪声
    
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
        loc_ = [SMI_API_RESAMPLE(loc,int_) dir];
        fringeData = [fringeData;loc_];
    end
end
