%% 代码源自AO 二次包络提取
%% 二次包络提取，使HT使用更广范围（实际上就是把它从大C变成小C）

function [p1,p2] = SMI_API_evenlopEXTRACT_HT_PRO(p,N)  % 其中p1为第一次提取后的结果，p2为第二次提取后的结果
    % 第一次包络提取
    [top1,loc1] = findpeaks(p);
    [top2,loc2] = findpeaks(-p);
    en_top = interp1(loc1,top1,1:N,'spline');  % 插值，曲线更平滑
    en_bottom = interp1(loc2,-top2,1:N,'spline');
    en_median = (en_top - (en_top - en_bottom)/2); 
    % 第一次包络提取结果（实际上是把翘起来的给往下拉）
    p1 = p - en_median; 

    % 第二次包络提取，实际上是往下拉的部分短了一点，通过除法变长
    [top3,loc3] = findpeaks(p1);
    en_top2 = interp1(loc3,top3,1:N,'spline');
    % 第二次包络提取结果
    p2 = p1 ./ en_top2;
    % p2 = 2 * (p2 - min(p2))./(max(p2)- min(p2)) - 1;
end