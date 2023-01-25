%% cos函数相位翻倍函数，i为翻倍次数，p为输入的自混合信号,i=1时条纹个数翻一倍，i=2时条纹倍数再翻一倍...如 2 4 8..
%% 主要是Yn = 2 * Yn-1^2 -1的自迭代算法（迭代底为i1 = 2 * p.^2 - 1）
%% 源自熊彦彬学长论文,最后翻倍为2^i次方！！！！
function [p_doubled] = SMI_API_COSDOUBLE(i,p) 
    if(i==1)
        p_doubled = 2 * p.^2 -1;
    else
        p_doubled = 2 * SMI_API_COSDOUBLE(i-1,p).^2 - 1;
    end 
end