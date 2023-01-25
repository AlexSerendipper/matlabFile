%% 偶次幂函数！！！ 主要是Yn = Yn-1^2 + coff * Yn-1的自迭代算法（迭代底为i0 = p.^4-p.^2）
%% 输入i为对应论文的翻倍情况，i=-1,条纹个数×2，i=0,条纹个数×4，i=1,条纹个数×8以此类推 
%% times，即翻倍的倍数为2^i+2,很关键
%% 当C>0.5，偶次幂之后自混合信号的翻转点可能会由峰值变谷值，所以偶次幂局限性很强
function  p_evenpowered = SMI_API_EVENPOWER(i,p)
    if(i==-2)
         p_evenpowered = p;
    elseif(i==-1)
         p_evenpowered = p.^2;
    elseif(i==0)
         p_evenpowered = p.^4-p.^2;  % i0
    else
        coff = 2^(i+1) -2;
        p_evenpowered = SMI_API_EVENPOWER(i-1,p).^2 + 1 /2^coff *  SMI_API_EVENPOWER(i-1,p);
    end 
end

