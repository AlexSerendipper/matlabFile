%% 源自熊彦彬学长论文，取出一次谐波和二次谐波，作为sin和cos输入
%% 主要为一次谐波迭代为：P1(n) = 2*P1(n-1).*P2(n-1)
%       二次谐波迭代为：P2(n) = P2(n-1).^2 - P1(n-1)^2
%% p1p2为输入的一次二次谐波（已去除贝塞尔分量），i为翻倍次数：0不翻倍，1翻一倍，times = 2^i
%% 源自熊彦彬学长论文,最后翻倍为2^i次方！！！！
function [phiF_wrapped] = SMI_API_TANDOUBLE(i,p1,p2)
        phiF_wrapped = atan(solve1(i,p1,p2)./solve2(i,p1,p2));     
end


function [Pn] = solve1(i,p1,p2)
    if(i==0)
        Pn = p1;
    else
        Pn = 2*solve1(i-1,p1,p2).*solve2(i-1,p1,p2);
    end
end

function [Pnn] = solve2(i,p1,p2)
    if(i==0)
        Pnn = p2;
    else
        Pnn = solve2(i-1,p1,p2).^2 - solve1(i-1,p1,p2).^2;
    end
    
end