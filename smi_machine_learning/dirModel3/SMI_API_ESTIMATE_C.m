%% 输入phiF_reconstruct，实现C和alpha的估计（alpha估计有些许偏差，但不影响重构结果）
function [C_reconstruct, alpha_reconstruct] = SMI_API_ESTIMATE_C(phiF_reconstruct)  
    C_reconstruct = [];
    alpha_reconstruct = [];
    xf = phiF_reconstruct;
    argmin = 1000000;  % 找最小值，所以初值设置大一点
    for i = 0.001: 0.001: 5  % 这里是C的范围  
        for j = 0.888:0.1:1.446  % 这里是α的范围，这里已经包含了arctan(1)~arctan(8)      
             x0 = xf + i .* sin(xf + j);  
             temp = sum( diff(x0) .^2 ); % 这里实际上是求差分的的最小值，即，[x0(i+1)-x0(i)]^2 ，理论上就是让phi0最连续的取值
             if temp < argmin
                 argmin = temp;
                 C_reconstruct = i;
                 alpha_reconstruct = tan(j); 
             end  
        end
    end
end
    