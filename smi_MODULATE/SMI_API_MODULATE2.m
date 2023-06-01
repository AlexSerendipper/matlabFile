%% 该函数接受自混合信号的phiF,h为调制深度，fm为调制频率，gamma为调制信号初始相位
%% 输出添加了调制的自混合信号,默认添加正弦信号
function [phi0,h] = SMI_API_MODULATE2(phi0,h,fm,gamma,t)
    phiM = h * sin(2*pi*fm*t + gamma);  % 调制信号
    phi0 = phi0 + 2*phiM;
    % 先后两次穿过置于外腔中的 EOM, 因此，实际的相位调制部分为2phiM
end

