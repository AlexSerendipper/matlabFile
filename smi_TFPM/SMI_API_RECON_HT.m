%% 希尔伯特变换重构
function  [phiF_wrapped,phiF_reconstruct,Lt_reconstruct] = SMI_API_RECON_HT(p,direction,N,lambda,c,alpha)
    %% hilbert transform
    inverse_hb = (imag(hilbert(p))) .* direction;  % 得到sin 
    phiF_wrapped = atan(inverse_hb./p);
    
    %% arctan相位解包裹
    for i = 2:N
        if (phiF_wrapped(i) - phiF_wrapped(i-1) > pi/2)
            phiF_wrapped(i:end) = phiF_wrapped(i:end) -  pi;
        elseif (phiF_wrapped(i) - phiF_wrapped(i-1) < -pi/2)  % 找峰值点，碰到峰值则加pi
            phiF_wrapped(i:end) = phiF_wrapped(i:end) +  pi;            
        end
    end

    phiF_reconstruct = phiF_wrapped;

    %% 参数估算
    if(c~=0 && alpha~=0)
        C_reconstruct = c;
        alpha_reconstruct = alpha;
    else
        [C_reconstruct, alpha_reconstruct] = SMI_API_ESTIMATE_C(phiF_reconstruct); 
    end
     
    %% 重构
    phi0_reconstrut = phiF_reconstruct + C_reconstruct*sin(phiF_reconstruct+atan(alpha_reconstruct));
    Lt_reconstruct = phi0_reconstrut * lambda / (4 * pi);

    Lt_reconstruct = Lt_reconstruct - mean(Lt_reconstruct); % 简写振动和余弦调制振动，重构后加上幅值A
    % Lt_reconstruct = Lt_reconstruct + 1.5 * lambda;  % 重构后的随机振动信号要加上幅值1.5的波长，这是为啥我页不知道
end    
