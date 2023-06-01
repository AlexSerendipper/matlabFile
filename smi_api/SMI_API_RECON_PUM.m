%% 注意：传入的需要为处理后的纯净的峰谷值
function [phiF_reconstruct,Lt_reconstruct] = SMI_API_RECON_PUM(p,loc_p,loc_v,direction,N,lambda,c,alpha)
    %% 解包裹重构
    acosp = acos(p);
    plot(acosp);
    title("acosp");
    hold on;

    % 设定初值,固定
    if sign(acosp(2) - acosp(1)) == direction(1)
        init = 1;
    else
        init = -1;
    end

    % 遇到峰谷值乘-1
    mul_op = init * ones(1,N);
    for i = 1:N
        if ismember(i, [loc_v,loc_p]) == 1  % 当遇到翻转点，变一个方向（折叠）
            mul_op(i:end) = -mul_op(i:end);
        end
    end

    acosp_op1 =  acosp .* mul_op;

    add_op = zeros(1,N);  % 累加阶梯

    for i = 2:N  
        if ismember(i, loc_v) && (direction(i)==-1)
            add_op(i:end) = add_op(i:end) -  2 * pi;
        elseif ismember(i, loc_v) && (direction(i)==1)
            add_op(i:end) = add_op(i:end) +  2 * pi;
        end
    end

    % add_op = add_op - (max(add_op) + min(add_op))/2;
    phiF_reconstruct = acosp_op1 + add_op;


    %% 参数估算
    if(c~=0 && alpha~=0)
        C_reconstruct = c;
        alpha_reconstruct = alpha;
    else
        [C_reconstruct, alpha_reconstruct] = SMI_API_ESTIMATE_C(phiF_reconstruct); 
    end

    %% 重构
    phi0_reconstrut = phiF_reconstruct + C_reconstruct*sin(phiF_reconstruct+atan(alpha_reconstruct));  % 这里的alpha如果用估算的，就会引入蛮大的误差
    Lt_reconstruct = phi0_reconstrut * lambda / (4 * pi);
    % Lt_reconstruct = Lt_reconstruct - mean(Lt_reconstruct); % 简写振动和余弦调制振动，重构后加上幅值A
    % Lt_reconstruct = Lt_reconstruct + 1.5 * lambda;  % 重构后的随机振动信号要加上幅值1.5的波长，这是为啥我页不知道

end


