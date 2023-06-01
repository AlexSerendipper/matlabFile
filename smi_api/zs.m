%% 重构。 loc_v 就是 valley谷值点的位置，同理loc_p就是峰值点的位置。direction是方向信息，用包络那个方法可以得到

subplot(5, 1, 3);
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

plot(acosp_op1);
hold on;
title("翻转点×-1");


add_op = zeros(1,N);  % 累加阶梯

for i = 2:N  
    if ismember(i, loc_v) && (direction(i)==-1)
        add_op(i:end) = add_op(i:end) -  2 * pi;
    elseif ismember(i, loc_v) && (direction(i)==1)
        add_op(i:end) = add_op(i:end) +  2 * pi;
    end
end

phiF_reconstruct = acosp_op1 + add_op;

%% 参数估算（这里是我给你发过的 估算C的程序）
[C_reconstruct, alpha_reconstruct] = SMI_API_ESTIMATE_C(phiF_reconstruct);  

%% 重构
figure(1);
subplot(5,1,4);
phi0_reconstrut = phiF_reconstruct + C_reconstruct*sin(phiF_reconstruct+atan(alpha_reconstruct));  % 这里的alpha如果用估算的，就会引入蛮大的误差
Lt_reconstruct = phi0_reconstrut * lambda / (4 * pi);

plot(Lt,'k')
hold on;
plot(Lt_reconstruct,'r')
title("重构后的信号");

%% 误差分析
subplot(5,1,5);
plot(Lt-Lt_reconstruct)
RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
title(['绝对误差，RMSE=', num2str(RMSE)])