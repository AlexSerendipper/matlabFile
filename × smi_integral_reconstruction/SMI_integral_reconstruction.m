% Self-mixing interference displacement measurement under very weak feedback regime based on integral reconstruction method,2019,Optics Communications

%% 产生自混合信号
subplot(5, 1, 1);
fs = 200000;  % 采样周期
N = 4000;  % 采样点——采样点可能会报错，原因在于某些采样点让direction出错了，如4000。原因未知
fv = 100;  % 震动频率
C = 0.08;
alpha = 4.6;
[t, lambda, L0, Lt, phi0, p] = SMI_API(fs, N, fv, C, alpha);
plot(Lt);
title(['外部简谐振动,C=',num2str(C)]);
subplot(5, 1, 2);
plot(t, p);
title("自混合信号");

%% 得到重构所需的相关信息
[top_p, loc_p, top_v, loc_v, top_r, loc_r, direction] = SMI_FRINGE(p,N);
direction = -direction;  % 如果初始震动用的cos，那方向是反的，一定要注意
plot(p);
hold on;
scatter(loc_p,top_p);
scatter(loc_v,top_v);
scatter(loc_r,top_r);
plot(direction);
title("pks, vls, rev")

%% 积分重构
subplot(5, 1, 3);
diffp = [0, diff(p)];
diff_phi0 = sqrt( (diffp).^2 ./(1 - p.^2) ) .* direction;
phi0_reconstruct = cumsum(diff_phi0) ;
Lt_reconstruct = (lambda .* phi0_reconstruct) ./ (4*pi);
Lt_reconstruct = Lt_reconstruct - mean(Lt_reconstruct);
plot(Lt_reconstruct);
title("重构后的信号")
subplot(5, 1, 4);
plot(Lt - Lt_reconstruct);
title("误差")
