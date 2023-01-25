% Displacement Measurements Using a Self-Mixing Laser Diode Under Moderate Feedback,2006  以及 王路毕业论文 
% PUM目前来看，适用于 0.05 <= C <= 3.6,上限仅优化到3.6，如还需优化，需要先优化包络检测....
% 包络判断法在实际使用中最稳定范围是：0.05<= C <=3.6,详见envelop_discriminate_direction.m 
% ①还需要模拟随机运动，代码不会写②实际中，C与alpha的联合估计也需要做，但我没做
% 如能准确识别出跳变点，程序可以进一步优化
%% global variance_1
T = 20e-3; %simulation time 
N = 4000; %numble of samples
t= 0:T/N:T-T/N; %根据采样点数，确定采样时间及间隔（T=1，N=10为例，以0.1为步进，保证t中个数与N相等）
lambda = 650e-9;
%---构筑外部简谐振动------------------
A = 2.25e-6;
L0 = 0.9;
f = 200;
Lt = L0 + A.*sin(2*pi*f*t); %L0为标称位置，外腔长度
%---global variance_2-----------------
C = 1.3;
alpha = 4.6;
beta = 1; %the amplitude of selfmixing signal
phi0 = 4*pi*Lt/lambda;

%% main function1 ,creat self mixing power
p = zeros(1,N);
for i = 1:N 
    p(i) = beta * cos(solve_phiF(C,phi0(i),alpha)); %遍历所有的phi0
end
subplot(7,1,1);
plot(t,Lt);
title(['外部简谐振动,C= ',num2str(C)]);
subplot(7,1,2);
plot(t,p);
title("自混合信号");
hold on;

%% main function3(envelop)
% 局部峰谷值
[top_p,location_p] = findpeaks(p);
location_p_convert = N2T(location_p,N,T); 
scatter(location_p_convert,top_p);
[top_v,location_v] = findpeaks(-p);  
location_v_convert = N2T(location_v,N,T); 
scatter(location_v_convert,-top_v);
subplot(7,1,3);
% 求导后数量比N少一个,补一个,diff是相邻两个数的差分，补第一个,这个地方可能会有问题..
diffp = diff(p);  
diffp = [diffp(1),diffp];
diff_acosp = diff(acos(p));
diff_acosp = [diff_acosp(1),diff_acosp];
plot(t,diffp,'b');
hold on;
[top1,location1] = findpeaks(diffp);  % 拿到极大值和索引值
[top2,location2] = findpeaks(-diffp);  % 拿到极小值和索引值
location1 = N2T(location1,N,T);  % 转换为时域       
location2 = N2T(location2,N,T);
en_top = interp1(location1,top1,t,'spline');  % 三次样条插值，曲线更平滑
en_bottom = interp1(location2,-top2,t,'spline');
plot(t,en_top,"r");  % 画包络线
hold on;
plot(t,en_bottom,"g");
envelop_median = (en_top - (en_top - en_bottom)/2);
plot(t,envelop_median,"k");
hold on;
direction = -sign(envelop_median);
plot(t,direction);
title("diffp及其包络确定条纹方向");

%% main function4(correct_direction_fringe_counting) —— 难点1(首先处理)：direction跳变点与自混合信号跳变点不完全重合，需要拓宽来去除自混合跳变点，引出direction_append1    
%                                                          —— 难点2：direction_append1有N个点，要将把location_p以及location_v中其余点按顺序补NaN才能和direction_append1点乘                                                                                                                              
subplot(7,1,4);
direction_append1 = direction;
for i = N:-1:2  % 实现拓宽跳变点,从右往左遇到跳变点设为nan（难点1）   
    if direction_append1(i) ~= direction_append1(i-1)
        direction_append1(i) = nan;
    end
end
for i = 1:N-1  % 从左往左再次拓宽跳变点，增加鲁棒性，拓宽的范围与N的取值密切相关，如果拓宽过大可能会丢失条纹，目前设置来看拓宽13在N=2000有效，26在N=4000有效。
    if isnan(direction_append1(i)) == 1
        direction_append1(i-1:-1:i-13*N/2000) = nan;
    end
end
for i = N:-1:2  % 从右往左再次拓宽跳变点，增加鲁棒性
    if isnan(direction_append1(i)) == 1
        direction_append1(i+1:i+13*N/2000) = nan;
    end
end

% 实现补全nan（难点2）
location_p_insert = zeros(1,N)*nan;
top_p_insert = zeros(1,N)*nan;
for i = 1:length(location_p)  
    a = location_p(i);
    top_p_insert(a) = top_p(i);
    location_p_insert(a) = a;
end
location_p_convert = N2T(location_p_insert,N,T); 
top_p_convert = direction_append1 .* top_p_insert;
stem(location_p_convert,top_p_convert);
hold on;
title("条纹数量及方向");

%% main function 5(PUM_reconstruction)1       注意：此前仅在top_p中去除跳变点！，location_p_insert中仍包含着跳变点，引出direction_append2（跳变点附近为nan，其余均为1）,用于去location_p_insert上的跳变点
%                                                        location_v_insert也也需要去除跳变点  ！ 若将跳变点计算成条纹，会出现黄贞论文中，整体重构向下偏移的情况！！！！
subplot(7,1,5);
direction_append2 = zeros(1,N);
for i = 1:N
    direction_append2(i) = direction_append1(i);
    if direction_append1(i) == -1
       direction_append2(i) = 1;
    end
end           
location_v_insert = zeros(1,N) * nan;
for i = 1:length(location_v)  % 实现补全nan（难点2）
    a = location_v(i);
    location_v_insert(a) = a;
end
location_p_insert = location_p_insert .* direction_append2; 
location_v_insert = location_v_insert .* direction_append2;

% 设定初值，准备开始重构
acosp = acos(p);
j = 0;
if (p(1)-p(2))*(acosp(1)-acosp(2)) > 0
    i = 0;
else
    i = 1;
end

% PUM重构
phiF_reconstruct = zeros(1,N);
for k = 1:N
    if (~isnan(location_p_insert(k)) ==  1) || (~isnan(location_v_insert(k)) ==  1)  % 判断是否在P,V点
        i = i + 1;
    end
    if direction(k) > 0 && (~isnan(location_v_insert(k)) ==  1)  % 位移增大 且 到达V点
        j = j + 1;
    end
    if direction(k) < 0 && (~isnan(location_v_insert(k)) ==  1)
        j = j - 1;
    end
    phiF_reconstruct(k) = (-1)^i * acosp(N) + 2*j*pi;
end
phi0_reconstrut = phiF_reconstruct + C*sin(phiF_reconstruct+atan(alpha));  % 此时重构的结果呈阶梯型，尝试使用寻找突变点来画
Lt_reconstruct = phi0_reconstrut * lambda / (4 * pi) + L0;  % 我也不懂为啥这里重构出来的是不包括L0的
plot(t,Lt_reconstruct);
hold on;
plot(t,Lt);
title("PUM重构1");

%% main function 5(PUM_reconstruction)2: 尝试使用寻找突变点插值画出平滑图像
subplot(7,1,6)
phi0_reconstrut_smooth = [];
phi0_reconstrut_smooth_x = [];

for i = 1:N-1
    if phi0_reconstrut(i+1) > phi0_reconstrut(i)  % 上升、下降阶段都取突变值
       phi0_reconstrut_smooth = [phi0_reconstrut_smooth,phi0_reconstrut(i+1)];
       phi0_reconstrut_smooth_x = [phi0_reconstrut_smooth_x,i+1];
    end
    if phi0_reconstrut(i+1) < phi0_reconstrut(i)
       phi0_reconstrut_smooth = [phi0_reconstrut_smooth,phi0_reconstrut(i)];
       phi0_reconstrut_smooth_x = [phi0_reconstrut_smooth_x,i];
    end 
%     if direction(i) ~= direction(i+1)  % 保留跳变点处的一些值，否则插值在跳变点处有误差
%        phi0_reconstrut_smooth = [phi0_reconstrut_smooth,phi0_reconstrut(i-13*N/2000:i+13*N/2000)];
%        phi0_reconstrut_smooth_x = [phi0_reconstrut_smooth_x,i-13*N/2000:i+13*N/2000]; 
%     end
end
% scatter(phi0_reconstrut_smooth_x,phi0_reconstrut_smooth);  % 取突变值没有发现问题
% hold on;
phi0_reconstrut_smooth = interp1(phi0_reconstrut_smooth_x,phi0_reconstrut_smooth,1:N,'spline');
Lt_reconstruct = phi0_reconstrut_smooth * lambda / (4 * pi) + L0;  % 我也不懂为啥这里重构出来的是不包括L0的
plot(t,Lt_reconstruct);
hold on;
plot(t,Lt);
hold on;
title("PUM重构2");

%% (补充)判断谷点检测是否出现误差
% subplot(7,1,7);
% stem(location4_insert .* direction);
% title("确认谷值点位置");


























%% subfunction 1（selfmixing-power）
function phiF = solve_phiF(C,phi0,alpha) %求解出每一个phi0对应的phiF
    if C<=1 %每个phi0的解区间
        [phiF_min, phiF_max] = bounds1(C,phi0);
    else
        [phiF_min, phiF_max] = bounds2(C,phi0,alpha);
    end
    
    excessphaze_equation = @(phiF)phiF-phi0+C*sin(phiF+atan(alpha));
    if (excessphaze_equation (phiF_min)>0) %文章的解释为,phiF_min值可能刚好比零点大一点点点，这时候取phiF_min为近似零点
        phiF = phiF_min;
    elseif (excessphaze_equation (phiF_max)<0)
    	phiF = phiF_max;
    else
    	phiF = fzero(excessphaze_equation,[phiF_min,phiF_max]);
    end  
end
%---C<1时的解区间函数------------------------------------------
function [phiF_min, phiF_max] = bounds1(C,phi0) 
    phiF_min = phi0-C;
    phiF_max = phi0+C;
end
%---C>1时的解区间函数------------------------------------------
function [phiF_min, phiF_max] = bounds2(C,phi0,alpha)
persistent m; 
if isempty (m); m = 0; end
mlower = ceil ((phi0 + atan (alpha) + acos (1/C)- sqrt (C*C- 1))/(2*pi)- 1.5);
mupper = floor ((phi0 + atan (alpha)- acos (1/C) + sqrt (C*C- 1))/(2*pi)- 0.5);
if (m < mlower); m = mlower; end
if (m > mupper); m = mupper; end
phiF_min = (2*m+1)*pi + acos (1/C)- atan (alpha); 
phiF_max = (2*m+3)* pi- acos (1/C)- atan (alpha); 
end


%% subfunction2(reconstruction)

% 该函数实现，将N转换为对应的t，也就是从T=1，N=10这些简单的时候推出来N(i)与t的关系
% 其中N为需要转换的点，N为采样点数，T为模拟时间
function  temp = N2T(temp1,N,T)
temp = zeros(1,length(temp1));    
    for i = 1:length(temp1) 
        temp(i) = T*(temp1(i)-1)/N;
    end
end

