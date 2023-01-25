% 基于envelop_discriminate_direction_improve的PUM重构，在1.5<C<4的范围内实用性较好,当C > 4会丢条纹啊！！！！！！！！！

%% global variance_1
T = 20e-3;  % simulation time 
N = 4000;  % numble of samples
t= 0:T/N:T-T/N;  % 根据采样点数，确定采样时间及间隔（T=1，N=10为例，以0.1为步进，保证t中个数与N相等）
lambda = 650e-9;
%---构筑外部简谐振动------------------
A = 2.25e-6;
L0 = 0.9;
f = 200;
Lt = L0 + A.*cos(2*pi*f*t);  % L0为标称位置，外腔长度
%---global variance_2-----------------
C = 1;
alpha = 4.6;
beta = 1;  % the amplitude of selfmixing signal
phi0 = 4*pi*Lt/lambda;

%% main function1 ,creat self mixing power & find local peaks
subplot(7,1,1);
p = zeros(1,N);
for i = 1:N 
    p(i) = beta * cos(solve_phiF(C,phi0(i),alpha));  % 遍历所有的phi0
end
plot(t,Lt);
title(['外部简谐振动,C= ', num2str(C)]);

% find local peaks
subplot(7,1,2);
plot(t,p);
title("自混合信号及峰谷值");
hold on;
[top_p,location_p] = findpeaks(p);
location_p_convert = N2T(location_p,N,T); 
scatter(location_p_convert,top_p);
[top_v,location_v] = findpeaks(-p);
location_v_convert = N2T(location_v,N,T); 
scatter(location_v_convert,-top_v);

%% main function3(envelop)
subplot(7,1,3);
diffp = diff(p);  % 这两玩意求导后数量比N少一个,补一个,diff是相邻两个数的差分，补第一个
diffp = [diffp(1),diffp];
diff_acosp = diff(acos(p));
diff_acosp = [diff_acosp(1),diff_acosp];
plot(t,diffp,'b');
hold on;
[top1,location1] = findpeaks(diffp);  % 拿到极大值和索引值
[top2,location2] = findpeaks(-diffp);  % 拿到极小值和索引值
location1 = N2T(location1,N,T); % 转换为时域  
location2 = N2T(location2,N,T);
%scatter(location1, top1)
en_top = interp1(location1,top1,t,'spline');  % 三次样条插值，曲线更平滑
en_bottom = interp1(location2,-top2,t,'spline');
plot(t,en_top,"r");  % 画包络线
hold on;
plot(t,en_bottom,"g");
envelop_median = (en_top - (en_top - en_bottom)/2);
plot(t,envelop_median,"k");
hold on;
direction = -sign(envelop_median);  % 让direction指明方向
plot(t,direction);
title("diffp及其包络确定条纹方向");
                                                                                                                            
%% 原先利用平均峰值去除翻转点，效果不佳，现利用方向信息，求出最合适的找翻转点的范围，大大增加了鲁棒性
subplot(7,1,4)
direction_seg1 = [];
for i = 1:length(direction)-1
    if direction(i) ~=  direction(i+1)
        direction_seg1 = [direction_seg1, i];
    end
end

direction_seg2 = direction_seg1 + 1;
direction_seg = [1, sort([direction_seg1,direction_seg2]), length(direction)];  % 这就是恒定正负一的区间

top_diffp_seg = [];
loc_diffp_seg = [];
for i = 1 : 2 : length(direction_seg)
    if direction(direction_seg(i):direction_seg(i+1)) < 0   % 在大于0的时候求两端的峰值！！！！这样就不用缩小范围了！！！！！！！！！！！！！！！！！！！！！！
        [tem1, tem2] = findpeaks(-diffp(direction_seg(i):direction_seg(i+1)));
        top_diffp_seg = [top_diffp_seg, -tem1(1), -tem1(end)];  
        loc_diffp_seg = [loc_diffp_seg, tem2(1)+ direction_seg(i) - 1 , tem2(end) + direction_seg(i) - 1];   % 这里返回的索引是个难点！由于我是分段进行峰值寻找，索引值应当加上初始值（假设在[x,y]中找极值，返回索引是2，第2是极值点,实际上对应的索引值为x+2-1）
        
    else  
        [tem3, tem4] = findpeaks(diffp(direction_seg(i):direction_seg(i+1)));
        top_diffp_seg = [top_diffp_seg, tem3(1), tem3(end)];
        loc_diffp_seg = [loc_diffp_seg, tem4(1) + direction_seg(i) - 1, tem4(end) + direction_seg(i) - 1];
    end          
end
top_diffp_seg([1,end]) = [];
loc_diffp_seg([1,end]) = [];
plot(t, diffp);
hold on;
loc_diffp_seg_convert = N2T(loc_diffp_seg,N,T);
scatter(loc_diffp_seg_convert, top_diffp_seg);
title("在该区间内找跳变点");
hold on;

%% 取出翻转点，与上述方法相同！
subplot(7, 1, 5)
top_r = [];
loc_r = [];
for i = 1:2:length(loc_diffp_seg)
    [temp1,temp2] = findpeaks(p(loc_diffp_seg(i):loc_diffp_seg(i+1))); 
    [temp3,temp4] = findpeaks(-p(loc_diffp_seg(i):loc_diffp_seg(i+1)));
    top_r = [top_r, temp1, -temp3];
    loc_r = [loc_r, temp2 + loc_diffp_seg(i) - 1, temp4 + loc_diffp_seg(i) - 1 ]
end
loc_r_convert = N2T(loc_r,N,T);
plot(t, p);
hold on;
scatter(loc_r_convert, top_r);
hold on;
title("翻转点")

%% 挖去峰值中的跳变点，返回，峰值、谷值、跳变点！
subplot(7, 1, 6);

for i = 1:length(location_p_convert)
    for j = 1:length(loc_r_convert)
        if location_p_convert(i) == loc_r_convert(j)
            location_p_convert(i) = nan;
            top_p(i) = nan;
        end
    end
end

% 挖去谷值中的跳变点
for i = 1:length(location_v_convert)
    for j = 1:length(loc_r_convert)
        if location_v_convert(i) == loc_r_convert(j)
            location_v_convert(i) = nan;
            top_v(i) = nan;
        end
    end
end

% 删除数组中的nan
location_p_convert(isnan(location_p_convert))=[];  
location_v_convert(isnan(location_v_convert))=[]; 
top_p(isnan(top_p))=[]; 
top_v(isnan(top_v))=[]; 

plot(t,p);
hold on;
scatter(location_v_convert,-top_v);
hold on;
scatter(location_p_convert,top_p);
hold on;
title("自混合信号峰谷值（不包含跳变点）")


%%  main functino 6_PUM重构
subplot(7,1,7);
% 在峰谷点中补全nan
location_p_nan = zeros(1,N) * nan;
location_v_nan = zeros(1,N) * nan;

for i = 1:length(location_p)  % 实现补全nan（难点）
    a = location_p(i);
    location_p_nan(a) = a;
end

for i = 1:length(location_v)  % 实现补全nan（难点）
    a = location_v(i);
    location_v_nan(a) = a;
end

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
    if (~isnan(location_p_nan(k)) ==  1) || (~isnan(location_v_nan(k)) ==  1)  % 判断是否在P,V点
        i = i + 1;
    end
    if direction(k) > 0 && (~isnan(location_v_nan(k)) ==  1)  % 位移增大 且 到达V点
        j = j + 1;
    end
    if direction(k) < 0 && (~isnan(location_v_nan(k)) ==  1)
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


%% subfunction2(reconstruction_T_N间转换)
% 该函数实现，将N转换为对应的t，也就是从T=1，N=10这些简单的时候推出来N(i)与t的关系
% 其中N为需要转换的点，N为采样点数，T为模拟时间
function  temp = N2T(temp1,N,T)
temp = zeros(1,length(temp1));    
    for i = 1:length(temp1) 
        temp(i) = T*(temp1(i)-1)/N;
    end
end

function temp = T2N(temp1,N,T)
temp = zeros(1,length(temp1));    
    for i = 1:length(temp1) 
        temp(i) = N*temp1(i)/T + 1;
    end
end

