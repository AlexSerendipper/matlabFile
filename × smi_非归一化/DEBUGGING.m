clc;
clear all;
close all;
fs = 200000;  
N = 4000;  
fv = 200;  
alpha = 4;

%% 产生自混合信号
figure(1);
subplot(7,1,1);
t = (0:N-1)/fs; 
lambda = 650e-9; 
A = 2 * lambda;  
L0 = 20 * lambda; 
Lt = A.* sin(2*pi*fv*t); 
beta = 1;  % the amplitude of selfmixing signal
phi0 = 4*pi*(L0+Lt)/lambda;

p = zeros(1,N);
C = [0.5, 0.5];
C_lower = C(1);
C_upper = C(2);

x = linspace(0, 3*pi, N);
c = (C_upper-C_lower)/2 * cos(x) + (C_upper - (C_upper-C_lower)/2);
% plot(x,c);

for i = 1:N 
    C = c(i);
    p(i) = beta * cos(solve_phiF(C, phi0(i), alpha));  % 遍历所有的phi0
end
plot(p);

title("自混合信号")

%% 找到所有的峰谷值
[top_p,loc_p] = findpeaks(p);
[top_v,loc_v] = findpeaks(-p);
top_v = -top_v;

%% 包络确定方向！
subplot(7,1,2);
diffp = diff(p); 
diffp = [0,diffp];
diff_acosp = diff(acos(p));
diff_acosp = [0,diff_acosp];
plot(diffp);
hold on;

[top_diffp_p,loc_diffp_p] = findpeaks(diffp);  
[top_diffp_v,loc_diffp_v] = findpeaks(-diffp);
en_top = interp1(loc_diffp_p,top_diffp_p,1:N,'spline');
en_bottom = interp1(loc_diffp_v,-top_diffp_v,1:N,'spline');

en_median = (en_top - (en_top - en_bottom)/2);
dir = -sign(en_median); 

plot(en_top);
plot(en_bottom);
plot(dir);
title("diffp及其包络确定方向")

%% 利用方向信息，求出最合适的找翻转点的范围
subplot(7,1,3);
direction_seg1 = [];  % 方向发生变化的点(ˇ∀ˇ)
for i = 1:length(dir)-1
    if dir(i) ~=  dir(i+1)
        direction_seg1 = [direction_seg1, i];
    end
end

direction_seg2 = direction_seg1 + 1;
direction_seg = [1, sort([direction_seg1,direction_seg2]), length(dir)]; 
top_diffp_seg = [];
loc_diffp_seg = []; 

for i = 1 : 2 : length(direction_seg)
    if dir(direction_seg(i):direction_seg(i+1)) < 0  
        [tem1, tem2] = findpeaks(-diffp(direction_seg(i):direction_seg(i+1)));
        top_diffp_seg = [top_diffp_seg, -tem1(1), -tem1(end)];
        loc_diffp_seg = [loc_diffp_seg, tem2(1)+ direction_seg(i) - 1 , tem2(end) + direction_seg(i) - 1];   
        
    else
        [tem3, tem4] = findpeaks(diffp(direction_seg(i):direction_seg(i+1))); 
        top_diffp_seg = [top_diffp_seg, tem3(1), tem3(end)];
        loc_diffp_seg = [loc_diffp_seg, tem4(1) + direction_seg(i) - 1, tem4(end) + direction_seg(i) - 1];
    end
end

top_diffp_seg([1,end]) = [];
loc_diffp_seg([1,end]) = [];  
plot(p);
hold on;
scatter(direction_seg,0,"r");
title("diffp极值区间（红点）")

subplot(7,1,4);
plot(diffp);
hold on;
scatter(loc_diffp_seg,0,"g");
title("(绿点)内不再包含条纹~~~")

%% 取出翻转点，与上述方法相同！
subplot(7,1,5);
top_r = [];
loc_r = [];
for i = 1:2:length(loc_diffp_seg)
    [temp1,temp2] = findpeaks(p(loc_diffp_seg(i):loc_diffp_seg(i+1)));
    [temp3,temp4] = findpeaks(-p(loc_diffp_seg(i):loc_diffp_seg(i+1)));
    top_r = [top_r, temp1, -temp3];
    loc_r = [loc_r, temp2 + loc_diffp_seg(i) - 1, temp4 + loc_diffp_seg(i) - 1 ];
end
plot(p);
hold on;
scatter(loc_r,top_r);
title("在上述区间内求出翻转点")

%% 挖去峰值中的跳变点
for i = 1:length(loc_p)
    for j = 1:length(loc_r)
        if loc_p(i) == loc_r(j)
            loc_p(i) = nan;
            top_p(i) = nan;
        end
    end
end

% 挖去谷值中的跳变点
for i = 1:length(loc_v)
    for j = 1:length(loc_r)
        if loc_v(i) == loc_r(j)
            loc_v(i) = nan;
            top_v(i) = nan;
        end
    end
end

% 删除数组中的nan
loc_p(isnan(loc_p))=[];
loc_v(isnan(loc_v))=[];
top_p(isnan(top_p))=[];
top_v(isnan(top_v))=[];


%% 根据求出的翻转点修正一下dir信息！
subplot(7,1,6);
direction = zeros(1,N);
direction(1) = dir(1);
k = 1;
for i = 1:N
    direction(i) = k;
    if ismember(i, loc_r)
        k = k * -1;
    end
end

plot(p);
hold on;
scatter(loc_p,top_p);
scatter(loc_v,top_v);
scatter(loc_r,top_r);
plot(direction);
title("得到纯净的峰谷值和完全正确的方向")

subplot(7,1,7);
plot(Lt);
hold on

%% 非归一化函数实现
figure(2);
fringe = [[loc_p, loc_r, loc_v]', [top_p, top_r, top_v]'];
fringe = sortrows(fringe,1);
loc_fringe = fringe(:,1)'; 
top_fringe = fringe(:,2)';
loc_fringe = [1, loc_fringe, N]; 
top_fringe = [p(1), top_fringe, p(N)]; 
distance = 0; 
Lt_reconstruct = []; 

for i=1:length(loc_fringe)-1
    if i-1==0 
       fringe = p(loc_fringe(i):loc_fringe(i+1));
       Lt_reconstruct = [Lt_reconstruct, p(loc_fringe(i):loc_fringe(i+1)) + distance];

    
    elseif (mode(direction(loc_fringe(i)+1:loc_fringe(i+1)) > 0)) && (mode(diff(p(loc_fringe(i)+1:loc_fringe(i+1))) < 0))  % 他这个得返回单个标量值我才能判断,所以加了个众数
            last_fringe = fringe;
            fringe = -p(loc_fringe(i)+1:loc_fringe(i+1));
            distance = distance + last_fringe(end) - fringe(1);
            Lt_reconstruct = [Lt_reconstruct, fringe + distance];
            
    
    elseif (mode(direction(loc_fringe(i):loc_fringe(i+1)) > 0)) && (mode(diff(p(loc_fringe(i):loc_fringe(i+1))) > 0))
            last_fringe = fringe;
            fringe = p(loc_fringe(i)+1:loc_fringe(i+1));           
            distance = distance + last_fringe(end) - fringe(1);
            Lt_reconstruct = [Lt_reconstruct, fringe + distance];
            
    
    elseif (mode(direction(loc_fringe(i)+1:loc_fringe(i+1)) < 0)) && (mode(diff(p(loc_fringe(i)+1:loc_fringe(i+1))) > 0))
            last_fringe = fringe;
            fringe = -p(loc_fringe(i)+1:loc_fringe(i+1));           
            distance = distance + last_fringe(end) - fringe(1);
            Lt_reconstruct = [Lt_reconstruct, fringe + distance];
        
    
    else (mode(direction(loc_fringe(i)+1:loc_fringe(i+1)) < 0)) && (mode(diff(p(loc_fringe(i)+1:loc_fringe(i+1))) < 0));
            last_fringe = fringe;
            fringe = p(loc_fringe(i)+1:loc_fringe(i+1));           
            distance = distance + last_fringe(end) - fringe(1);
            Lt_reconstruct = [Lt_reconstruct, fringe + distance];
    end 
end

plot(1:N,Lt_reconstruct);












%% subfunction （selfmixing-power）
function phiF = solve_phiF(C,phi0,alpha)  % 求解出每一个phi0对应的phiF
    if C<=1  % 每个phi0的解区间
        [phiF_min, phiF_max] = bounds1(C,phi0);
    else
        [phiF_min, phiF_max] = bounds2(C,phi0,alpha);
    end
    
    excessphaze_equation = @(phiF)phiF-phi0+C*sin(phiF+atan(alpha));
    
    if (excessphaze_equation (phiF_min)>0)  % 文章的解释为,phiF_min值可能刚好比零点大一点点点，这时候取phiF_min为近似零点
        phiF = phiF_min;
    elseif (excessphaze_equation (phiF_max)<0)
    	phiF = phiF_max;
    else
    	phiF = fzero(excessphaze_equation,[phiF_min,phiF_max]);
    end  
end

%---C < 1时的解区间函数------------------------------------------
function [phiF_min, phiF_max] = bounds1(C,phi0) 
    phiF_min = phi0-C;
    phiF_max = phi0+C;
end

%---C > 1时的解区间函数------------------------------------------
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


