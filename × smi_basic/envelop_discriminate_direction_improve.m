% 对于前述利用包络法实现条纹倾斜方向的判断以及实现正确的条纹计数出现的翻转点去除问题
% 2011_applid_optics_Improving the measurement performance for a self-mixing interferometry，
% 该文章提出了定位翻转点的解决方案,但这也有一些问题，①就是当C非常小的时候，diffp的图像呈包络状，很难搞~~~~~②取出的翻转点区间实际上包含了条纹，若适当缩小，将仅适用于C > 1.5的场合
% 在C<0.05时，sign函数出现问题（envelop_middle）(可能是由于diff函数我补的那个值？)，该问题仍然没有解决

%% global variance_1
T = 20e-3;  % simulation time 
N = 6000;  % numble of samples
t= 0:T/N:T-T/N;  % 根据采样点数，确定采样时间及间隔（T=1，N=10为例，以0.1为步进，保证t中个数与N相等）
lambda = 650e-9;
%---构筑外部简谐振动------------------
A = 2.25e-6;
L0 = 0.9;
f = 200;
Lt = L0 + A.*sin(2*pi*f*t);  % L0为标称位置，外腔长度
%---global variance_2-----------------
C = 5;
alpha = 4.6;
beta = 1;  % the amplitude of selfmixing signal
phi0 = 4*pi*Lt/lambda;

%% main function1 ,creat self mixing power
p = zeros(1,N);
for i = 1:N 
    p(i) = beta * cos(solve_phiF(C,phi0(i),alpha));  % 遍历所有的phi0
end
subplot(7,1,1);
plot(t,Lt);
title("外部简谐振动");
subplot(7,1,2);
plot(t,p);
title("自混合信号及峰谷值");
hold on;

%% main function2, find local peaks（top3,location3）
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

%% main function4(fringe_counting)                                                                                                                                 
subplot(7,1,4);
plot(t,diffp)
hold on;
scatter(location2,-top2);
title("diffp的局部极值（错误）")
average1 = mean(top1(:));  % 求出平均峰值
average2 = mean(-top2(:));
for i = 1:length(top1)  % 保留diffp中所有小于平均峰值的点
    if top1(i) < average1
        top1(i) = nan;
        location1(i) = nan;
    end
end
for i = 1:length(top2)
    if -top2(i) > average2
        top2(i) = nan;
        location2(i) = nan;
    end
end

%% main function 5(toss_reverse_point)        
% 我要把location1和location2以及他们对应的top1和top2组合成一个location4,top4,方便取出跳变区间!!!!!!!!!!!!!!
location1(isnan(location1))=[];  % 将矩阵location1中的nan全部删除
location2(isnan(location2))=[]; 
top1(isnan(top1))=[];
top2(isnan(top2))=[];
location4 = [location1,location2;top1,-top2].';  % 将横坐标合并并且排序
location4 = sortrows(location4);
top4 = location4(:,2).';  % 取出第二列的元素并转置
location4 = location4(:,1).';
subplot(7,1,5)
plot(t,diffp)
hold on;
scatter(location4,top4);
hold on;
location5 = [];  % 用来存储跳变区间，但该区间实际上包括了条纹，后续还需要处理
top5 = [];
for i = 1:length(location4)-1
    if top4(i) * top4(i+1) < 0
        top5 = [top5,top4(i),top4(i+1)];
        location5 = [location5,location4(i),location4(i+1)];
    end
end
scatter(location5,top5,'g');   
title("diffp的正确极值")

% 根据location5求reverse point,想法是将t转换为N，然后在p里求局部极值
subplot(7,1,6);
location5 = T2N(location5,N,T);  % 该区间实际上包括了条纹
location5 = round(location5);  % 不知道为啥有些是小数，我按照最接近取整
location6 = [];  % 用location6存储峰值跳变点
location7 = [];  % 用location7存储谷值跳变点
top6 = [];
top7 = [];
for i = 1:2:length(location5)
[temp1,temp2] = findpeaks(p(location5(i):location5(i+1)-N/2000*5));  % 适当缩小了找峰值的区间，因为跳变点区间包含了条纹
[temp3,temp4] = findpeaks(-p(location5(i):location5(i+1)-N/2000*5));
location6 = [location6,temp2+location5(i)-1];  % 这里返回的索引是个难点！由于我是分段进行峰值寻找，索引值应当加上初始值（假设在[x,y]中找极值，返回索引是2，第2是极值点,实际上对应p中的索引值为x+2-1）
top6 = [top6,temp1];
location7 = [location7,temp4+location5(i)-1];
top7 = [top7,temp3];
end
plot(t,p);
hold on;
location6_convert = N2T(location6,N,T); 
scatter(location6_convert,top6);
location7_convert = N2T(location7,N,T); 
scatter(location7_convert,-top7);
title("确定的峰值（谷值）点")


% 挖去峰值中的跳变点
for i = 1:length(location_p)
    for j = 1:length(location6)
        if location_p(i) == location6(j)
            location_p(i) = nan;
            top_p(i) = nan;
        end
    end
end
% 挖去谷值中的跳变点
for i = 1:length(location_v)
    for j = 1:length(location7)
        if location_v(i) == location7(j)
            location_v(i) = nan;
            top_v(i) = nan;
        end
    end
end        
subplot(7,1,7);
plot(t,p)
hold on;
location_p_convert = N2T(location_p,N,T); 
scatter(location_p_convert,top_p);
location_v_convert = N2T(location_v,N,T); 
scatter(location_v_convert,-top_v);
title("自混合信号峰谷值（不包含跳变点）")

































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

