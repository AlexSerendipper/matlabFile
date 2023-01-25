% 基于峰谷点检测算法的相位解包裹重构唯位移，文中描述适用于C∈(0.1,4),即弱反馈到适度反馈区域
% diffp与diff_acosp的交点

%% global variance_1
T = 20e-3; %simulation time 
N = 4000; %numble of samples
t= 0:T/N:T-T/N; %根据采样点数，确定采样时间及间隔（T=1，N=10为例，以0.1为步进，保证t中个数与N相等）
lambda = 650e-9;
%---构筑外部简谐振动------------------
A = 2.25e-6;
L0 = 0.9;
f = 100;
Lt = L0 + A.*sin(2*pi*f*t); %L0为标称位置，外腔长度
%---global variance_2-----------------
C = 0.1;
alpha = 4.6;
beta = 1; %the amplitude of selfmixing signal
phi0 = 4*pi*Lt/lambda;


%% main function1 ,creat self mixing power
p = zeros(1,N);
for i = 1:N 
    p(i) = beta * cos(solve_phiF(C,phi0(i),alpha)); %遍历所有的phi0
end
subplot(4,1,1);
plot(t,Lt);
title("外部简谐振动");
subplot(4,1,2);
plot(p);
title("自混合信号");


%% main function(peak_valley_point)
subplot(4,1,3);
diffp = diff(p); %这两玩意求导后数量比N少一个,补了一个
diffp(N) = diffp(end);
diff_acosp = diff(acos(p));
diff_acosp(N) = diff_acosp(end);
plot(diffp);
hold on;
plot(diff_acosp);
hold on;
y = diffp - diff_acosp;
zero1 = intersection(y,N);
zero2 = zeros(1,N);
plot(zero1,zero2,'o');
hold on;
title("diffp与diffacosp的交点")
















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
function zero1 = intersection(y,N) %求数值曲线交点，y为两组数的差，N为数据的个数
zero1 = zeros(1,N);
for i = 1:N-1 % 数值曲线找交点
    if y(i) * y(i+1) == 0   
        if y(i) == 0
            zero1(i) = i;
        end
        if y(i+1) == 0
            zero1(i+1) = i+1;
        end
    elseif y(i) * y(i+1) < 0  %一定有交点，用一次插值
        k = abs(y(i))/(abs(y(i))+abs(y(i+1))); %交点在i与i+1之间的比例,这个....没看懂，应该就是插值吧
        zero1(i) = i + k;
    else            
    end
    if zero1(i) == 0   %除掉不是交点的部分
    zero1(i) = nan;
    end
end
end
