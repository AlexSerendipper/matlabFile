% 该程序实现了已知绝对距离为L的自混合波形。通常通过三角波来调制激光器偏置电流从而实现频率调制,影响phi0而后同样通过超相位方程，求出功率。
% 利用激光线性频率调制（实际中做不到线性），可用于确认绝对距离，亦或是静止远程目标的复折射率探测。
   
%% global variance_1
c = 3e8; %speed of light
T = 20e-3; %simulation time 
N = 1000; %numble of samples
t= 0:T/N:T-T/N; %根据采样点数，确定采样时间及间隔（T=1，N=10为例，以0.1为步进，保证t中个数与N相等）
lambda = 850e-9;
%---频率调制参数------------------
deltaF = -46e9; %Frequency modulation coefficient(Hz)
rho = 2.9; %Power modulation coefficient
tri = 1 + sign(mod(t/T,1)-0.5).*(1-2*mod(t/T,1)); %三角波，采样点为N
%---global variance_2-----------------
C = 2;
alpha = 4.6;
beta = 0.1; %the amplitude of selfmixing signal
L = 0.024;
phi0 = 4*pi*L*(1/lambda + deltaF * tri/c); %phi0受频率调制的公式
%% main function,creat self mixing power
p = zeros(1,N);
for i = 1:N 
    p(i) = beta * cos(solve_phiF(C,phi0(i),alpha)) + rho * tri(i); %遍历所有的phi0
end
subplot(2,1,1);
plot(t,p);
title("绝对距离下的自混合信号");
subplot(2,1,2);
plot(diff(p));
title("微分后的自混合信号");
% 若只有图，可以利用公式，L = N*c/(2deltaF)来估算绝对距离，N为微分后条纹数

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