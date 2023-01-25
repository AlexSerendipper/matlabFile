% 求出所有phi0（Lt公式），便利所有的Phi0求出所有的phiF（超相位公式求零点），因为不同的c值他们的解区间不同，所以先确定不同C值下的解区间
% T=20e-3,A=2.25，L0=0.9，f=200,C=0.8，alpha=4，lamdba=650最接近王路学姐论文中的图
%% global variance_1
T = 10e-3;  % simulation time 
N = 2000;  % numble of samples
t= 0:T/N:T-T/N;  % 根据采样点数，确定采样时间及间隔（T=1，N=10为例，以0.1为步进，保证t中个数与N相等）
lambda = 650e-9;
%---构筑外部简谐振动------------------
A = 2.25e-6;
L0 = 0.9;
f = 200;
Lt = L0 + A.*sin(2*pi*f*t);  % L0为标称位置，外腔长度
%---global variance_2-----------------
C = 3;
alpha = 4.6;
beta = 1;  % the amplitude of selfmixing signal
phi0 = 4*pi*Lt/lambda;

%% main function,creat self mixing power
subplot(3,1,1);
p = zeros(1,N);
for i = 1:N 
    p(i) = beta * cos(solve_phiF(C,phi0(i),alpha)); %遍历所有的phi0
end
plot(t,Lt);
title("外部简谐振动");
subplot(3,1,2);
plot(t,p);
title("自混合信号");

%% for test
subplot(3,1,3);
diff_acosp = diff(acos(p));
diff_acosp = [diff_acosp(1),diff_acosp];
plot(t, diff_acosp);

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