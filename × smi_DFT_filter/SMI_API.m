%% 自混合信号源API,产生自混合信号并返回必要参数
function [t, lambda, L0, Lt, phi0, p] = SMI_API(fs, N, fv, C, alpha)   % fs为采样率(s)，每秒钟采样多少个点。采样率要比采样点数大的多才能不失真！
    % T = 1/fs;  % 采样周期（s）,几秒钟采一个点
    t = (0:N-1)/fs;  % 采样时间，设N=10, fs=200，即采样了0.05s，t为[0...0.045]
    lambda = 650e-9;  % 波长
    A = 2 * lambda;  % 幅值(该设置在频谱上不出错)
    L0 = 20 * lambda;  % 外腔距离 (该设置在频谱上不出错)
    Lt = A.* sin(2*pi*fv*t);  % L0为标称位置，外腔长度
    beta = 1;  % the amplitude of selfmixing signal
    phi0 = 4*pi*(L0+Lt)/lambda;

    p = zeros(1,N);
    for i = 1:N 
        p(i) = beta * cos(solve_phiF(C, phi0(i), alpha));  % 遍历所有的phi0
    end
end





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