%% 加载一些保存好的随机振动，并返回必要参数
%% 使用随机振动时，方向×负

function [t, lambda, L0, Lt, phi0, p, c] = MOVE_API_PPG_LOAD(fs, N, C, alpha)   % fs为采样率(s)，每秒钟采样多少个点。采样率要比采样点数大的多才能不失真！N需要为2次幂,否则频谱混叠（✔）N = KM（M为运动周期）
    % T = 1/fs;  % 采样周期（s）,几秒钟采一个点
    t = (0:N-1)/fs;  % 采样时间，设N=10, fs=200，即采样了0.05s，t为[0...0.045]
    
    lambda = 650e-9;  % 波长
    load('D:\matlab save\smi_脉搏波仿真\1.mat');
    L0 = 5 * lambda;  % 外腔距离（✔） 
    Lt = signals.PPG;
    
    
    padding = N - length(Lt);  % 把这个信号的长度补到N的长度
    Lt = (10^-5 * Lt)/3;  % 脉搏波信号大约是30个微米这样（3个微米！！！） 
    Lt = [Lt, zeros(1,padding)];
   
    
    beta = 1;  % the amplitude of selfmixing signal
    phi0 = 4*pi*(L0+Lt)/lambda;

    p = zeros(1,N);
    
    if length(C) == 1
        c = ones(1,N) * C;
        for i = 1:N 
            p(i) = beta * cos(solve_phiF(C, phi0(i), alpha));  % 遍历所有的phi0
        end
    elseif length(C) == 2
        % C的变化是一个正弦曲线
        C_lower = C(1);
        C_upper = C(2);
        % 这个乘和加保证了c的上下限
        x = linspace(0, 3*pi, N);
        c = (C_upper-C_lower)/2 * cos(x) + (C_upper - (C_upper-C_lower)/2);
        for i = 1:N 
            C = c(i);
            p(i) = beta * cos(solve_phiF(C, phi0(i), alpha));  % 遍历所有的phi0
        end
    elseif length(C) == 3
        load('D:\matlab save\self-mixing\smi_api\DATA_LLT_ALEATORY_4000(1).mat');  
        C_lower = C(1);
        C_upper = C(2);
        c = C_lower + (LLt - min(LLt))/(max(LLt)-min(LLt))*(C_upper-C_lower);
        for i = 1:N 
            C = c(i);
            p(i) = beta * cos(solve_phiF(C, phi0(i), alpha));  % 遍历所有的phi0
        end
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