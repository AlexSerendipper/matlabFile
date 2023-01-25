%Self-mixing signal generation for harmonic motion.

function [p, fs, N, lambda0, d, alpha] = harmonic_motion(C, K, f, num,j)

fs = 8000*f;                                                               % Sample rate0
N = 3*8000;                                                                % Number of samples
beta = 1;                                                                  % Laser power modulation coefficient
lambda0 = 532e-9;                                                          % Laser wavelength(m)
L0 = 20*lambda0;                                                           % Target nominal distance(m)
A = K*lambda0;                                                             % Motion amplitude(μm)
t = (0:N-1)/fs;                                                            % Sample times

if num == 1
    d = A.*cos(2*pi*f*t+pi/4);                                             % Displacement samples
elseif num==2
    d = A.*sin(2*pi*f*t+pi/4);  
elseif num==3
    d = A.*cos(2*pi*f*t);
elseif num==4
    d = A.*sin(2*pi*f*t);
elseif num==5
    d = A.*sin(2*pi*f*t + j*pi/100);
elseif num==6                                                              
    d = A.*cos(2*pi*f*t) + 0.5*A.*cos(2*pi*0.5*f*t) + 0.2*A.*cos(2*pi*1.2*f*t);
elseif num==7
    A1 = 8/3 .* lambda0; A2 = 4/3 .* lambda0;
    fv1 = 200; fv2 = 100;
    d = (A1 + A2 * sin(2*pi*fv2*t)) .* sin(2*pi*fv1*t);
elseif num==8
    x=[0 30 55 82 160 210 250 300 364 414 480 500]*16; 
    y=[0 0.5 1.2 1.6 0.7 0.88 0.65 0.40 0.20 0.08 0.01 0]*12*10^(-6); 
    xx = 1:8000;  yy=interp1(x,y,xx,'spline');
    d = [yy yy yy];     
end

phi0 = 4 * pi * (L0 + d)/lambda0;                                          % Round-trip phase samples
p = zeros(1, length(d));                                                   % Self-mixing signal samples
alpha = 4;                                                                


if j == 1
for i = 1 : length(d)                                                      % Generate the synthetic self-mixing signal
    p(i) = beta * cos(selfmixingpower(C(i), phi0(i), alpha));
end 
elseif j == 2
for i = 1 : length(d)                                                      % Generate the synthetic self-mixing signal
    p(i) = beta * sin(selfmixingpower(C(i), phi0(i), alpha));
end 
end
end

% Functions for solving self-mixing equations

function power = selfmixingpower(C, phi0, alpha)                           % Power level at a sample in time

if (C <= 1.0)
    [phimin, phimax] = boundsweak(C, phi0);
else
    [phimin, phimax] = boundsstrong(C, phi0, alpha);
end

excessphase = @(x)x - phi0 + C.*sin(x + atan(alpha));

%If the value at the left bound positive, then it will be very close to the solution.
%If the value at the upper bound is negative, it will be very close to the solution.

if (excessphase(phimin) > 0)
%     excessphase(phimin);
    phi = phimax;
elseif (excessphase(phimax) < 0)
%     excessphase(phimax);
    phi = phimax;
else
    phi = fzero(excessphase, [phimin, phimax]);
end

power = (phi);

end

function [phimin, phimax] = boundsweak(C, phi0)                            %Find search region when C <= 1.0

phimin = phi0 - C;
phimax = phi0 + C;

end

function [phimin, phimax] = boundsstrong(C, phi0, alpha)                   %Find search region when C >= 1.0

persistent m;   %Solution region number.

if isempty (m); m = 0; end

%Calculate upper & lower values of m where solutions exist then ensure m is between them

mlower = ceil( (phi0 + atan(alpha) + acos(1/C) - sqrt(C*C - 1))/(2*pi) - 1.5 );
mupper = floor( (phi0 + atan(alpha) - acos(1/C) + sqrt(C*C - 1))/(2*pi) - 0.5 );

if (m < mlower); m = mlower; end
if (m > mupper); m = mupper; end

phimin = (2*m + 1)*pi + acos(1/C) - atan(alpha); %Trough
phimax = (2*m + 3)*pi - acos(1/C) - atan(alpha); %Peak

end
