clear; clc; close all;

f = 10; K = 0.25; num = 2; C = 0.1+0*cos(2*pi*20*(1:24000)/80000); 
[p, fs, N, lambda0, d, alpha] = harmonic_motion(C, K, f, num,1); t = (1:N)/fs;
direction = sign(diff(d)); direction = [direction(1) direction];
subplot 511; plot(t,p); grid on; axis tight; 
%% PUM
acosp = acos(p); q = 1;
[~, ploc] = findpeaks(p, "MinPeakHeight", 0.9);
[~, vloc] = findpeaks(-p, "MinPeakHeight",0.999);

reverse = q*ones(1, N);
plus = zeros(1,N);
for i = 1:N
    if ismember(i, [ploc vloc])
        reverse(i:end) = -reverse(i:end);
    end
end
axf = acosp .* reverse;
for i = 2:N
    dif = axf(i) - axf(i-1);
    if dif > pi/2
        plus(i:end) = plus(i:end) - 2*pi;
    elseif dif < -pi/2
        plus(i:end) = plus(i:end) + 2*pi;
    end
end

xf = axf + plus; 
% xf = xf + (max(xf)-min(xf))/2;
x0 = xf + C.*sin(xf + atan(alpha));
L1 = lambda0 * x0 / (4*pi); L1 = L1 - (L1(1) - d(1));
subplot 523; plot(t,L1,t,d); grid on; axis tight; title("PUM Error = " + num2str(RMSE(L1,d)));
subplot 524; plot(t, L1-d); grid on; axis tight; 
%% Tan xf
[pf, ~, ~, ~, ~, ~] = harmonic_motion(C, K, f, num,2);
atanx = atan(pf./p); xf = atanx;
for i = 2:N
    dif = atanx(i)-atanx(i-1);
    if dif > pi/2
        xf(i:end) = xf(i:end) - pi;
    elseif dif < -pi/2
        xf(i:end) = xf(i:end) + pi;
    end
end
xf = xf - (max(xf)+min(xf))/2;
x0 = xf + C.*sin(xf + atan(alpha));
L2 = lambda0 * x0 / (4*pi); L2 = L2 - (L2(1) - d(1));
subplot 525; plot(t,L2,t,d); grid on; axis tight; title("Arctan Error = " + num2str(RMSE(L2,d)));
subplot 526; plot(t,L2-d); grid on; axis tight;
%% Even-power Unwrapping
p0 = p.^4 - p.^2;
p1 = p0.^2 + p0/(2.^2);
p2 = p1.^2 + p1/(2.^6);

times = 2.^4;
pn = (p2 + 1/(2.^15)) .* (2.^15); 
acospn = acos(pn); q = -1;

[pks, ploc] = findpeaks(pn, "MinPeakHeight", 0.9);
[vly, vloc] = findpeaks(-pn, "MinPeakHeight",0.9);

reverse = q*ones(1, N);
plus = zeros(1,N);
for i = 1:N
    if ismember(i, [ploc vloc])
        reverse(i:end) = -reverse(i:end);
    end
end
axfn = acospn .* reverse;
for i = 2:N
    dif = axfn(i) - axfn(i-1);
    if dif > pi/2
        plus(i:end) = plus(i:end) - 2*pi;
    elseif dif < -pi/2
        plus(i:end) = plus(i:end) + 2*pi;
    end
end
xf = axfn + plus; xf = xf + (max(xf)-min(xf))/2;
x0 = (xf + C.*sin(xf + atan(alpha)))/times; 
L3 = lambda0 * x0 / (4*pi); L3 = L3 - (L3(1) - d(1));
subplot 527; plot(t,L3,t,d); grid on; axis tight; title("Even-power PUM Error = " + num2str(RMSE(L3,d)));
subplot 528; plot(t,L3-d);grid on; axis tight;
%% Even-power Interpolation
p0 = p.^4 - p.^2;
p1 = p0.^2 + p0/(2.^2);
p2 = p1.^2 + p1/(2.^6);
p3 = p2.^2 + p2/(2.^14);

pn = (p3 + 1/(2.^31)) .* (2.^31);  
times = 2.^5;

dpx = (lambda0/2)/times; % 一根条纹对应的波长
[vtop, vloc] = findpeaks(-pn, "MinPeakHeight",0.9);
vtop = - vtop;
vtop(ismember(vloc, find(abs(diff(direction))==2))) = [];
vloc(ismember(vloc, find(abs(diff(direction))==2))) = [];
% 到这就拿到所有的谷值点，不包含峰谷点
value = zeros(1,length(vloc));

% 数条纹啊
for i = 2:length(value)
    if direction(vloc(i))==direction(vloc(i-1))
        value(i) = value(i-1) + dpx*direction(vloc(i));
    elseif direction(vloc(i))~=direction(vloc(i-1))
        value(i) = value(i-1);
    end
end
subplot 529;plot(value)
L4 = interp1(vloc, value, 1:N, "spline");
L4 = L4 - (L4(1)-d(1));
subplot 529; plot(t,L4,t,d); grid on; axis tight; title("Even-power Interpolation Error = " + num2str(RMSE(L4,d)));
subplot (5,2,10); plot(t,L4-d);grid on; axis tight;

%% Even-power Interpolation Comparison
p0 = p.^4 - p.^2;
p1 = p0.^2 + p0/(2.^2);
p2 = p1.^2 + p1/(2.^6);
p3 = p2.^2 + p2/(2.^14);

pn = (p3 + 1/(2.^31)) .* (2.^31); 
times = 2.^5;

dpx = (lambda0/2)/times;
[~, vloc] = findpeaks(-pn, "MinPeakHeight",0.9);
vloc(ismember(vloc, find(abs(diff(direction))==2))) = [];
value = zeros(1,length(vloc));

for i = 2:length(value)
    if direction(vloc(i))==direction(vloc(i-1))
        value(i) = value(i-1) + dpx*direction(vloc(i));
    elseif direction(vloc(i))~=direction(vloc(i-1))
        value(i) = value(i-1);
    end
end
L4 = interp1(vloc, value, 1:N, "spline");
L4 = L4 - (L4(1)-d(1));
subplot 529; plot(t,L4,t,d); grid on; axis tight; title("Even-power Interpolation Error = " + num2str(RMSE(L4,d)));
subplot (5,2,10); plot(t,L4-d);grid on; axis tight;
