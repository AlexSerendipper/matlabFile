%% 知道一些参数，可以用于定标。其中phi和offset用1 2 3..逐渐尝试即可，phi是左右，越大越左，offset是上下
function [Lt] = MOVE_API_STANDARD(mv, fv, fs, N, phi, offset)
    A = 81*mv/7*10^-9;  % 7mv对应80-82nm，输入mv得到对应A值
    t = (0:N-1)/fs;  % 采样时间，设N=10, fs=200，即采样了0.05s，t为[0...0.045]
    Lt = A.* sin(2*pi*fv*t + phi*pi/64)+offset*650*10^-9;  % 以pi/64为步进进行相移，慢慢调试
                                                         % 以波长为步进调试上下偏移量
end