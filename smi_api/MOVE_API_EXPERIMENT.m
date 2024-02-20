%% 本程序针对老电脑实验信号的csv文件的读取
%% 从第一列第三行开始往下都是时间（ns），第二列第三行开始往下都是对应数据值
function [t, p, fs] = MOVE_API_EXPERIMENT(M, N, path)   % fs为采样率(s)
    data = readmatrix(path);  % 输入文件路径读取文件
    t = data(M:M+N-1 ,1)';  % 第一列，从第M个点开始读取N个数据
    p = data(M:M+N-1 ,2)';  
    fs = 1 / ((t(2) - t(1))* 10^-9);  % 采样率，从采样时间(ns)中算出
    
    t = (0:N-1)/fs;  % 保证t是从0开始，要把之前的覆盖掉
end


