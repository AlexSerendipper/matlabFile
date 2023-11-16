clc;
clear all;
close all;

% 数据集
fringeData = [];
alpha = 5;
for C=0.05:0.01:1.5
    [data] = SMI_API_FringeData(C,alpha);
    fringeData = [fringeData;data];
end



