clc;
clear all;
close all;

% 数据集
fringeData = [];

for C=0.02:0.01:3
    for alpha=4:0.1:6
        [data] = SMI_API_FringeData(C,alpha);
        fringeData = [fringeData;data];
    end
end


