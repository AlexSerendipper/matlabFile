clc;
clear all;
close all;

% 数据集
fringeData = [];

for C=0.1:0.1:0.5
    for alpha=4:0.1:4
        [data] = SMI_API_FringeData(C,alpha);
        fringeData = [fringeData;data];
    end
end


