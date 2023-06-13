% 数据集
clc;
clear all;
close all;
fringeData = [];
SNR = 20;
for C=0.05:0.001:2.8
    if C<=1
        HP = 1.4;
        [data] = SMI_API_FringeData(C,SNR,HP);
        fringeData = [fringeData;data];
    else
        HP = 1.25;
        [data] = SMI_API_FringeData(C,SNR,HP);
        fringeData = [fringeData;data];
    end
end


