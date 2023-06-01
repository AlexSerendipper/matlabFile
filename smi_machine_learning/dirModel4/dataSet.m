% 数据集
fringeData = [];
for C=0.05:0.001:0.05
    [data] = SMI_API_FringeData(C,5)
    fringeData = [fringeData;data];
end


