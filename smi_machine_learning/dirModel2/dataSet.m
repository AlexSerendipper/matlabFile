% 数据集
fringeData = [];
for C=0.05:0.001:5
    [data] = SMI_API_FringeData(C);
    fringeData = [fringeData;data];
end

