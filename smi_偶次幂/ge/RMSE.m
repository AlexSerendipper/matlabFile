function x = RMSE(yr, yi)
% 计算向量的均方误差
% yr:需要被计算的向量
% yi:定标

x = sqrt( sum( (yr - yi).^2 ) / length( yr ) );

end