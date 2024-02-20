%% 四行一列图窗
function DRWA_API_4row1col(a,b,c,d)
    % 第1幅图
    axes('Position',[0.06,0.702,0.401,0.232]);  axis tight;  % 设置图像位置 + 坐标紧凑型
    
    load("D:\matlab save\self-mixing\smi_MODULATE\temp_phiF_wrapped_times1.mat");
    plot(phiF_wrapped_init); grid on;  % 作图 + 网格
    
    ylabel('Amp.(a.u)'); title('(a)', 'Units', 'normalized','FontSize', 16, 'Position', [0.035 0.7]);  % y轴注释 + 左上角角标，小标位置是相对当前图而言
%     set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    % legend("C = 0.03");
    
    % 第2幅图
    axes('Position',[0.538,0.702,0.401,0.232]);   % 设置图像位置 + 坐标紧凑型
    
    load("D:\matlab save\self-mixing\smi_MODULATE\temp_phiF_wrapped_times2.mat");
    plot(phiF_wrapped_init); grid on;  % 作图 + 网格
    
    axis tight;
    ylabel('Amp.(a.u)'); title('(b)', 'Units', 'normalized','FontSize', 16, 'Position', [0.035 0.7]);  % y轴注释 + 左上角角标，小标位置是相对当前图而言
%     set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    % legend("C = 0.03");
    
    % 第3幅图
    axes('Position',[0.06, 0.4 ,0.401,0.232]);   
    load("D:\matlab save\self-mixing\smi_MODULATE\temp_recon._times1.mat");
    plot(Lt);
    hold on;
    plot(Lt_reconstruct_times1); grid on;
    
%     axis tight;
    ylabel('Displacement.(um)'); title('(c)', 'Units', 'normalized','FontSize', 16, 'Position', [0.035 0.7]); 
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    legend("real","OSPM(N=1)");

    % 第4幅图
    axes('Position',[0.538 0.4 ,0.401,0.232]);   
    load("D:\matlab save\self-mixing\smi_MODULATE\temp_recon._times2.mat");
    plot(Lt);
    hold on;
    plot(Lt_reconstruct_times2); grid on;
    
    axis tight;
    ylabel('Displacement.(um)'); title('(d)', 'Units', 'normalized','FontSize', 16, 'Position', [0.035 0.7]); 
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    legend("REAL","OSPM(N=2)");
    
    % 第5幅图
    axes('Position',[0.06 0.1 0.401,0.232]);  
    
%     plot(Lt(padding:end-padding+1)-Lt_reconstruct_times1(padding:end-padding+1)); grid on;  % 作图 + 网格
    plot(Lt-Lt_reconstruct_times1); grid on;  % 作图 + 网格
    
    axis tight;
    ylabel('Error.(nm)'); title('(e)', 'Units', 'normalized', 'FontSize', 16, 'Position', [0.035 0.7]); 
    set(gca,'YTickLabel',10^9*get(gca,'YTick'))  % 乘10^3，让图像的纵坐标显示的是毫米量级
    % legend("Reconstruction","Reference");

    % 第6幅图
    axes('Position',[0.538 0.1 0.401,0.232]);
    
%     plot(Lt(padding:end-padding+1)-Lt_reconstruct_times2(padding:end-padding+1)); grid on;  % 作图 + 网格
    plot(Lt-Lt_reconstruct_times2); grid on;  % 作图 + 网格
    axis tight;
    ylabel('Error.(nm)'); title('(f)', 'Units', 'normalized', 'FontSize', 16, 'Position', [0.035 0.7]); 
    set(gca,'YTickLabel',10^9*get(gca,'YTick'))  % 乘10^3，让图像的纵坐标显示的是毫米量级
    axis tight;
    % legend("Reconstruction","Reference");
end

