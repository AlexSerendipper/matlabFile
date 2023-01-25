%% 四行一列图窗，论文图1，翻倍一次的时候跑这个程序
function DRWA_API_4row1col(a,b,c,d)
    % 第1幅图
    axes('Position',[0.06,0.775,0.401,0.16]);  axis tight;  % 设置图像位置 + 坐标紧凑型
    
    plot(p); grid on;  % 作图 + 网格
    
    ylabel('Amp.(a.u)'); title('(a)', 'Units', 'normalized','FontSize', 16, 'Position', [0.035 0.7]);  % y轴注释 + 左上角角标，小标位置是相对当前图而言
%     set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    % legend("C = 0.03");
    
    % 第2幅图
    axes('Position',[0.538,0.775,0.401,0.16]);  axis tight;  % 设置图像位置 + 坐标紧凑型
    
    plot(amp1); grid on;  % 作图 + 网格
    
    ylabel('Freq.(Khz)'); title('(b)', 'Units', 'normalized','FontSize', 16, 'Position', [0.035 0.7]);  % y轴注释 + 左上角角标，小标位置是相对当前图而言
    set(gca,'YTickLabel',10^-3*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    % legend("C = 0.03");
    
    % 第3幅图
    axes('Position',[0.06 0.55 ,0.401 0.16]);  axis tight; 
    
    plot(p1); grid on;
    
    ylabel('Amp.(a.u) '); title('(c)', 'Units', 'normalized','FontSize', 16, 'Position', [0.035 0.7]); 
    % set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    % legend("C = 0.03");

    % 第4幅图
    axes('Position',[0.538 0.55 ,0.401 0.16]);  axis tight; 
    
    plot(p2); grid on;
    
    ylabel('Amp.(a.u) '); title('(d)', 'Units', 'normalized','FontSize', 16, 'Position', [0.035 0.7]); 
    % set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    % legend("C = 0.03");
    
    % 第5幅图
    axes('Position',[0.06 0.325 0.401 0.16]);  axis tight;
    plot(Pn); grid on;
    ylabel('Amp.(a.u)'); title('(e)', 'Units', 'normalized', 'FontSize', 16, 'Position', [0.035 0.7]); 
    % set(gca,'YTickLabel',10^3*get(gca,'YTick'))  % 乘10^3，让图像的纵坐标显示的是毫米量级
    % legend("Reconstruction","Reference");

    % 第6幅图
    axes('Position',[0.538 0.325 0.401 0.16]);  axis tight;
    
    plot(Pnn); grid on;
    
    ylabel('Amp.(mm)'); title('(f)', 'Units', 'normalized', 'FontSize', 16, 'Position', [0.035 0.7]); 
    % set(gca,'YTickLabel',10^3*get(gca,'YTick'))  % 乘10^3，让图像的纵坐标显示的是毫米量级
    % legend("Reconstruction","Reference");
    
    % 第7幅图
    axes('Position',[0.06 0.1 0.401 0.16]);  axis tight;  % 设置图像位置 + 坐标紧凑型 
    plot(Lt);
    hold on;
    load("D:\matlab save\self-mixing\smi_MODULATE\temp_recon._times0.mat");
    plot(Lt_reconstruct_times0,'r'); grid on;  % 作图 + 网格
    hold on;
    load("D:\matlab save\self-mixing\smi_MODULATE\temp_recon._times1.mat");
    plot(Lt_reconstruct_times1,'g');
    legend("real","EOM","OSPM");
    ylabel('Displacement(um)'); title('(g)', 'Units', 'normalized', 'FontSize', 16, 'Position', [0.035 0.7]);  % y轴注释 + 左上角角标，小标位置是相对当前图而言
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    
    % 第8幅图
    axes('Position',[0.538 0.1 0.401 0.16]);   % 设置图像位置 + 坐标紧凑型 
    
    plot(Lt(padding:end-padding+1)-Lt_reconstruct_times0(padding:end-padding+1)); grid on;  % 作图 + 网格
    hold on;
    plot(Lt(padding:end-padding+1)-Lt_reconstruct_times1(padding:end-padding+1));
    
    axis tight;
    ylabel('Error(nm)'); title('(h)', 'Units', 'normalized', 'FontSize', 16, 'Position', [0.035 0.7]);  % y轴注释 + 左上角角标，小标位置是相对当前图而言
    set(gca,'YTickLabel',10^9*get(gca,'YTick'))  % 乘10^9，让图像的纵坐标显示的是微米量级
    legend("EOM","OSPM")
end

