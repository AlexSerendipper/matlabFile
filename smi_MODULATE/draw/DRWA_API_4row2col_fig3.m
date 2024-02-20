
    % 第1幅图
    axes('Position',[0.06,0.775,0.401,0.16]);  axis tight;  % 设置图像位置 + 坐标紧凑型
    
    plot(p(100:13400)); grid on;  % 作图 + 网格
    
    ylabel('Amp.(a.u)'); title('(a)', 'Units', 'normalized','FontSize', 16, 'Position', [0.035 0.7]);  % y轴注释 + 左上角角标，小标位置是相对当前图而言
%     set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    % legend("C = 0.03");
    
    % 第2幅图
    axes('Position',[0.538,0.775,0.401,0.16]);  axis tight;  % 设置图像位置 + 坐标紧凑型
    
    plot(amp1(18000:18000+13301)); grid on;  % 作图 + 网格
    
    ylabel('Freq.(Khz)'); title('(b)', 'Units', 'normalized','FontSize', 16, 'Position', [0.035 0.7]);  % y轴注释 + 左上角角标，小标位置是相对当前图而言
    set(gca,'YTickLabel',10^-3*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    % legend("C = 0.03");
    
    % 第3幅图
    axes('Position',[0.06 0.55 ,0.401 0.16]);  axis tight; 
    
    plot(p1(100:13400)); grid on;
    
    ylabel('Amp.(a.u) '); title('(c)', 'Units', 'normalized','FontSize', 16, 'Position', [0.035 0.7]); 
    % set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    % legend("C = 0.03");

    % 第4幅图
    axes('Position',[0.538 0.55 ,0.401 0.16]);  axis tight; 
    
    plot(p2(100:13400)); grid on;
    
    ylabel('Amp.(a.u) '); title('(d)', 'Units', 'normalized','FontSize', 16, 'Position', [0.035 0.7]); 
    % set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    % legend("C = 0.03");
    
    % 第5幅图
    axes('Position',[0.06 0.325 0.401 0.16]);  axis tight;
    
    load("D:\matlab save\self-mixing\smi_MODULATE\temp_experi._phiF_wrapped_times3.mat")
    plot(phiF_wrapped_init(100:13400)); grid on;
    
    ylabel('Amp.(a.u)'); title('(e)', 'Units', 'normalized', 'FontSize', 16, 'Position', [0.035 0.7]); 
    % set(gca,'YTickLabel',10^3*get(gca,'YTick'))  % 乘10^3，让图像的纵坐标显示的是毫米量级
    % legend("Reconstruction","Reference");

    % 第6幅图
    axes('Position',[0.538 0.325 0.401 0.16]);  axis tight;
    
    load("D:\matlab save\self-mixing\smi_MODULATE\temp_experi._phiF_wrapped_times4.mat")
    plot(phiF_wrapped_init(100:13400)); grid on;
    
    ylabel('Amp.(a.u)'); title('(f)', 'Units', 'normalized', 'FontSize', 16, 'Position', [0.035 0.7]); 
    % set(gca,'YTickLabel',10^3*get(gca,'YTick'))  % 乘10^3，让图像的纵坐标显示的是毫米量级
    % legend("Reconstruction","Reference");
    
    % 第7幅图
    axes('Position',[0.06 0.1 0.401 0.16]);  axis tight;  % 设置图像位置 + 坐标紧凑型 
    plot(Lt_reconstruct1(100:13400));
    hold on;
    load("D:\matlab save\self-mixing\smi_MODULATE\temp_experi._recon._times3.mat");
    plot(Lt_reconstruct_times3(100:13400),'Color','r'); grid on;  % 作图 + 网格
    hold on;
    load("D:\matlab save\self-mixing\smi_MODULATE\temp_experi._recon._times4.mat");
    plot(Lt_reconstruct_times4(100:13400),'Color','g');
    legend("reference","OSPM(N=3)","OSPM(N=4)");
    ylabel('Displacement(um)'); title('(g)', 'Units', 'normalized', 'FontSize', 16, 'Position', [0.035 0.7]);  % y轴注释 + 左上角角标，小标位置是相对当前图而言
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    
    % 第8幅图
    axes('Position',[0.538 0.1 0.401 0.16]);  
    load("D:\matlab save\self-mixing\smi_MODULATE\temp_experi._recon._times1.mat");
%   Lt_reconstruct_times2(4263:4890) = Lt_reconstruct_times2(1213:1213+627);
    plot(Lt_reconstruct1(100:13400)-Lt_reconstruct_times3(100:13400)); grid on;  % 作图 + 网格
    hold on;
    plot(Lt_reconstruct1(100:13400)-Lt_reconstruct_times4(100:13400));
    ylabel('Error(nm)'); title('(h)', 'Units', 'normalized', 'FontSize', 16, 'Position', [0.035 0.7]);  % y轴注释 + 左上角角标，小标位置是相对当前图而言
    set(gca,'YTickLabel',10^9*get(gca,'YTick'))  % 乘10^9，让图像的纵坐标显示的是微米量级
    legend("OSPM(N=3)","OSPM(N=4)");
    axis tight;  % 设置图像位置 + 坐标紧凑型 


