%% 五行一列图窗,扁图A_API_6row1col(a,b,c,d,e,f)
    fontName = 'Arial';
    fontSize = 8;
    fontWeight = 'bold';
    lineWidth = 1.5; 
    set(gcf,'Position',[259.4,243.4,840,500])  % 设置最舒服的图窗大小
    %% 第1幅图
    axes('Position',[0.057 0.8595 0.88 0.118]);  axis tight;  % 设置图像位置 + 坐标紧凑型 
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(a)', 'Units', 'normalized', 'Position', [0.0222,0.566,0],'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    
    
    %% 第2幅图
    axes('Position',[0.057 0.666 0.88 0.118]);  axis tight; 
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(a)', 'Units', 'normalized', 'Position', [0.0222,0.566,0],'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级

    %% 第3幅图
    axes('Position',[0.057 0.479 0.88 0.118]);  axis tight;
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(a)', 'Units', 'normalized', 'Position', [0.0222,0.566,0],'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级

    %% 第4幅图
    axes('Position',[0.057 0.275 0.88 0.118]);  axis tight;  % 设置图像位置 + 坐标紧凑型 
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(a)', 'Units', 'normalized', 'Position', [0.0222,0.566,0],'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    
    %% 第5幅图
    axes('Position',[0.06 0.076 0.88 0.118]);  axis tight;
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(a)', 'Units', 'normalized', 'Position', [0.0222,0.566,0 ],'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级

end

