%% 六行一列图窗
function DRWA_API_6row1col(a,b,c,d,e,f)
    fontName = 'Arial';
    fontSize = 8;
    fontWeight = 'bold';
    lineWidth = 1.5; 
    %% 第1幅图
    axes('Position',[0.057 0.8595 0.88 0.118]);  axis tight;  % 设置图像位置 + 坐标紧凑型 
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(a)', 'Units', 'normalized', 'Position', [0.0278,0.6259,0],'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    
    
    %% 第2幅图
    axes('Position',[0.057 0.7007 0.88 0.118]);  axis tight; 
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(a)', 'Units', 'normalized', 'Position', [0.0278,0.6259,0],'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级

    %% 第3幅图
    axes('Position',[0.057 0.5367 0.88 0.118]);  axis tight;
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(a)', 'Units', 'normalized', 'Position', [0.0278,0.6259,0],'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级

    %% 第4幅图
    axes('Position',[0.057 0.3729 0.88 0.118]);  axis tight;  % 设置图像位置 + 坐标紧凑型 
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(a)', 'Units', 'normalized', 'Position', [0.0278,0.6259,0],'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    
    %% 第5幅图
    axes('Position',[0.06 0.2048 0.88 0.118]);  axis tight;
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(a)', 'Units', 'normalized', 'Position', [0.0278,0.6259,0],'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级

    %% 第6幅图
    axes('Position',[0.06 0.0415 0.88 0.118]);  axis tight;  % 设置图像位置 + 坐标紧凑型 
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(a)', 'Units', 'normalized', 'Position', [0.0278,0.6259,0],'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
end

