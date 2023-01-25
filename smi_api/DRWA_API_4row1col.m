%% 四行一列图窗
function DRWA_API_4row1col(a,b,c,d)
    fontName = 'Arial';
    fontSize = 8;
    fontWeight = 'bold';
    lineWidth = 1.5;   
    %% 第1幅图
    axes('Position',[0.06 0.775 0.88 0.16]);  axis tight;  % 设置图像位置 + 坐标紧凑型
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(a)', 'Units', 'normalized', 'Position', [0.0342,0.652,0],'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级

    %% 第2幅图
    axes('Position',[0.06 0.55 0.88 0.16]);  axis tight; 
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(a)', 'Units', 'normalized', 'Position', [0.0342,0.652,0],'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级

    %% 第3幅图
    axes('Position',[0.06 0.325 0.88 0.16]);  axis tight;
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(a)', 'Units', 'normalized', 'Position', [0.0342,0.652,0],'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级

    %% 第四幅图
    axes('Position',[0.06 0.1 0.88 0.16]);  axis tight;  % 设置图像位置 + 坐标紧凑型 
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(a)', 'Units', 'normalized', 'Position', [0.0342,0.652,0],'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
end

