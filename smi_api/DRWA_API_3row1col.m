%% 这个要放小比较好看
function DRWA_API_2row2col(c)
    fontName = 'Arial';
    fontSize = 8;
    fontWeight = 'bold';
    lineWidth = 1.5;
    %% 第1幅图
    axes('Position',[0.06,0.73,0.878,0.22]);
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;  % 作图 + 网格
    axis tight;  % 坐标紧凑型
    title('(a)', 'Units', 'normalized', 'Position', [0.033,0.656,0],'FontSize', 16);  % 左上角角标，小标位置是相对当前图而言
 
    %% 第2幅图
    axes('Position',[0.06,0.413,0.878,0.22]);
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(c)', 'Units', 'normalized', 'Position', [0.033,0.656,0],'FontSize', 16);
    
    %% 第3幅图
    axes('Position',[0.06,0.098,0.878,0.22]);
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(c)', 'Units', 'normalized', 'Position', [0.033,0.656,0],'FontSize', 16);
end

