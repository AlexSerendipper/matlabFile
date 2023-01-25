%% 一行两列图窗
function DRWA_API_1row2col(a,b)
    fontName = 'Arial';
    fontSize = 8;
    fontWeight = 'bold';
    lineWidth = 1.5;
    %% 画图：第1幅图
    axes('Position',[0.088 0.145 0.383 0.746]);   
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;  % 作图 + 网格
    axis tight;  % 坐标紧凑型
    title('(a)', 'Units', 'normalized', 'Position', [0.097,0.893,0],'FontSize', 16);  % 左上角角标，小标位置是相对当前图而言
    % set(gca,'YTickLabel',10^-3*get(gca,'YTick'))  % 乘10^3，让图像的纵坐标显示的是Khz
    % set(gca, 'xtick', 1:5);
    % legend("C = 0.03");
    
    %% 第2幅图
    axes('Position',[0.555 0.145 0.383 0.746]);
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(b)', 'Units', 'normalized', 'Position', [0.097,0.893,0],'FontSize', 16);
end

