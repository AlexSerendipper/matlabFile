%% 六行一列图窗
function DRWA_API_6row1col(a,b,c,d,e,f)
    fontName = '微软雅黑';
    fontSize = 11;
    fontWeight = 'bold';
    lineWidth = 1.5; 
    set(gcf,'Position',[1,41,1536,749.6])  % 设置最舒服的图窗大小
    %% 位置信息    
    x = 0.88;  % 大小
    y = 0.12;  % 大小
    m1 = [0.057,0.8595];
    m2 = [0.057,0.700668943436499];
    m3 = [0.057,0.542005122732124];
    m4 = [0.057,0.383575773745998];
    m5 = [0.06,0.225077481323373];
    m6 = [0.0615625,0.069248132337247];
    %% 第1幅图
    axes('Position',[m1(1) m1(2) x y]);  axis tight;  % 设置图像位置 + 坐标紧凑型 
    
    plot([1,2,3],'LineWidth',lineWidth);
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
%     xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(a)', 'Units', 'normalized', 'Position', [0.0278,0.6259,0],'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    
    
    %% 第2幅图
    axes('Position',[m2(1) m2(2) x y]);  axis tight; 
    
    plot([1,2,3],'LineWidth',lineWidth);
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
%     xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(b)', 'Units', 'normalized', 'Position', [0.0278,0.6259,0],'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级

    %% 第3幅图
    axes('Position',[m3(1) m3(2) x y]);  axis tight;
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
%     xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(a)', 'Units', 'normalized', 'Position', [0.0278,0.6259,0],'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级

    %% 第4幅图
    axes('Position',[m4(1) m4(2) x y]);  axis tight;  % 设置图像位置 + 坐标紧凑型 
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
%     xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(d)', 'Units', 'normalized', 'Position', [0.0278,0.6259,0],'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    
    %% 第5幅图
    axes('Position',[m5(1) m5(2) x y]);  axis tight;
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
%     xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(e)', 'Units', 'normalized', 'Position', [0.0278,0.6259,0],'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级

    %% 第6幅图
    axes('Position',[m6(1) m6(2) x y]);  axis tight;  % 设置图像位置 + 坐标紧凑型 
    plot([1,2,3],'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
     xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(f)', 'Units', 'normalized', 'Position', [0.0278,0.6259,0],'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
end

