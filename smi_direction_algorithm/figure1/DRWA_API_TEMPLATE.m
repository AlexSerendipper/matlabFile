function DRWA_API_TEMPLATE(varargin)
    % 获取当前文件夹下的所有文件
    files = dir('*.mat');
    for i = 1:length(files)
        % 加载所有.mat文件
        load(files(i).name);
    end
    
    %% 五行一列图窗,扁图A_API_6row1col(a,b,c,d,e,f)
    fontName = '微软雅黑';
    fontSize = 10;
    fontWeight = 'bold';
    lineWidth = 1.5; 
    %% 位置信息
    % 仅一张图
    if length(varargin)-3 == 1
        x = varargin{1} ;  % 大小
        y = varargin{2};  % 大小
        z = varargin{3}; % 角标位置        
        m1 = cell2mat(varargin(4));
    end
    if length(varargin)-3 == 2
        x = varargin{1} ;  % 大小
        y = varargin{2};  % 大小
        z = varargin{3}; % 角标位置        
        m1 = cell2mat(varargin(4));
        m2 = cell2mat(varargin(5));
    end
    if length(varargin)-3 == 3
        x = varargin{1} ;  % 大小
        y = varargin{2};  % 大小
        z = varargin{3}; % 角标位置        
        m1 = cell2mat(varargin(4));
        m2 = cell2mat(varargin(5));
        m3 = cell2mat(varargin(6));
    end
    if length(varargin)-3 == 4
        x = varargin{1} ;  % 大小
        y = varargin{2};  % 大小
        z = varargin{3}; % 角标位置        
        m1 = cell2mat(varargin(4));
        m2 = cell2mat(varargin(5));
        m3 = cell2mat(varargin(6));
        m4 = cell2mat(varargin(7));
    end
    if length(varargin)-3 == 5
        x = varargin{1} ;  % 大小
        y = varargin{2};  % 大小
        z = varargin{3}; % 角标位置        
        m1 = cell2mat(varargin(4));
        m2 = cell2mat(varargin(5));
        m3 = cell2mat(varargin(6));
        m4 = cell2mat(varargin(7));
        m5 = cell2mat(varargin(8));
    end
    if length(varargin)-3 == 6
        x = varargin{1} ;  % 大小
        y = varargin{2};  % 大小
        z = varargin{3}; % 角标位置        
        m1 = cell2mat(varargin(4));
        m2 = cell2mat(varargin(5));
        m3 = cell2mat(varargin(6));
        m4 = cell2mat(varargin(7));
        m5 = cell2mat(varargin(8));
        m6 = cell2mat(varargin(9));
    end    
    if length(varargin)-3 == 7
        x = varargin{1} ;  % 大小
        y = varargin{2};  % 大小
        z = varargin{3}; % 角标位置        
        m1 = cell2mat(varargin(4));
        m2 = cell2mat(varargin(5));
        m3 = cell2mat(varargin(6));
        m4 = cell2mat(varargin(7));
        m5 = cell2mat(varargin(8));
        m6 = cell2mat(varargin(9));
        m7 = cell2mat(varargin(10));
    end 
    if length(varargin)-3 == 8
        x = varargin{1} ;  % 大小
        y = varargin{2};  % 大小
        z = varargin{3}; % 角标位置        
        m1 = cell2mat(varargin(4));
        m2 = cell2mat(varargin(5));
        m3 = cell2mat(varargin(6));
        m4 = cell2mat(varargin(7));
        m5 = cell2mat(varargin(8));
        m6 = cell2mat(varargin(9));
        m7 = cell2mat(varargin(10));
        m8 = cell2mat(varargin(11));
    end 
    %% 第1幅图
    axes('Position',[m1(1) m1(2) x y]);  axis tight;  % 设置图像位置 + 坐标紧凑型 
    plot(t,p_init,'LineWidth',lineWidth);
    hold on;
%     scatter(loc_op0/fs,top_op0','r','LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Amp.(a.u)')
    ylabel('Amp.(a.u)')
    grid on;
%     axis tight; 
    title('(a)', 'Units', 'normalized','FontWeight','bold','color', 'black','Position', z ,'FontSize', 16);
    set(gca,'YTickLabel')  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 1
        return;
    end
    %% 第2幅图
    axes('Position',[m2(1) m2(2) x y]);  axis tight; 
    plot(c,'LineWidth',lineWidth);
    hold on;
    scatter(loc_c,top_c);
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Amp.(a.u)')
    ylabel('Amp.(a.u)')
    grid on;
%     axis tight; 
    title('(b)', 'Units', 'normalized','FontWeight','bold','color', 'black','Position', z ,'FontSize', 16);
    set(gca,'YTickLabel')  % 乘10^6，让图像的纵坐标显示的是微米量级
     set(gca,'XTick',[0,800,800*2,800*3,800*4,800*5,800*6,800*7,800*8,800*9,800*10])
    set(gca,'XTickLabel',[0,0.002,0.002*2,0.002*3,0.002*4,0.002*5,0.002*6,0.002*7,0.002*8,0.002*9,0.002*10])
    if nargin - 3== 2
        return;
    end
    %% 第3幅图
    axes('Position',[m3(1) m3(2) x y]);  axis tight;
    plot(t,p_init,'LineWidth',lineWidth);
    hold on;
    scatter(loc_op1/fs,top_op1,'p','filled','r','LineWidth',lineWidth);

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Amp.(a.u)')
    ylabel('Amp.(a.u)')
    grid on;
%     axis tight; 
    title('(c)', 'Units', 'normalized','FontWeight','bold','color', 'black','Position', z ,'FontSize', 16);
    set(gca,'YTickLabel')  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 3
        return;
    end
    %% 第4幅图
    axes('Position',[m4(1) m4(2) x y]);  axis tight;  % 设置图像位置 + 坐标紧凑型 
    plot(t,p_init,'LineWidth',lineWidth);
    hold on;
    scatter(loc_op2/fs,top_op2,'p','filled','r','LineWidth',lineWidth);

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Amp.(a.u)')
    ylabel('Amp.(a.u)')
    grid on;
%     axis tight; 
    title('(d)', 'Units', 'normalized','FontWeight','bold','color', 'black','Position', z ,'FontSize', 16);
    set(gca,'YTickLabel')  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 4
        return;
    end
    %% 第5幅图
    axes('Position',[m5(1) m5(2) x y]);  axis tight;  % 设置图像位置 + 坐标紧凑型 
    plot(t,p_init,'LineWidth',lineWidth);
    hold on;
    scatter(loc_op3/fs,top_op3,'p','filled','r','LineWidth',lineWidth);

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Amp.(a.u)')
    ylabel('Amp.(a.u)')
    grid on;
%     axis tight; 
    title('(e)', 'Units', 'normalized','FontWeight','bold','color', 'black','Position', z ,'FontSize', 16);
    set(gca,'YTickLabel')  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 5
        return;
    end
    %% 第6幅图
    axes('Position',[m6(1) m6(2) x y]);  axis tight;
    plot(t,p_init,'LineWidth',lineWidth);
    hold on;
    scatter(loc_op4/fs,top_op4,'p','filled','r','LineWidth',lineWidth);

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Amp.(a.u)')
    grid on;
%     axis tight; 
    title('(f)', 'Units', 'normalized','FontWeight','bold','color', 'black','Position', z ,'FontSize', 16);
    set(gca,'YTickLabel')  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 6
        return;
    end
    %% 第7幅图
    axes('Position',[m7(1) m7(2) x y]);  axis tight;
    plot([1,2,3],'LineWidth',lineWidth);

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(g)', 'Units', 'normalized', 'Position', z ,'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 7
        return;
    end
    %% 第8幅图
    axes('Position',[m8(1) m8(2) x y]);  axis tight;
    plot([1,2,3],'LineWidth',lineWidth);

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(h)', 'Units', 'normalized', 'Position', z ,'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 8
        return;
    end
end
