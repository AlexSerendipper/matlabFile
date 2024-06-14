%% 画频谱图时，fshift如果0点不在最中心，可以调整xlabel，加一行再减一行后可以调整到中心
function DRWA_API_TEMPLATE(varargin)
    % 获取当前文件夹下的所有文件
    files = dir('*.mat');
    for i = 1:length(files)
        % 加载所有.mat文件
        load(files(i).name);
    end
    
     % 根据需要设置最舒服的图窗大小
%      set(gcf,'Position',[200.2,193,1149.6,463.2]);
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
    if length(varargin)-3 == 10
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
        m9 = cell2mat(varargin(12));
        m10 = cell2mat(varargin(13));
    end 
    z1 = [0.062564826175869,0.601803614457831,0];
    %% 第1幅图
    axes('Position',[m1(1) 0.8763 x y]);  % 设置图像位置 
    plot(t,p_init,'LineWidth',lineWidth);
    hold on;
    scatter(loc_op./fs,top_op,'p','filled','k','LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Amp.(a.u.)')
    grid on;
%   box off;  %关闭右面和上面的坐标轴
%   axis tight; % 设置坐标紧凑型
%   axis([0 0.08 -20*10^-9 20*10^-9]);  % 设置具体横坐标和纵坐标的范围[x1,x2,y1,y2]
    title('(a)', 'Units', 'normalized', 'FontWeight','bold', 'color', 'black', 'Position', z ,'FontSize', 16);
    set(gca,'YTickLabel',1*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
%   set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')  % 不显示x轴和y轴
    if nargin - 3  == 1
        return;
    end
    %% 第2幅图
    axes('Position',[0.076,0.710871970624234,0.403012345679012,0.1025]);
    len = length(F);
    TF = abs(TF);
    mesh(T,F(len/2-100:len/2+100,:),TF(len/2-100:len/2+100,:));
%     mesh(abs(TF_curb)); 
    view(0,90); % 设置初始视角为俯视角

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
%   box off;  %关闭右面和上面的坐标轴
  axis tight; % 设置坐标紧凑型
    title('(b)', 'Units', 'normalized','FontWeight','bold','color', 'white','Position', z1 ,'FontSize', 16);
    set(gca,'YTickLabel',10^-3*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3== 2
        return;
    end
    %% 第3幅图
    axes('Position',[m3(1) m3(2) x y]);
    plot(t,p,'LineWidth',lineWidth);
    hold on;
    scatter(loc_p./fs,top_p,'p','filled','r','LineWidth',lineWidth);
    hold on;
    plot(t,ddir,'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Amp.(a.u.)')
    grid on;
%   axis tight; % 设置坐标紧凑型
    title('(d)', 'Units', 'normalized','FontWeight','bold','color', 'black', 'Position', z ,'FontSize', 16);
    set(gca,'YTickLabel',1*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 3
        return;
    end
    %% 第4幅图
    axes('Position',[m4(1) m4(2) x y]);  % 设置图像位置
    plot(t,direction,'LineWidth',lineWidth);

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Amp.(a.u.)')
    grid on;
%   axis tight; % 设置坐标紧凑型
    axis([0 0.04 -1.2 1.2]);
    title('(e)', 'Units', 'normalized', 'FontWeight','bold','color', 'black','Position',z ,'FontSize', 16);
    set(gca,'YTickLabel',1*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 4
        return;
    end
    %% 第5幅图
    axes('Position',[m5(1) m5(2) x y]);
    [Lt] = MOVE_API_STANDARD(220, 100, 125000, 5000, 11.5, 0);
    plot(t,Lt,'LineWidth',lineWidth);
    hold on;
    plot(t,Lt_reconstruct,'LineStyle','--','LineWidth',lineWidth);
    legend("参考信号","重构信号")
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Lt(\mum)')
    grid on;
%   axis tight; % 设置坐标紧凑型
    title('(f)', 'Units', 'normalized','FontWeight','bold','color', 'black', 'Position',z ,'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 5
        return;
    end
    %% 第6幅图
    axes('Position',[m6(1) m6(2) x y]);
    RMSE = sqrt(mean((Lt-Lt_reconstruct).^2));
    plot(t,Lt_reconstruct-Lt,'LineWidth',lineWidth);
    hold on;
    plot(t,31.74*ones(1,N)*10^-9,'LineWidth',lineWidth);
%     title(['RMSE=', num2str(RMSE)])
    legend("","RMSE");

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Error(nm)')
    grid on;
%   axis tight; % 设置坐标紧凑型
    title('(g)', 'Units', 'normalized','FontWeight','bold','color', 'black', 'Position', z ,'FontSize', 16);
    set(gca,'YTickLabel',10^9*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 7
        return;
    end
    %% 第7幅图
    axes('Position',[0.546502057613169,0.710871970624234,0.403012345679012,0.1025]);
    len = length(F);
    TF_curb = abs(TF_curb);
    mesh(T,F(len/2-100:len/2+100,:),TF_curb(len/2-100:len/2+100,:));
%     mesh(abs(TF_curb)); 
    view(0,90); % 设置初始视角为俯视角

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
  axis tight; % 设置坐标紧凑型
    title('(c)', 'Units', 'normalized','FontWeight','bold','color', 'white', 'Position', z1 ,'FontSize', 16);
    set(gca,'YTickLabel',10^-3*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 6
        return;
    end
    %% 第8幅图
    axes('Position',[m8(1) m8(2) x y]);
    plot([1,2,3],'LineWidth',lineWidth);

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
%   axis tight; % 设置坐标紧凑型
    title('(h)', 'Units', 'normalized', 'FontWeight','bold','color', 'black', 'Position', z ,'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 8
        return;
    end
    %% 第9幅图
    axes('Position',[m9(1) m9(2) x y]);
    plot([1,2,3],'LineWidth',lineWidth);

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
%   axis tight; % 设置坐标紧凑型
    title('(i)', 'Units', 'normalized', 'FontWeight','bold','color', 'black', 'Position', z ,'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 8
        return;
    end
    %% 第10幅图
    axes('Position',[m10(1) m10(2) x y]);
    plot([1,2,3],'LineWidth',lineWidth);

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
%   axis tight; % 设置坐标紧凑型
    title('(j)', 'Units', 'normalized', 'FontWeight','bold','color', 'black', 'Position', z ,'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 8
        return;
    end
end

