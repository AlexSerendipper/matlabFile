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
    %% 第1幅图
    axes('Position',[m1(1) m1(2) x y]);  % 设置图像位置 
    plot(t,p_init,'LineWidth',lineWidth);
    hold on;
    scatter(loc_op./fs,top_op,'p','filled','k','LineWidth',lineWidth);
%     plot(t,direc,'LineWidth',lineWidth);
    
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
    axes('Position',[m2(1) m2(2) x y]);
    plot(t,diffp,'LineWidth',lineWidth);

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Amp.(a.u.)')
    grid on;
%   box off;  %关闭右面和上面的坐标轴
%   axis tight; % 设置坐标紧凑型
    title('(b)', 'Units', 'normalized','FontWeight','bold','color', 'black','Position', z ,'FontSize', 16);
    set(gca,'YTickLabel',1*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3== 2
        return;
    end
    %% 第3幅图
    axes('Position',[m3(1) m3(2) x y]);
    plot(t,p,'LineWidth',lineWidth);
    hold on;
    plot(t,temp1,'LineWidth',lineWidth);
    scatter(loc_p./fs,top_p,'p','filled','r','LineWidth',lineWidth);
%     plot(t,en_median,'LineWidth',lineWidth);
%     hold on;
%     plot(t,direc/6,'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Amp.(a.u.)')
    grid on;
%   axis tight; % 设置坐标紧凑型
    title('(c)', 'Units', 'normalized','FontWeight','bold','color', 'black', 'Position', z ,'FontSize', 16);
    set(gca,'YTickLabel',1*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 3
        return;
    end
    %% 第4幅图
    axes('Position',[0.55,0.579,0.382,0.2028]);  % 设置图像位置
    si = 200;
    len = length(F);
    TF = abs(TF);
    mesh(T,F(len/2-si:len/2+si,:),TF(len/2-si:len/2+si,:));
    view(0,90); % 设置初始视角为俯视角
%     plot([1,2,3],'LineWidth',lineWidth);
%     colormap(jet);
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
  axis tight; % 设置坐标紧凑型
    title('(d)', 'Units', 'normalized', 'FontWeight','bold','color', 'white','Position',[0.054090837723508,0.753477954272838,0] ,'FontSize', 16);
    set(gca,'YTickLabel',10^-3*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 4
        return;
    end
     %% 第5幅图
    axes('Position',[m5(1) m5(2) x y]);
    
    plot(t,temp2,'LineWidth',lineWidth);
   axis([0 0.025 -1.2 1.2]);  % 设置具体横坐标和纵坐标的范围[x1,x2,y1,y2]
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Amp.(a.u.)')
    grid on;
%   axis tight; % 设置坐标紧凑型
    title('(f)', 'Units', 'normalized','FontWeight','bold','color', 'black', 'Position',z ,'FontSize', 16);
    set(gca,'YTickLabel',1*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 5
        return;
    end
    %% 第6幅图
%     axes('Position',[m6(1) m6(2) x y]);
%     plot([1,2,3],'LineWidth',lineWidth);
% 
%     % 横纵坐标轴+图框的粗细、大小、线宽
%     set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
%     xlabel('Time(ms)')
%     ylabel('Freq.(kHz)')
%     grid on;
% %   axis tight; % 设置坐标紧凑型
%     title('(f)', 'Units', 'normalized','FontWeight','bold','color', 'black', 'Position', z ,'FontSize', 16);
%     set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
%     if nargin - 3  == 6
%         return;
%     end
    %% 第7幅图
    axes('Position',[m7(1) m7(2) x y]);
    [Lt] = MOVE_API_STANDARD(200, 100, 125000, 3000, 70, 0.4);
    plot(t,Lt,'LineWidth',lineWidth);
    hold on;
    plot(t,Lt_reconstruct,'color','#EDB120','LineStyle','--','LineWidth',lineWidth);
    hold on;
    legend("原信号","重构信号");
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Lt(um)')
    grid on;
%   axis tight; % 设置坐标紧凑型
    title('(g)', 'Units', 'normalized','FontWeight','bold','color', 'black', 'Position', z ,'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 7
        return;
    end
    %% 第8幅图
    axes('Position',[0.55,0.280,0.382,0.2028]);

    len = length(F);
    TF_curb = abs(TF_curb);
    mesh(T,F(len/2-si:len/2+si,:),TF_curb(len/2-si:len/2+si,:));
    view(0,90); % 设置初始视角为俯视角
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; % 设置坐标紧凑型
    title('(e)', 'Units', 'normalized', 'FontWeight','bold','color', 'white', 'Position', [0.054090837723508,0.753477954272838,0] ,'FontSize', 16);
    set(gca,'YTickLabel',10^-3*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 8
        return;
    end
    %% 第9幅图
    axes('Position',[m9(1) m9(2) 0.873 0.1146]);
    plot(t,Lt-Lt_reconstruct,'LineWidth',lineWidth);
    hold on;
    plot(t,ones(1,N)*34.54*10^-9,'LineWidth',lineWidth);
    
    axis([0 0.025 -100*10^-9 50*10^-9]); 
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Error(nm)')
    legend("","RMSE")
    grid on;
%   axis tight; % 设置坐标紧凑型
    title('(h)', 'Units', 'normalized', 'FontWeight','bold','color', 'black', 'Position', [0.020683265340434,0.612063812858697,0] ,'FontSize', 16);
    set(gca,'YTickLabel',10^9*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 8
        return;
    end
    %% 第10幅图
%     axes('Position',[m10(1) m10(2) x y]);
%     plot([1,2,3],'LineWidth',lineWidth);
% 
%     % 横纵坐标轴+图框的粗细、大小、线宽
%     set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
%     xlabel('Time(ms)')
%     ylabel('Freq.(kHz)')
%     grid on;
% %   axis tight; % 设置坐标紧凑型
%     title('(j)', 'Units', 'normalized', 'FontWeight','bold','color', 'black', 'Position', z ,'FontSize', 16);
%     set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 8
        return;
    end
end

