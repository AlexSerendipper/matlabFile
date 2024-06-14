function DRWA_API_TEMPLATE(varargin)
    % 获取当前文件夹下的所有文件
    files = dir('*.mat');
    for i = 1:length(files)
        % 加载所有.mat文件
        load(files(i).name);
    end
    
     % 根据需要设置最舒服的图窗大小
      set(gcf,'Position',[285,50,1202,872]);
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
    axes('Position',[m1(1) 0.5488 0.401 0.385149425287356]);  % 设置图像位置 
    
plot(phi01,phiF,'LineWidth',lineWidth);
hold on;

plot(phi02,phiF,'LineWidth',lineWidth,'LineStyle',':');
hold on;

plot(phi03,phiF,'LineWidth',lineWidth,'LineStyle','-.');
legend('α=2','α=4','α=7');

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel("\Phi_0(rad)")
    ylabel('\Phi_F(rad)')
    grid on;
%   box off;  %关闭右面和上面的坐标轴
    axis tight; % 设置坐标紧凑型
%   axis([0 0.08 -20*10^-9 20*10^-9]);  % 设置具体横坐标和纵坐标的范围[x1,x2,y1,y2]
    title('(a)', 'Units', 'normalized', 'FontWeight','bold', 'color', 'black', 'Position', [0.061790825688073,0.875631901840491,0] ,'FontSize', 16);
    set(gca,'YTickLabel',1*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
%   set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')  % 不显示x轴和y轴
    if nargin - 3  == 1
        return;
    end
    %% 第2幅图
    axes('Position',[0.538 0.5488 0.401 0.385149425287356]);
    
plot(phi11,phiF,'LineWidth',lineWidth);
hold on;

plot(phi12,phiF,'LineWidth',lineWidth,'LineStyle',':');
hold on;

plot(phi13,phiF,'LineWidth',lineWidth,'LineStyle','-.');
legend('α=2','α=4','α=7');

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel("\Phi_0(rad)")
    ylabel('\Phi_F(rad)')
    grid on;
%     box off;  %关闭右面和上面的坐标轴
    axis tight; 
    title('(b)', 'Units', 'normalized','FontWeight','bold','color', 'black','Position', [0.061790825688073,0.875631901840491,0] ,'FontSize', 16);
    set(gca,'YTickLabel',1*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3== 2
        return;
    end
    %% 第3幅图
    axes('Position',[0.06 0.306537865818863 0.877544867193109 0.167457406048749]);
plot(t,p1,'LineWidth',lineWidth);
hold on;

plot(t,p2,'LineWidth',lineWidth,'LineStyle',':');
hold on;

plot(t,p3,'LineWidth',lineWidth,'LineStyle','-.');
legend('α=2','α=4','α=7');
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Amp.(a.u.)')
    grid on;
%     axis tight; 
    title('(c)', 'Units', 'normalized','FontWeight','bold','color', 'black', 'Position', [0.026844589128933,0.727631205673759,0] ,'FontSize', 16);
    set(gca,'YTickLabel',1*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 3
        return;
    end
    %% 第4幅图
%     axes('Position',[m4(1) m4(2) x y]);  % 设置图像位置
%     plot([1,2,3],'LineWidth',lineWidth);
% 
%     % 横纵坐标轴+图框的粗细、大小、线宽
%     set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
%     xlabel('Time(ms)')
%     ylabel('Freq.(kHz)')
%     grid on;
%     axis tight; 
%     title('(d)', 'Units', 'normalized', 'FontWeight','bold','color', 'black','Position',z ,'FontSize', 16);
%     set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
%     if nargin - 3  == 4
%         return;
%     end
    %% 第5幅图
    axes('Position',[0.058529411764706 0.082758620689655 0.877544867193109 0.167457406048749]);
plot(t,p4,'LineWidth',lineWidth);
hold on;

plot(t,p5,'LineWidth',lineWidth,'LineStyle',':');
hold on;

plot(t,p6,'LineWidth',lineWidth,'LineStyle','-.');
legend('α=2','α=4','α=7');

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Amp.(a.u.)')
    grid on;
%     axis tight; 
    title('(d)', 'Units', 'normalized','FontWeight','bold','color', 'black', 'Position',[0.026844589128933,0.727631205673759,0] ,'FontSize', 16);
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
%     axis tight; 
%     title('(f)', 'Units', 'normalized','FontWeight','bold','color', 'black', 'Position', z ,'FontSize', 16);
%     set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 6
        return;
    end
    %% 第7幅图
    axes('Position',[m7(1) m7(2) x y]);
    plot([1,2,3],'LineWidth',lineWidth);

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
    axis tight; 
    title('(g)', 'Units', 'normalized','FontWeight','bold','color', 'black', 'Position', z ,'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 7
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
    axis tight; 
    title('(h)', 'Units', 'normalized', 'FontWeight','bold','color', 'black', 'Position', z ,'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 8
        return;
    end
end

