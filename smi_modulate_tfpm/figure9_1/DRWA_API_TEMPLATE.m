function DRWA_API_TEMPLATE(varargin)
    % 获取当前文件夹下的所有文件
    files = dir('*.mat');
    for i = 1:length(files)
        % 加载所有.mat文件
        load(files(i).name);
    end
    
     % 根据需要设置最舒服的图窗大小
     set(gcf,'Position',[278,217,963,537]);
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
    axes('Position',[m1(1) m1(2) 0.382 0.618991683991683]);  % 设置图像位置 
    mesh(T,F,abs(TF));
    colormap Jet;
    view(0,90);

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Freq.(kHz)')
    grid on;
%   box off;  %关闭右面和上面的坐标轴
    axis tight; % 设置坐标紧凑型
%   axis([0 0.08 -20*10^-9 20*10^-9]);  % 设置具体横坐标和纵坐标的范围[x1,x2,y1,y2]
    title('(b)', 'Units', 'normalized', 'FontWeight','bold', 'color', 'white', 'Position', [0.057322915448572,0.900722850013309,0] ,'FontSize', 16);
    set(gca,'YTickLabel',10^-3*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 1
        return;
    end
    %% 第2幅图
    axes('Position',[m2(1) m2(2) x y]);
    
    p1 = -1 + (p1-min(p1))./(max(p1)-min(p1))*2;
    p11 = -1 + (p11-min(p11))./(max(p11)-min(p11))*2;
    
    % 计算信噪比
    noise = p_ini2 - p_ini22';
    ps = var(p_ini22)/length(p_ini22');
    pn = var(noise)/length(noise);
    SNR = 10*log10(ps/pn);
    disp(SNR);
    
    plot(t(2000:6000),p11(2000:6000),'color','#D95319','LineStyle','-');
    hold on;
    plot(t(2000:6000),p1(2000:6000),'color','#0072BD','LineStyle','-');
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Amp.(a.u)')
    grid on;
%     box off;  %关闭右面和上面的坐标轴
    axis tight; 
    title('(d)', 'Units', 'normalized','FontWeight','bold','color', 'black','Position', z ,'FontSize', 16);
    set(gca,'YTickLabel',1*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3== 2
        return;
    end
    %% 第3幅图
    axes('Position',[m3(1) m3(2) x y]);

    p2 = -1 + (p2-min(p2))/(max(p2)-min(p2))*2;
    p22 = -1 + (p22-min(p22))/(max(p22)-min(p22))*2;
    plot(t,p22,'color','#D95319','LineStyle','-');
    hold on;
    plot(t,p2,'color','#0072BD','LineStyle','-');
    legend("before TFIP","after TFIP")

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Amp.(a.u)')
    grid on;
%     axis tight; 
    title('(e)', 'Units', 'normalized','FontWeight','bold','color', 'black', 'Position', z ,'FontSize', 16);
    set(gca,'YTickLabel',1*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 3
        return;
    end
    %% 第4幅图
    axes('Position',[m4(1) m4(2) x y]);  % 设置图像位置
    plot(t,Ltt,'LineWidth',lineWidth,'color','red');
    hold on;
    plot(t,Lt_reconstruct,'color','black','LineStyle','--','LineWidth',lineWidth);

%     plot(45.4*10^-6*ones(1,length(Ltt)));
    legend("RMS=45.4nm");

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Amp.(a.u)')
    grid on;
%     axis tight; 
    title('(f)', 'Units', 'normalized', 'FontWeight','bold','color', 'black','Position',z ,'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 4
        return;
    end
    %% 第5幅图
%     axes('Position',[m5(1) m5(2) x y]);
%     plot(t,p2);
% 
%     % 横纵坐标轴+图框的粗细、大小、线宽
%     set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
%     xlabel('Time(ms)')
%     ylabel('Freq.(kHz)')
%     grid on;
%     axis tight; 
%     title('(e)', 'Units', 'normalized','FontWeight','bold','color', 'black', 'Position',z ,'FontSize', 16);
%     set(gca,'YTickLabel',1*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
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
%     if nargin - 3  == 6
%         return;
%     end
    %% 第7幅图
    axes('Position',[m7(1) m8(2) 0.382 y]);
    plot(t,p_ini2,'LineWidth',lineWidth);

    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
%     xlabel('Time(ms)')
    ylabel('Amp.(a.u)')
    grid on;
    axis tight; 
    title('(a)', 'Units', 'normalized','FontWeight','bold','color', 'black', 'Position', z ,'FontSize', 16);
    set(gca,'YTickLabel',1*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 7
        return;
    end
    %% 第8幅图
    axes('Position',[m8(1) m8(2) x y]);
    
    plot(t,p11,'color','#D95319','LineStyle','-');
    hold on;
    plot(t,p1,'color','#0072BD','LineStyle','-');
    legend("before TFIP","after TFIP")

%     legend("before TFIP","after harmonic")
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
%     xlabel('Time(ms)')
    ylabel('Amp.(a.u)')
    grid on;
%     box off;
%     axis tight; 
    title('(c)', 'Units', 'normalized', 'FontWeight','bold','color', 'black', 'Position', z ,'FontSize', 16);
    set(gca,'YTickLabel',1*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 8
        return;
    end
end

