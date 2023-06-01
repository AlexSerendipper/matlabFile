function DRWA_API_TEMPLATE_fig3(varargin)
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
    load("D:\matlab save\self-mixing\smi_MODULATE\fig3_relate.mat");
    axes('Position',[m1(1) m1(2) x y]);  axis tight;  % 设置图像位置 + 坐标紧凑型

    plot(t(100:13400),p(100:13400),'LineWidth',lineWidth);

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
	ylabel('Amp.(a.u)')
    grid on;
    axis tight; 
    title('(a)', 'Units', 'normalized', 'Position', z ,'FontSize', 16);
%     set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 1
        return;
    end
    %% 第2幅图
    axes('Position',[m2(1) m2(2) x y]);  axis tight; 
     plot(t(100:100+13301),amp1(18000:18000+13301),'LineWidth',lineWidth)
%     plot(t(100:13400),amp1(100:13400),'LineWidth',lineWidth)
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
	ylabel('Amp.(a.u)')
    grid on;
    axis tight; 
    title('(b)', 'Units', 'normalized', 'Position', z ,'FontSize', 16);
    set(gca,'YTickLabel',10^-3*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3== 2
        return;
    end
    %% 第3幅图
    axes('Position',[m3(1) m3(2) x y]);  axis tight;
    plot(t(100:13400),p1(100:13400),'LineWidth',lineWidth);

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Amp.(a.u) ')
    grid on;
    axis tight; 
    title('(c)', 'Units', 'normalized', 'Position', z ,'FontSize', 16);
%     set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 3
        return;
    end
    %% 第4幅图
    axes('Position',[m4(1) m4(2) x y]);  axis tight;  % 设置图像位置 + 坐标紧凑型 
    plot(t(100:13400),p2(100:13400),'LineWidth',lineWidth);

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Amp.(a.u) ')
    grid on;
    axis tight; 
    title('(d)', 'Units', 'normalized', 'Position',z ,'FontSize', 16);
%     set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 4
        return;
    end
    %% 第5幅图
    axes('Position',[m5(1) m5(2) x y]);  axis tight;  % 设置图像位置 + 坐标紧凑型 
    load("D:\matlab save\self-mixing\smi_MODULATE\temp_experi._phiF_wrapped_times3.mat")
    plot(t(100:13400),phiF_wrapped_init(100:13400),'LineWidth',lineWidth);

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Amp.(a.u) ')
    grid on;
    axis tight; 
    title('(e)', 'Units', 'normalized', 'Position',z ,'FontSize', 16);
%     set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 5
        return;
    end
    %% 第6幅图
    axes('Position',[m6(1) m6(2) x y]);  axis tight;
    load("D:\matlab save\self-mixing\smi_MODULATE\temp_experi._phiF_wrapped_times4.mat")
    plot(t(100:13400),phiF_wrapped_init(100:13400),'LineWidth',lineWidth);

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Amp.(a.u) ')
    grid on;
    axis tight; 
    title('(f)', 'Units', 'normalized', 'Position', z ,'FontSize', 16);
%     set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 6
        return;
    end
    %% 第7幅图
    axes('Position',[m7(1) m7(2) x y]);  axis tight;
    plot(t(100:13400),Lt_reconstruct1(100:13400),'LineWidth',lineWidth);
    hold on;
    load("D:\matlab save\self-mixing\smi_MODULATE\temp_experi._recon._times3.mat");
    plot(t(100:13400),Lt_reconstruct_times3(100:13400),'r','LineWidth',lineWidth)
    hold on;
    load("D:\matlab save\self-mixing\smi_MODULATE\temp_experi._recon._times4.mat");
    plot(t(100:13400),Lt_reconstruct_times4(100:13400),'g','LineWidth',lineWidth);
    legend("reference","OSPM(N=3)","OSPM(N=4)");

    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Displacement(um)')
    grid on;
    axis tight; 
    title('(g)', 'Units', 'normalized', 'Position', z ,'FontSize', 16);
    set(gca,'YTickLabel',10^6*get(gca,'YTick'))  % 乘10^6，让图像的纵坐标显示的是微米量级
    if nargin - 3  == 7
        return;
    end
    %% 第8幅图
    axes('Position',[m8(1) m8(2) x y]);  axis tight;
    load("D:\matlab save\self-mixing\smi_MODULATE\temp_experi._recon._times1.mat");
    plot(t(100:13400),Lt_reconstruct1(100:13400)-Lt_reconstruct_times3(100:13400),'LineWidth',lineWidth); grid on;  % 作图 + 网格
    hold on;
    plot(t(100:13400),(Lt_reconstruct1(100:13400)-Lt_reconstruct_times4(100:13400)),'LineWidth',lineWidth);
    
    % 横纵坐标轴+图框的粗细、大小、线宽
    set(gca,'LineWidth',lineWidth,'FontName',fontName,'FontWeight',fontWeight,'FontSize',fontSize);
    xlabel('Time(ms)')
    ylabel('Error(nm)');
    grid on;
    axis tight; 
    title('(h)', 'Units', 'normalized', 'Position', z ,'FontSize', 16);
    set(gca,'YTickLabel',10^9*get(gca,'YTick'))  % 乘10^9，让图像的纵坐标显示的是微米量级
    legend("OSPM(N=3)","OSPM(N=4)");
    if nargin - 3  == 8
        return;
    end
end

