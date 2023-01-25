%% 画图函数，设置图窗，输入图像大小，输出无边框，无坐标的图形
%% x,y 为输入的图片大小
function DRAW_API_drawClear(x,y)  
    x = x * 1.70667 / 2.66667;
    y = y * 1.70667 / 2.66667;
    set(0,'defaultfigurecolor','w')   % 设置背景颜色为白色
    set(gcf,'Position',[0,0,x,y]);
    set(gca,'Position',[0,0,1,1]);  %去除白边（用来设置绘制的图像距离画板Figure边界的距离，ab代表绘图起始坐标，cd代表宽度高度。取值范围都为0-1.）
%     plot(z,'Linewidth', 1.5);
    axis tight;
    axis off;
    saveas(gcf,['D://桌面存放//',num2str(1),'.png'],'png');
end