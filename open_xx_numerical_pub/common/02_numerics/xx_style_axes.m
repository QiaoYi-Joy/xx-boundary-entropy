function xx_style_axes(ax)
%XX_STYLE_AXES 统一图像样式；输入坐标轴句柄，无数值输出，不改变数据和坐标尺度。
% 所有主、次刻度均朝内；图名采用物理量，不写算法调试措辞。
set(ax,'FontName','Times New Roman','FontSize',10,'LineWidth',.7, ...
    'TickDir','in','Box','on','XMinorTick','on','YMinorTick','on', ...
    'TickLabelInterpreter','latex');
ax.XLabel.FontSize=11; ax.YLabel.FontSize=11; ax.Title.FontSize=11;
ax.Title.FontWeight='normal';
end
