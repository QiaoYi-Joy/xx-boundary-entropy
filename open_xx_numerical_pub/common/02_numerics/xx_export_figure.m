function xx_export_figure(fig,folder,stem)
%XX_EXPORT_FIGURE 保存并关闭图像；输入图句柄、目录、无扩展名文件名，无返回值。
% 输出矢量 PDF 和 260 dpi PNG；只处理画图文件，不保存计算诊断。
if ~exist(folder,'dir'), mkdir(folder); end
drawnow;
exportgraphics(fig,fullfile(folder,[stem '.pdf']),'ContentType','vector');
exportgraphics(fig,fullfile(folder,[stem '.png']),'Resolution',260);
close(fig);
end
