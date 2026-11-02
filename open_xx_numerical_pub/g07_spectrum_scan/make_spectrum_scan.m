function data=make_spectrum_scan(profile,recompute)
%MAKE_SPECTRUM_SCAN 输入 quick/paper 和重算开关；保存数据并绘制 v13r6 Fig.S1 和 Fig.S2(d)。
if nargin<1, profile='paper'; end
if nargin<2, recompute=true; end
profile=validatestring(profile,{'quick','paper'});
here=fileparts(mfilename('fullpath')); root=fileparts(here);
addpath(fullfile(root,'common','01_physics'),fullfile(root,'common','02_numerics'));
out=fullfile(here,'data'); figs=fullfile(here,'figures'); if ~exist(out,'dir'), mkdir(out); end
file=fullfile(out,sprintf('spectrum_scan_%s.mat',profile));
cfg=struct('profile',profile,'parameters',struct('J',1,'aLeft',1.3,'aRight',.2, ...
    'bLeft',.7,'bRight',.6,'phiLeft',.1,'phiRight',.6), ...
    'classes',{{'ND-ND','ND-D','D-ND','D-D','ND-F','F-ND','D-F','F-D','F-F'}}, ...
    'Nlist',[4 6 8],'rootOptions',struct);
if strcmp(profile,'quick'), cfg.Nlist=[4 6]; end
if recompute
    data=calc_spectrum_scan(cfg);
else
    stored=load(file,'data'); data=stored.data;           % 重画版只读本包数据，不依赖旧目录。
end
assert(all(data.checkData.stateCount==2.^data.plotData.Nlist(:),'all'));
assert(max(data.plotData.spectrumError,[],'all')<1e-10,'Spectrum accuracy failed.');
assert(data.checkData.maxStringNewtonCorrection<2e-13,'String refinement failed.');
if recompute, save(file,'data','-v7'); end              % 核验失败不写入新数据。
d=data.plotData;                                        % 以下所有图只读取 plotData。

% Fig.S1(a,b)：同一次求根/ED 扫描的代表能谱和九类最大误差，沿用 r6 的科学量标题。
fig=figure('Visible','off','Color','w','Units','inches','Position',[1 1 7.05 3.05]);
p=d.panel; ax=axes(fig,'Position',[.080 .23 .37 .65]); hold(ax,'on');
plot(ax,p.spinSpectrum,p.spinSpectrum,'-','Color',[.35 .35 .35],'LineWidth',.65);
plot(ax,p.zeroSpectrum,p.spinSpectrum,'o','Color',[.478 .243 .616],'MarkerSize',2.8,'LineWidth',.55);
xlabel(ax,'Projected zero spectrum $E/J$','Interpreter','latex');
ylabel(ax,'Spin diagonalization $E/J$','Interpreter','latex');
title(ax,sprintf('(a) ND-ND, N=%d',p.N),'FontWeight','normal'); xx_style_axes(ax); ax.FontSize=8;
text(ax,.055,.93,sprintf('$2^{%d}/2=2^{%d}$ levels',p.N+1,p.N), ...
    'Units','normalized','Interpreter','latex','FontSize',8);
inset=axes(fig,'Position',[.285 .335 .135 .185]);
levelError=abs(p.zeroSpectrum-p.spinSpectrum);
semilogy(inset,0:numel(levelError)-1,max(levelError,1e-16),'Color',[.478 .243 .616],'LineWidth',.6);
set(inset,'YLim',[1e-16 1e-11],'YTick',[1e-15 1e-13],'XTick',[0 floor(numel(levelError)/2) numel(levelError)-1]);
xlabel(inset,'level'); ylabel(inset,'$|\Delta E|/J$','Interpreter','latex');
xx_style_axes(inset); inset.FontSize=6;
ax=axes(fig,'Position',[.555 .23 .34 .65]);
heatmapErrors(ax,d.spectrumError,d.classes,d.Nlist,[-15 -12]);
title(ax,'(b) Maximum spectral error','FontWeight','normal');
cb=colorbar(ax); cb.Label.Interpreter='latex';
cb.Label.String='$\log_{10}(\max|E^{\rm zeros}-E^{\rm ED}|/J)$'; cb.FontSize=7;
ax.Position=[.555 .23 .345 .65]; cb.Position=[.913 .23 .017 .65]; cb.Ticks=-15:-12;
xx_export_figure(fig,figs,sprintf('figS1_spectrum_validation_%s',profile));

% Fig.S2(d)：同一组直线根的原方程残差，不包含已单独检验的离轴边界弦。
fig=figure('Visible','off','Color','w','Units','inches','Position',[1 1 4.6 3.25]);
ax=axes(fig,'Position',[.12 .26 .72 .62]);
heatmapErrors(ax,d.lineRootResidual,d.classes,d.Nlist,[-16 -9]);
title(ax,'(d) Line-root residuals','FontWeight','normal');
cb=colorbar(ax); cb.Label.String='log_{10} residual';
xx_export_figure(fig,figs,sprintf('figS2d_line_root_residuals_%s',profile));
end

%% 第二类：辅助数值计算与作图
function heatmapErrors(ax,values,classes,sizes,limits)
%HEATMAPERRORS 输入最终误差矩阵，画十进对数热图；显示截断不修改数据。
v=log10(max(values,1e-16)); imagesc(ax,v); clim(ax,limits); colormap(ax,parula);
set(ax,'XTick',1:numel(classes),'XTickLabel',classes,'YTick',1:numel(sizes), ...
    'YTickLabel',compose('N=%d',sizes),'YDir','reverse');
xtickangle(ax,45); xlabel(ax,'Boundary class'); xx_style_axes(ax); ax.FontSize=8;
for i=1:size(v,1)                                     % 每行对应一个 N，逐格标注误差量级。
    for j=1:size(v,2)                                 % 每列对应一个有序边界组合。
        color=[0 0 0]; if v(i,j)<limits(1)+.25*diff(limits), color=[1 1 1]; end
        text(ax,j,i,sprintf('%.1f',v(i,j)),'HorizontalAlignment','center', ...
            'FontSize',5.5,'FontName','Times New Roman','Color',color);
    end
end
end
