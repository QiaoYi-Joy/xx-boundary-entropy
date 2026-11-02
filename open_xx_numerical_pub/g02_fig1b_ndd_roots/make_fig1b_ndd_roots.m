function data=make_fig1b_ndd_roots(profile,recompute)
%MAKE_fig1b_ndd_roots 输入 quick/paper 和是否重算；读取或保存三层数据，绘制 fig1b_ndd_roots。
if nargin<1, profile='paper'; end
if nargin<2, recompute=true; end
profile=validatestring(profile,{'quick','paper'});
here=fileparts(mfilename('fullpath')); root=fileparts(here);
addpath(fullfile(root,'common','01_physics'),fullfile(root,'common','02_numerics'));
out=fullfile(here,'data'); figs=fullfile(here,'figures'); if ~exist(out,'dir'), mkdir(out); end
file=fullfile(out,sprintf('fig1b_ndd_roots_%s.mat',profile));
p=struct('N',8,'J',1,'aLeft',1.3,'aRight',1.1,'bLeft',.7,'bRight',4,'phiLeft',.1,'phiRight',.6);
cfg=struct('profile',profile,'parameters',p);
nGrid=4097; if strcmp(profile,'quick'), nGrid=513; end
cfg.rootOptions=struct('xMin',1e-8,'xMax',8,'nGrid',nGrid,'acceptResidual',1e-10,'acceptPhaseCycles',1e-10);
if recompute
    data=calc_fig1b_ndd_roots(cfg);
else
    stored=load(file,'data'); data=stored.data;       % 只重画时无需重新求根。
end
M=data.plotData.parameters.N;
assert(data.checkData.modeCount==M && data.checkData.positionCount==4*M);
assert(data.checkData.maxLogResidual<2e-10);
if recompute, save(file,'data','-v7'); end            % 全部核验通过后才覆盖本组数据。
d=data.plotData;                                    % 以下画图只读取最终作图数据。

blue=[0 .447 .741]; gray=[.52 .52 .52];
fig=figure('Visible','off','Color','w','Units','inches','Position',[1 1 4.05 3.55]);
ax=axes(fig,'Position',[.16 .16 .79 .76]); hold(ax,'on');
yline(ax,[pi/4 -pi/4],':','Color',gray,'HandleVisibility','off');
xline(ax,0,'--','Color',[.72 .72 .72],'HandleVisibility','off');
z=d.allowedPositions(:);
scatter(ax,real(z),imag(z),24,'filled','MarkerFaceColor',blue, ...
    'MarkerEdgeColor','w','DisplayName','allowed zero positions'); % 每点均为允许位置，并非全被同一态占据
xmax=max(abs(real(z)))+.38; xlim(ax,[-xmax xmax]); ylim(ax,[-pi/2 pi/2]);
yticks(ax,[-pi/2 -pi/4 0 pi/4 pi/2]);
yticklabels(ax,{'$-\pi/2$','$-\pi/4$','$0$','$\pi/4$','$\pi/2$'});
xlabel(ax,'$\mathrm{Re}\,z$','Interpreter','latex'); ylabel(ax,'$\mathrm{Im}\,z$','Interpreter','latex');
title(ax,sprintf('strict ND--D: $M=N=%d$ allowed quartets',d.parameters.N),'Interpreter','latex','FontWeight','normal');
legend(ax,'Location','southwest','Box','off','Interpreter','latex');
set(ax,'TickLabelInterpreter','latex'); axis(ax,'square'); % 字体和朝内刻度统一交给公共函数。
xx_style_axes(ax); xx_export_figure(fig,figs,sprintf('fig1b_ndd_roots_%s',profile));
end
