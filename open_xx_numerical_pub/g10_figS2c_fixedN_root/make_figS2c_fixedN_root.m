function data=make_figS2c_fixedN_root(profile,recompute)
%MAKE_FIGS2C_FIXEDN_ROOT 输入 quick/paper 和重算开关；计算、存储或只重画 Fig.S2(c)。
if nargin<1, profile='paper'; end
if nargin<2, recompute=true; end
profile=validatestring(profile,{'quick','paper'});
here=fileparts(mfilename('fullpath')); root=fileparts(here);
addpath(fullfile(root,'common','01_physics'),fullfile(root,'common','02_numerics'));
out=fullfile(here,'data'); figs=fullfile(here,'figures'); if ~exist(out,'dir'), mkdir(out); end
file=fullfile(out,sprintf('figS2c_fixedN_root_%s.mat',profile));
cfg=struct('profile',profile,'parameters',struct('N',10,'J',1,'aLeft',1.3,'aRight',.2, ...
    'bLeft',.7,'phiLeft',.1,'phiRight',.6),'b',2:8);
nGrid=4097; if strcmp(profile,'quick'), cfg.b=2:4; nGrid=1025; end
cfg.rootOptions=struct('xMin',1e-10,'xMax',8,'nGrid',nGrid,'tolEndpointX',1e-7, ...
    'acceptResidual',1e-8,'acceptPhaseCycles',1e-9);
if recompute
    data=calc_figS2c_fixedN_root(cfg);
else
    stored=load(file,'data'); data=stored.data;
end
assert(data.checkData.maxLogResidual<2e-7 && all(data.checkData.modeCount==data.plotData.parameters.N+1));
d=data.plotData;                                       % 以下只读最终作图数据。
assert(all(diff(d.xMax)>0));
if recompute, save(file,'data','-v7'); end             % 所有检查通过后才保存。
fig=figure('Visible','off','Color','w','Units','inches','Position',[1 1 4.2 3.2]);
ax=axes(fig); hold(ax,'on');
plot(ax,d.b,d.offset,'-o','LineWidth',1.35,'MarkerSize',4);
yline(ax,d.lastPointGuide,':','Color',[.35 .35 .35]);    % 最大 b 的末点只作为渐近参考线。
xlabel(ax,'Regularizing parameter b'); ylabel(ax,'$x_{\max}-b/2$','Interpreter','latex');
title(ax,sprintf('(c) Fixed N=%d farthest root',d.parameters.N),'FontWeight','normal');
xx_style_axes(ax); xx_export_figure(fig,figs,sprintf('figS2c_fixedN_root_%s',profile));
end
