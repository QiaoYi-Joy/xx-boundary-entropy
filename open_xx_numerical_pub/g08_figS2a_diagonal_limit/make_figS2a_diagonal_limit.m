function data=make_figS2a_diagonal_limit(profile,recompute)
%MAKE_FIGS2A_DIAGONAL_LIMIT 输入 quick/paper、是否重算；保存并画 Fig. S2(a)。
% false 从本组数据目录读取；绘图仅使用 plotData。
if nargin<1, profile='paper'; end
if nargin<2, recompute=true; end
profile=validatestring(profile,{'quick','paper'}); nB=80;
if strcmp(profile,'quick'), nB=30; end
here=fileparts(mfilename('fullpath')); root=fileparts(here);
addpath(fullfile(root,'common','01_physics'),fullfile(root,'common','02_numerics'));
cfg=struct('profile',profile,'b',linspace(.6,7,nB),'T',[.1 .25 .5]);
stem=['figS2a_diagonal_limit_' profile]; folder=fullfile(here,'data');
if ~exist(folder,'dir'), mkdir(folder); end
file=fullfile(folder,[stem '.mat']);
if recompute
    data=calc_figS2a_diagonal_limit(cfg);
else
    saved=load(file,'data'); data=saved.data;
end
assert(all(isfinite(data.plotData.deltaF),'all') && data.checkData.maxAbsAtLargestB<2e-4);
assert(all(diff(data.plotData.groundEnergy)>0));
if recompute, save(file,'data','-v7'); end
draw(data.plotData,fullfile(here,'figures'),stem);
fprintf('g08 [%s]: largest-b magnitude %.3e\n',profile,data.checkData.maxAbsAtLargestB);
end

%% 第二类：辅助绘图函数
function draw(d,folder,stem)
%DRAW 输入 plotData 和输出位置；绘制对角极限下自由能差，输出 PDF、PNG。
fig=figure('Visible','off','Color','w','Units','inches','Position',[1 1 4.3 3.2]);
ax=axes(fig); hold(ax,'on');
plot(ax,d.b,d.groundEnergy,'k--','LineWidth',1.15,'DisplayName','T=0');
for j=1:numel(d.parameters.temperature)               % 每个固定温度画一条自由能差曲线。
    plot(ax,d.b,d.deltaF(:,j),'LineWidth',1.2,'DisplayName',sprintf('T=%.2f',d.parameters.temperature(j)));
end
xlabel(ax,'Boundary parameter b'); ylabel(ax,'Boundary free-energy difference');
title(ax,'(a) Diagonal limit','FontWeight','normal'); legend(ax,'Location','best','Box','off');
xx_style_axes(ax); xx_export_figure(fig,folder,stem);
end
