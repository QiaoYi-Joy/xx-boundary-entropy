function data=make_fig2a_boundary_classes(profile,recompute)
%MAKE_FIG2A_BOUNDARY_CLASSES 输入 quick/paper 和是否重算；保存数据并绘制 Fig. 2(a)。
% 重算默认 true；false 读取本组已有同配置数据。绘图只读取 plotData。
if nargin<1, profile='paper'; end
if nargin<2, recompute=true; end
profile=validatestring(profile,{'quick','paper'}); nT=64;
if strcmp(profile,'quick'), nT=24; end
here=fileparts(mfilename('fullpath')); root=fileparts(here);
addpath(fullfile(root,'common','01_physics'),fullfile(root,'common','02_numerics'));
cfg=struct('profile',profile,'aLeft',1.3,'aRight',.2,'bLeft',.7,'bRight',.6, ...
    'temperature',logspace(log10(.01),log10(.4),nT).');
stem=['fig2a_boundary_classes_' profile]; folder=fullfile(here,'data');
if ~exist(folder,'dir'), mkdir(folder); end
file=fullfile(folder,[stem '.mat']);
if recompute
    data=calc_fig2a_boundary_classes(cfg);
else
    saved=load(file,'data'); data=saved.data;
end
assert(all(isfinite(data.plotData.boundaryRatio),'all'));
assert(data.checkData.maxPhiConvergenceDifference<1e-9,'Phi 尾项收敛检查未通过。');
if recompute, save(file,'data','-v7'); end
draw(data.plotData,fullfile(here,'figures'),stem);
fprintf('g04 [%s]: Phi convergence %.3e, finite-T offset %.3e\n', ...
    profile,data.checkData.maxPhiConvergenceDifference,max(abs(data.checkData.lowestTemperatureOffset)));
end

%% 第二类：辅助绘图函数
function draw(d,folder,stem)
%DRAW 输入 plotData 和输出位置；绘制三类边界比值，输出 PDF、PNG。
colors=[.15 .36 .72;.84 .33 .19;.25 .58 .36];
fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 10.8 7.2]);
ax=axes(fig); hold(ax,'on');
for j=1:3                                             % 每类一条曲线及对应的低温极限参考线。
    semilogx(ax,d.temperature,d.boundaryRatio(:,j),'LineWidth',1.45, ...
        'Color',colors(j,:),'DisplayName',d.classes{j});
    yline(ax,d.targets(j),':','Color',colors(j,:),'HandleVisibility','off');
end
xlabel(ax,'$T/J$','Interpreter','latex'); ylabel(ax,'$(F_{\rm b}-E_{\rm b})/T$','Interpreter','latex');
legend(ax,'Units','normalized','Position',[.65 .28 .24 .18],'Box','off'); % 图例放在两条曲线之间的空白区。
xlim(ax,d.temperature([1 end])); ylim(ax,[-.08 .82]);
set(ax,'XScale','log'); xx_style_axes(ax); xx_export_figure(fig,folder,stem);
end
