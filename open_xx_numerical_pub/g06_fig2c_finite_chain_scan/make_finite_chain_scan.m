function data=make_finite_chain_scan(profile,recompute)
%MAKE_FINITE_CHAIN_SCAN 计算或读取共用扫描，核验并保存，再分别绘制 Fig. 2(c)、S3(a,b)。
% 输入 profile='quick'/'paper'、recompute=true/false；返回 meta/plotData/checkData。
% false 只读取本组已保存数据；画图函数仅接受 plotData，不读取 checkData。
if nargin<1, profile='quick'; end
if nargin<2, recompute=true; end
profile=validatestring(profile,{'quick','paper'});
here=fileparts(mfilename('fullpath')); root=fileparts(here);
addpath(fullfile(root,'common','01_physics'),fullfile(root,'common','02_numerics'));
folder=fullfile(here,'data'); figures=fullfile(here,'figures');
file=fullfile(folder,['finite_chain_scan_' profile '.mat']);
if recompute
    if strcmp(profile,'paper')
        Nlist=[16 32 64 128 256 512 1024 2048 4096];
        displayN=[128 512 1024 2048 4096]; nT=240;
    else
        Nlist=[16 32 64]; displayN=Nlist; nT=80;
    end
    rootOptions=struct('xMin',1e-10,'xMax',8,'tolDuplicateX',1e-11, ...
        'tolRefineX',1e-14,'acceptResidual',2e-7,'acceptPhaseCycles',1e-9);
    cfg=struct('profile',profile,'J',1,'Nlist',Nlist,'displayN',displayN, ...
        'temperature',logspace(-5,log10(.6),nT).','fixedTemperature',[.025 .05 .10].', ...
        'aLeft',1.3,'aRight',1.1,'bLeft',.7,'phiLeft',.1,'phiRight',.6, ...
        'target',log(2)/2,'rootOptions',rootOptions);
    data=calc_finite_chain_scan(cfg);                 % 所有图共用这一遍计算，不重复求根与对角化。
else
    loaded=load(file,'data'); data=loaded.data;
end
data.meta.panels={'Fig. 2(c)','Fig. S3(a)','Fig. S3(b)'}; % 读取旧 MAT 时仅在内存中统一新图号。
p=data.plotData; c=data.checkData;
assert(all(c.modeCountByN==repmat(p.Nlist,1,2),'all'),'模数检查失败。');
assert(max(c.maxLogResidualByN,[],'all')<2e-7,'对数求根残差超限。');
assert(max(p.modeEnergyDifference,[],'all')<5e-8,'普通零点与独立 Majorana 谱不一致。');
assert(max(c.maxRatioDifferenceByN)<2e-7,'独立谱所得 R_N 不一致。');
if recompute
    if ~exist(folder,'dir'), mkdir(folder); end
    save(file,'data','-v7');                          % 只保存最终绘图量和必要的少量核验量。
end
draw_finite_chains(p,figures,profile);
draw_supplementary(p,figures,profile);
fprintf('g05 [%s]: residual %.3e, mode difference %.3e, R_N difference %.3e\n', ...
    profile,max(c.maxLogResidualByN,[],'all'),max(p.modeEnergyDifference,[],'all'), ...
    max(c.maxRatioDifferenceByN));
end

%% 第二类：辅助数值计算函数（以下仅负责显示 plotData）
function draw_finite_chains(p,folder,profile)
%DRAW_FINITE_CHAINS 输入 plotData/输出目录/配置名，绘制 Fig. 2(c) 的有限链和精确热力学参照。
[fig,ax]=new_axes; q=p.fivePercentInterval;
if numel(q)==2
    patch(ax,q([1 2 2 1]),p.target*[.95 .95 1 1],[.87 .96 .88], ...
        'EdgeColor',[.36 .64 .40],'LineWidth',.55,'HandleVisibility','off');
    text(ax,sqrt(prod(q)),p.target+.018,'5%','HorizontalAlignment','center','FontSize',9);
end
colors=[.53 .53 .53;0 .447 .698;.835 .369 0;0 .620 .451;.478 .243 .616];
for j=1:numel(p.displayN)                            % 只画文章指定的链长，quick 则画三个测试链长。
    semilogx(ax,p.temperature,p.finiteRatio(:,j),'LineWidth',1.25,'Color',colors(j,:), ...
        'DisplayName',sprintf('$N=%d$',p.displayN(j)));
end
semilogx(ax,p.temperature,p.thermodynamicRatio,'k--','LineWidth',1.35,'DisplayName','$N\to\infty$');
yline(ax,p.target,':','HandleVisibility','off');
set(ax,'XScale','log','XLim',[1e-4,.65],'YLim',[-.012,.382]);
xlabel(ax,'$T$','Interpreter','latex'); ylabel(ax,'$R_N(T)$','Interpreter','latex');
title(ax,'(c) Finite chains');
lg=legend(ax,'Location','northwest','Box','off','Interpreter','latex','FontSize',8);
lg.ItemTokenSize=[16 12]; drawnow; lg.Position(2)=lg.Position(2)-.05; % 图例避开热力学虚线与上升曲线。
xx_style_axes(ax); xx_export_figure(fig,folder,['fig2c_finite_chains_' profile]);
end

function draw_supplementary(p,folder,profile)
%DRAW_SUPPLEMENTARY 输入 plotData/目录/配置名，输出合并 Fig. S3；沿用 v13r6 两个物理标题。
fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 23 8]);
layout=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(layout,1); hold(ax,'on');
colors=[0 .447 .698;.835 .369 0;0 .620 .451]; markers={'o','s','^'};
for j=1:numel(p.fixedTemperature)                   % 三条曲线分别对应三个固定温度。
    loglog(ax,p.Nlist,abs(p.finiteSizeDifference(:,j)),'LineWidth',1.25, ...
        'Color',colors(j,:),'Marker',markers{j},'MarkerFaceColor','w', ...
        'DisplayName',sprintf('$T=%.3g$',p.fixedTemperature(j)));
end
set(ax,'XScale','log','YScale','log');
xlabel(ax,'$N$','Interpreter','latex'); ylabel(ax,'$|R_N(T)-R(b,T)|$','Interpreter','latex');
title(ax,'(a) Fixed-temperature convergence'); legend(ax,'Location','southwest','Box','off','Interpreter','latex');
xx_style_axes(ax);
ax=nexttile(layout,2); hold(ax,'on'); labels={'ND-D','D-D'};
for j=1:2                                         % 分别画 ND-D 与 D-D 的最大逐模能差。
    loglog(ax,p.Nlist,p.modeEnergyDifference(:,j),'LineWidth',1.25,'Color',colors(j,:), ...
        'Marker',markers{j},'MarkerFaceColor','w','DisplayName',labels{j});
end
set(ax,'XScale','log','YScale','log');
xlabel(ax,'$N$','Interpreter','latex');
ylabel(ax,'$\max_j|\varepsilon_j^{\rm roots}-\varepsilon_j^{\rm Maj}|$','Interpreter','latex');
title(ax,'(b) Mode-energy differences'); legend(ax,'Location','south','Box','off');
xx_style_axes(ax); xx_export_figure(fig,folder,['figS3_finite_size_validation_' profile]);
end

function [fig,ax]=new_axes
%NEW_AXES 创建隐藏的白底画布；输出图与轴句柄，不读取任何计算数据。
fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 12 8]);
ax=axes(fig); hold(ax,'on');
end
