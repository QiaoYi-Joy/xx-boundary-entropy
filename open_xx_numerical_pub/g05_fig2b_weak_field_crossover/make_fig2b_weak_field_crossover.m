function data=make_fig2c_weak_field_crossover(profile,recompute)
%MAKE_FIG2C_WEAK_FIELD_CROSSOVER 输入 quick/paper、是否重算；保存数据并画新 Fig. 2(b)。
% 函数名、目录及 MAT 文件名沿用旧名，避免破坏历史数据路径；输出图片按新编号命名。
% false 读取本组已有数据；绘图函数仅使用 plotData。
if nargin<1, profile='paper'; end
if nargin<2, recompute=true; end
profile=validatestring(profile,{'quick','paper'}); nS=76;
if strcmp(profile,'quick'), nS=32; end
here=fileparts(mfilename('fullpath')); root=fileparts(here);
addpath(fullfile(root,'common','01_physics'),fullfile(root,'common','02_numerics'));
cfg=struct('profile',profile,'bValues',3:6,'target',log(2)/2, ...
    'scaledTemperature',logspace(-2,log10(5e2),nS).','absTol',2e-12,'relTol',2e-11);
stem=['fig2c_weak_field_crossover_' profile]; folder=fullfile(here,'data');
if ~exist(folder,'dir'), mkdir(folder); end
file=fullfile(folder,[stem '.mat']);
if recompute
    data=calc_fig2c_weak_field_crossover(cfg);
else
    saved=load(file,'data'); data=saved.data;
end
data.meta.panels={'Fig. 2(b)'}; % 读取旧 MAT 时仅在内存中统一新图号。
assert(all(isfinite(data.plotData.boundaryRatio),'all'));
assert(all(data.plotData.boundaryRatio>=-2e-10 & ...
    data.plotData.boundaryRatio<=data.plotData.target+2e-10,'all')); % 数值范围，不把特定弱场中点写成通用条件。
assert(data.checkData.maximumMonotonicIncrease<=2e-10);
if recompute, save(file,'data','-v7'); end
draw(data.plotData,fullfile(here,'figures'),['fig2b_weak_field_crossover_' profile]);
fprintf('g06 [%s]: T1/2 exp(2b)=%s\n',profile,mat2str(data.checkData.halfScaledTemperature,7));
end

%% 第二类：辅助绘图函数
function draw(d,folder,stem)
%DRAW 输入 plotData 和输出位置；绘制弱场缩放曲线，输出 PDF、PNG。
fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 10.8 7.2]);
ax=axes(fig); hold(ax,'on'); colors=lines(numel(d.parameters.bValues));
for j=1:numel(d.parameters.bValues)                    % 各 b 的曲线共用缩放横坐标。
    semilogx(ax,d.scaledTemperature,d.boundaryRatio(:,j),'LineWidth',1.4, ...
        'Color',colors(j,:),'DisplayName',sprintf('$b=%g$',d.parameters.bValues(j)));
end
yline(ax,d.target,':','HandleVisibility','off'); yline(ax,d.halfTarget,'--','HandleVisibility','off');
xlabel(ax,'$T e^{2b}/J$','Interpreter','latex'); ylabel(ax,'$R(b,T)$','Interpreter','latex');
legend(ax,'Interpreter','latex','Location','southwest','Box','off');
xlim(ax,d.scaledTemperature([1 end])); ylim(ax,[-.015 .38]);
set(ax,'XScale','log'); xx_style_axes(ax); xx_export_figure(fig,folder,stem);
end
