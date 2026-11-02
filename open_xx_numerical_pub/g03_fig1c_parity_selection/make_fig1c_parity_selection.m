function data=make_fig1c_parity_selection(profile,recompute)
%MAKE_fig1c_parity_selection 输入 quick/paper 和是否重算；读取或保存三层数据，绘制 fig1c_parity_selection。
if nargin<1, profile='paper'; end
if nargin<2, recompute=true; end
profile=validatestring(profile,{'quick','paper'});
here=fileparts(mfilename('fullpath')); root=fileparts(here);
addpath(fullfile(root,'common','01_physics'),fullfile(root,'common','02_numerics'));
out=fullfile(here,'data'); figs=fullfile(here,'figures'); if ~exist(out,'dir'), mkdir(out); end
file=fullfile(out,sprintf('fig1c_parity_selection_%s.mat',profile));
p=struct('N',8,'J',1,'aLeft',1.3,'aRight',1.1,'bLeft',.7,'bRight',4,'phiLeft',.1,'phiRight',.6);
cfg=struct('profile',profile,'parameters',p);
cfg.exampleFlips={[],[2 p.N+1],2};                 % 随 N 改变最后一列，不写死第9个四元组。
if recompute
    data=calc_fig1c_parity_selection(cfg);
else
    stored=load(file,'data'); data=stored.data;       % 只重画时无需重新求根。
end
d=data.plotData;
assert(all(data.checkData.selectedPairCount==d.M));
assert(isequal(d.parity,prod(d.selectedSigns,2)) && isequal(d.isPhysical,d.parity==d.kappa));
assert(numel(d.rowLabels)==size(d.selectedSigns,1));
if recompute, save(file,'data','-v7'); end            % 先核验选择规则，再保存数据。
% 以下画图只读取最终作图数据。

blue=[0 .447 .741]; orange=[.85 .325 .098]; red=[.62 .10 .16];
gray=[.55 .55 .55]; light=[.80 .83 .86]; ink=[.12 .15 .18];
nRows=size(d.selectedSigns,1); rowY=(nRows-1:-1:0)*.94+.54; headerY=rowY(1)+.49;
statusX=d.M+1.85; offset=.18;
fig=figure('Visible','off','Color','w','Units','inches', ...
    'Position',[1 1 max(5,.46*d.M+2.96) max(2.1,.94*nRows+.33)]);
ax=axes(fig,'Position',[.035 .06 .94 .89]); hold(ax,'on');
set(ax,'TickDir','in','XMinorTick','on','YMinorTick','on'); % 即使隐藏坐标轴也统一向内刻度设置
xlim(ax,[-2.05 d.M+3.25]); ylim(ax,[-.12 headerY+1.17]); axis(ax,'off');
text(ax,-1.96,headerY+1,'Physical-state selection from the allowed quartets', ...
    'Interpreter','latex','FontSize',11,'Color',ink);
text(ax,-1.96,headerY+.64,'Each column is one quartet: both slots are allowed, but a state selects exactly one pair.', ...
    'Interpreter','latex','FontSize',8.5,'Color',gray);
text(ax,-1.96,headerY+.34,'$\mathrm{Im}\,z=+\pi/4\Rightarrow s_j=-1$ (upper), $\mathrm{Im}\,z=-\pi/4\Rightarrow s_j=+1$ (lower).', ...
    'Interpreter','latex','FontSize',8.5,'Color',gray);

rectangle(ax,'Position',[d.M-.42 .20 .84 rowY(1)+.20],'EdgeColor',orange,'FaceColor',[1 .975 .95]);
for j=1:d.M                                            % 循环目的：标出每个允许四重组编号
    color=ink; if j==d.M, color=orange; end
    text(ax,j,headerY,sprintf('$%d$',j),'Interpreter','latex', ...
        'HorizontalAlignment','center','Color',color);
end
text(ax,statusX,headerY,'$P=\prod_j s_j$','Interpreter','latex','HorizontalAlignment','center');
for row=1:nRows                                        % 每个已保存的示例一行；行数和宇称标签都随数据改变。
    y=rowY(row); color=blue; if ~d.isPhysical(row), color=red; end
    text(ax,-1.96,y,d.rowLabels{row},'Interpreter','latex','Color',ink);
    scatter(ax,1:d.M,repmat(y+offset,1,d.M),23,'MarkerEdgeColor',light);
    scatter(ax,1:d.M,repmat(y-offset,1,d.M),23,'MarkerEdgeColor',light);
    selectedY=y-offset*d.selectedSigns(row,:);          % s=-1 取上槽，s=+1 取下槽
    scatter(ax,1:d.M,selectedY,25,'filled','MarkerFaceColor',color);
    text(ax,statusX,y+.11,d.status{row},'HorizontalAlignment','center','Color',color);
    text(ax,statusX,y-.17,sprintf('$P=%+d$',d.parity(row)), ...
        'Interpreter','latex','HorizontalAlignment','center','Color',color);
end
text(ax,-1.96,.03,sprintf('A physical state selects $M=N+1=%d$ pairs; $P=\\kappa=%+d$.',d.M,d.kappa), ...
    'Interpreter','latex','FontSize',8.4);
xx_export_figure(fig,figs,sprintf('fig1c_parity_selection_%s',profile));
end
