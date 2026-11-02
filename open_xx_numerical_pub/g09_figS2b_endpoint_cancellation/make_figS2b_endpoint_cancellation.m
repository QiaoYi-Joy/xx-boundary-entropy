function data=make_figS2b_endpoint_cancellation(profile,recompute)
%MAKE_FIGS2B_ENDPOINT_CANCELLATION 输入 quick/paper、是否重算；保存并画 Fig. S2(b)。
% false 读取本组已存数据；作图仅使用 plotData。
if nargin<1, profile='paper'; end
if nargin<2, recompute=true; end
profile=validatestring(profile,{'quick','paper'}); nB=70;
if strcmp(profile,'quick'), nB=30; end
here=fileparts(mfilename('fullpath')); root=fileparts(here);
addpath(fullfile(root,'common','01_physics'),fullfile(root,'common','02_numerics'));
cfg=struct('profile',profile,'b',linspace(.5,7,nB),'beta',4,'absTol',2e-12,'relTol',2e-11);
stem=['figS2b_endpoint_cancellation_' profile]; folder=fullfile(here,'data');
if ~exist(folder,'dir'), mkdir(folder); end
file=fullfile(folder,[stem '.mat']);
if recompute
    data=calc_figS2b_endpoint_cancellation(cfg);
else
    saved=load(file,'data'); data=saved.data;
end
assert(data.checkData.maxIdentityError<2e-8 && data.checkData.lastIntegralError<2e-4);
if recompute, save(file,'data','-v7'); end
draw(data.plotData,fullfile(here,'figures'),stem);
fprintf('g09 [%s]: identity %.3e, endpoint %.3e\n',profile, ...
    data.checkData.maxIdentityError,data.checkData.lastIntegralError);
end

%% 第二类：辅助绘图函数
function draw(d,folder,stem)
%DRAW 输入 plotData 和输出位置；画移动核与端点扣除，输出 PDF、PNG。
fig=figure('Visible','off','Color','w','Units','inches','Position',[1 1 4.3 3.2]);
ax=axes(fig); hold(ax,'on');
plot(ax,d.b,d.integral,'LineWidth',1.3,'DisplayName','moving-kernel integral');
plot(ax,d.b,d.subtracted,'LineWidth',1.3,'DisplayName','integral - ln2/2');
yline(ax,[d.target 0],':','HandleVisibility','off');       % 原积分趋 ln2/2，扣除后趋零。
xlabel(ax,'Boundary parameter b'); ylabel(ax,'Endpoint weight');
title(ax,'(b) Endpoint cancellation','FontWeight','normal'); legend(ax,'Location','best','Box','off');
xx_style_axes(ax); xx_export_figure(fig,folder,stem);
end
