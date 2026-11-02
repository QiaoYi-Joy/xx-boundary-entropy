function data=calc_figS2c_fixedN_root(cfg)
%CALC_FIGS2C_FIXEDN_ROOT 输入固定 N 及 b 扫描，输出最远普通根轨迹；正文(8)、SM(S19)。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'common','01_physics'));
b=cfg.b(:)'; xMax=zeros(size(b)); residual=xMax; count=xMax;
for k=1:numel(b)                                     % 逐个右端 b 求解同一固定 N 的有限零方程。
    p=cfg.parameters; p.bRight=b(k);
    z=xx_solve_zeros('ND-ND',p,cfg.rootOptions);
    xMax(k)=max(real(z.quartetRepresentatives));       % 最远代表元，不用热力学密度峰代替。
    residual(k)=z.diagnostics.maxLogEquationResidual;
    count(k)=numel(z.epsilon);
end
data.meta=struct('version','v13r6','revision','r1','panels',{{'Fig. S2(c)'}}, ...
    'profile',cfg.profile,'units','dimensionless','equations',{{'(8),(17); (S19),(S24)'}},'numerics',z.numerics);
data.plotData=struct('parameters',cfg.parameters,'b',b, ...
    'xMax',xMax,'offset',xMax-b/2,'lastPointGuide',xMax(end)-b(end)/2);
data.checkData=struct('maxLogResidual',max(residual),'modeCount',count);
end
