function data=calc_fig1b_ndd_roots(cfg)
%CALC_FIG1B_NDD_ROOTS 输入配置，计算严格 ND--D 允许零点；输出三层数据。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'common','01_physics'));
p=cfg.parameters; s=xx_solve_zeros('ND-D',p,cfg.rootOptions); % 正文(8)、SM(S21)的严格退化方程。
z=xx_complete_quartets(s.quartetRepresentatives); % 只在位置图中补全四元组，不让每次求谱都做此工作。
data.meta=struct('version','v13r6','revision','r1','panels',{{'Fig. 1(b)'}}, ...
    'profile',cfg.profile,'units','J=1','equations',{{'(8)-(10),(17); (S21),(S24)'}},'numerics',s.numerics);
data.plotData=struct('parameters',p,'allowedPositions',z);
data.checkData=struct('modeCount',numel(s.epsilon),'positionCount',numel(z), ...
    'maxLogResidual',s.diagnostics.maxLogEquationResidual);
end
