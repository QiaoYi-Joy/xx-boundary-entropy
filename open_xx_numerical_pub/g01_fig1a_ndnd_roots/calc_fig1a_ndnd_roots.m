function data=calc_fig1a_ndnd_roots(cfg)
%CALC_FIG1A_NDND_ROOTS 输入配置，计算 ND--ND 允许零点；输出 meta/plotData/checkData。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'common','01_physics'));
p=cfg.parameters; s=xx_solve_zeros('ND-ND',p,cfg.rootOptions); % 正文(8)求根，(9)补全允许位置。
z=xx_complete_quartets(s.quartetRepresentatives); [xStar,kStar]=max(real(s.quartetRepresentatives));
escapeMask=abs(abs(real(z))-xStar)<5e-7; assert(nnz(escapeMask)==4,'Incomplete quartet.');
data.meta=struct('version','v13r6','revision','r1','panels',{{'Fig. 1(a)'}}, ...
    'profile',cfg.profile,'units','J=1','equations',{{'(8)-(10); (S19),(S24)'}},'numerics',s.numerics);
data.plotData=struct('parameters',p,'allowedPositions',z,'escapeMask',escapeMask,'escapeX',xStar);
data.checkData=struct('modeCount',numel(s.epsilon),'positionCount',numel(z), ...
    'maxLogResidual',s.diagnostics.maxLogEquationResidual,'escapeEnergy',s.epsilon(kStar));
end
