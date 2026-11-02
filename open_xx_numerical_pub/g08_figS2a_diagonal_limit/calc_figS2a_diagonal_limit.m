function data=calc_figS2a_diagonal_limit(cfg)
%CALC_FIGS2A_DIAGONAL_LIMIT 输入 b 网格、固定温度 T；输出 Fig. S2(a) 的 Delta F_b。
% v13r6 (5),(16),(19)/(S47),(S64)；Delta F_b=-T PhiHat，T=0 为 e(b)。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'common','01_physics'),fullfile(root,'common','02_numerics'));
b=cfg.b(:); T=cfg.T(:).'; deltaF=zeros(numel(b),numel(T));
parameters=struct('temperature',T);                   % 固定温度集中记录，扫描的 b 单独作为横坐标。
opt=struct('nCut',2500,'tailOrder',7,'minTermsPerBeta',80);
ground=xx_boundary_energy(b);                         % 同一固定 a 比较中的基态能差。
for ib=1:numel(b)                                     % 增大 b，逐步逼近严格对角边界。
    for it=1:numel(T)                                 % 每个固定温度下求精确自由能差。
        deltaF(ib,it)=-T(it)*xx_phi_hat(b(ib),1/T(it),opt);
    end
end
data.meta=struct('version','v13r6','revision','r1','panels',{{'Fig. S2(a)'}}, ...
    'profile',cfg.profile,'units','J=k_B=1','equations',{{'(5),(16),(19); (S47),(S64)'}},'numerics',opt);
data.plotData=struct('parameters',parameters,'b',b,'deltaF',deltaF,'groundEnergy',ground);
data.checkData=struct('maxAbsAtLargestB',max(abs([deltaF(end,:),ground(end)])));
end
