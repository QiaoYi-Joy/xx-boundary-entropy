function data=calc_fig2c_weak_field_crossover(cfg)
%CALC_FIG2C_WEAK_FIELD_CROSSOVER 输入 b 与缩放温度网格；输出新 Fig. 2(b) 的 R(b,T)。
% v13r6 (19)-(20)/(S64),(S66),(S68),(S87)；只用正核积分，不做低温大数相减。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'common','01_physics'),fullfile(root,'common','02_numerics'));
s=cfg.scaledTemperature(:); b=cfg.bValues(:).'; parameters=struct('bValues',b);
opt=struct('absTol',cfg.absTol,'relTol',cfg.relTol);
ratio=zeros(numel(s),numel(b));
for ib=1:numel(b)                                     % 比较不同弱横场参数 b。
    for it=1:numel(s)                                 % 在相同缩放温度 s=T exp(2b) 上取值。
        ratio(it,ib)=xx_boundary_ratio(b(ib),s(it)*exp(-2*b(ib)),opt);
    end
end
halfScale=zeros(size(b));
for ib=1:numel(b)                                     % 求 R=ln(2)/4 的中点，检查共同交叉能标。
    f=@(x) xx_boundary_ratio(b(ib),exp(x-2*b(ib)),opt)-cfg.target/2;
    halfScale(ib)=exp(fzero(f,log(s([1 end]))));
end
data.meta=struct('version','v13r6','revision','r1','panels',{{'Fig. 2(b)'}}, ...
    'profile',cfg.profile,'units','J=k_B=1','equations',{{'(19)-(20); (S64),(S66)-(S70),(S87)'}},'numerics',opt);
data.plotData=struct('parameters',parameters,'scaledTemperature',s, ...
    'target',cfg.target,'halfTarget',cfg.target/2,'boundaryRatio',ratio);
data.checkData=struct('halfScaledTemperature',halfScale, ...
    'maximumMonotonicIncrease',max(diff(ratio,1,1),[],'all'));
end
