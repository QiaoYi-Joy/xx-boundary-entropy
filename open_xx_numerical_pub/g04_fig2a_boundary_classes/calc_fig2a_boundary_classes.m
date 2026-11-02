function data=calc_fig2a_boundary_classes(cfg)
%CALC_FIG2A_BOUNDARY_CLASSES 输入温度和边界参数；输出三类 (F_b-E_b)/T，正文(4)-(6),(15)-(16)。
% 每个温度只算4个独立 Phi 块，再按 (S47)-(S49),(S58),(S63) 组合；不重复计算。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'common','01_physics'),fullfile(root,'common','02_numerics'));
T=cfg.temperature(:); assert(all(T>0) && issorted(T));
parameters=struct('aLeft',cfg.aLeft,'aRight',cfg.aRight,'bLeft',cfg.bLeft,'bRight',cfg.bRight);
p=[cfg.bLeft,1i*(pi/2-cfg.aLeft),cfg.bRight,1i*(pi/2-cfg.aRight)];
sets={[1 2 3 4],[1 2 4],[2 4]};                    % ND-ND、ND-D、D-D 的参数多重集。
opt=struct('nCut',2500,'tailOrder',7,'minTermsPerBeta',80);
phi=zeros(numel(T),4); energy=xx_boundary_energy(p); % e(p) 不随温度变化，只计算一次。
for it=1:numel(T)                                  % 每个温度的4个独立块只求一次。
    for k=1:4
        phi(it,k)=xx_phi_hat(p(k),1/T(it),opt);
    end
end
ratio=zeros(numel(T),3); targets=[log(2),log(2)/2,0];
for ic=1:3                                        % 按边界类别组合已计算的块。
    indices=sets{ic};
    ratio(:,ic)=real(-sum(phi(:,indices),2)-sum(energy(indices))./T ...
        +log1p(exp(-4./T))-log(2));                 % (15)-(16)：公共 2/T 已解析抵消。
end
fine=struct('nCut',5000,'tailOrder',9,'minTermsPerBeta',160); difference=zeros(1,4);
for k=1:4                                         % 最低温时加长前缀、提高尾阶，检查实际求和收敛。
    difference(k)=abs(phi(1,k)-xx_phi_hat(p(k),1/T(1),fine));
end
data.meta=struct('version','v13r6','revision','r1','panels',{{'Fig. 2(a)'}}, ...
    'profile',cfg.profile,'units','J=k_B=1','equations',{{'(4)-(6),(15)-(16); (S47)-(S49),(S58),(S63)'}}, ...
    'numerics',struct('phi',opt,'convergenceCheck',fine));
data.plotData=struct('parameters',parameters,'temperature',T,'classes',{{'ND--ND','ND--D','D--D'}}, ...
    'targets',targets,'boundaryRatio',ratio);
data.checkData=struct('maxPhiConvergenceDifference',max(difference), ...
    'lowestTemperatureOffset',ratio(1,:)-targets);    % 有限温度相对极限的偏移，不当作数值误差。
end
