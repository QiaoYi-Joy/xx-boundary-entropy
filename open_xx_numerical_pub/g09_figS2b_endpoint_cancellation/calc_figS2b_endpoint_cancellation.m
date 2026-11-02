function data=calc_figS2b_endpoint_cancellation(cfg)
%CALC_FIGS2B_ENDPOINT_CANCELLATION 输入 b、beta；输出 Fig. S2(b) 移动核及扣除端点后的值。
% v13r6 (18)/(S44)-(S45),(S50),(S53)；checkData 只保留最大恒等式误差及末端偏差。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'common','01_physics'),fullfile(root,'common','02_numerics'));
b=cfg.b(:); target=log(2)/2; integralValue=zeros(size(b)); identityError=0;
parameters=struct('beta',cfg.beta);                    % 固定逆温与扫描横坐标 b 分开保存。
opt=struct('absTol',cfg.absTol,'relTol',cfg.relTol);
for k=1:numel(b)                                      % 跟踪核峰随 b 移动，独立级数核验积分恒等式。
    [integralValue(k),check]=xx_endpoint_kernel(b(k),cfg.beta,opt);
    identityError=max(identityError,abs(check.identityError));
end
data.meta=struct('version','v13r6','revision','r1','panels',{{'Fig. S2(b)'}}, ...
    'profile',cfg.profile,'units','J=k_B=1','equations',{{'(18); (S44)-(S45),(S50),(S53)'}},'numerics',opt);
data.plotData=struct('parameters',parameters,'b',b,'integral',integralValue, ...
    'subtracted',integralValue-target,'target',target);
data.checkData=struct('maxIdentityError',identityError,'lastIntegralError',abs(integralValue(end)-target));
end
