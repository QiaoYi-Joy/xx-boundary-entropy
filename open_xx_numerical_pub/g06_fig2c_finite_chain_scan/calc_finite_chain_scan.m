function data=calc_finite_chain_scan(cfg)
%CALC_FINITE_CHAIN_SCAN 一次根谱扫描供 Fig. 2(c)、S3(a,b) 共用；不保存或画图。
% 输入 cfg：Nlist、temperature、fixedTemperature、参数及求根选项；输出三层 data。
% v13r6 (19),(22),(S66),(S73),(S88)-(S90)：热力学参照、同 N 基态扣除、独立谱核验。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'common','01_physics'),fullfile(root,'common','02_numerics'));
assert(cfg.J==1,'本图组的温度和热力学参照采用 J=1 单位。');
Nlist=cfg.Nlist(:); T=cfg.temperature(:); fixedT=cfg.fixedTemperature(:);
selectedN=intersect(cfg.displayN(:),Nlist,'stable');
nN=numel(Nlist); nT=numel(T);
Rfinite=NaN(nT,numel(selectedN)); correction=NaN(nN,numel(fixedT));
modeError=NaN(nN,2); residual=modeError; counts=modeError; ratioError=NaN(nN,1);
Rthermo=zeros(nT,1); fixedThermo=zeros(numel(fixedT),1); gridSize=zeros(nN,1);
integralOptions=struct('absTol',2e-12,'relTol',2e-11);
for k=1:nT                                         % 逐温度求精确热力学参照，不依赖有限链根谱。
    Rthermo(k)=xx_boundary_ratio(cfg.bLeft,T(k),integralOptions);
end
for k=1:numel(fixedT)                               % S3(a) 的三个固定温度使用同一核积分。
    fixedThermo(k)=xx_boundary_ratio(cfg.bLeft,fixedT(k),integralOptions);
end
parameters=struct('J',cfg.J,'aLeft',cfg.aLeft,'aRight',cfg.aRight, ...
    'bLeft',cfg.bLeft,'bRight',Inf,'phiLeft',cfg.phiLeft,'phiRight',cfg.phiRight);
interval=[]; peakT=NaN; peakR=NaN;
for i=1:nN                                         % 每个 N 只求一次两类边界根谱，同时产生三幅图的数据。
    N=Nlist(i); p=parameters; p.N=N;
    nd=xx_solve_zeros('ND-D',p,cfg.rootOptions); gridSize(i)=nd.numerics.nGrid;
    dd=xx_solve_zeros('D-D',p,cfg.rootOptions);
    eND=sort(nd.epsilon(:)); eDD=sort(dd.epsilon(:));
    assert(numel(eND)==N && numel(eDD)==N,'有限链正能模数必须等于 N。');
    counts(i,:)=[numel(eND),numel(eDD)];
    residual(i,:)=[nd.diagnostics.maxLogEquationResidual,dd.diagnostics.maxLogEquationResidual];
    R=xx_finite_ratio(eND,eDD,[T;fixedT]);            % (S89)：先解析扣除各自同 N 基态，避免大数相消。
    j=find(selectedN==N,1);
    if ~isempty(j), Rfinite(:,j)=R(1:nT); end        % 只保留 Fig. 2(c) 真正使用的温度曲线。
    correction(i,:)=R(nT+1:end).'-fixedThermo.';     % 物理有限尺寸差，不是求解器残差。
    mND=xx_majorana_modes(xx_make_model(N,cfg.J,'ND-D',p));
    mDD=xx_majorana_modes(xx_make_model(N,cfg.J,'D-D',p));
    modeError(i,:)=[max(abs(eND-mND.epsilon)),max(abs(eDD-mDD.epsilon))]; % S3(b) 逐模比较。
    Rmajorana=xx_finite_ratio(mND.epsilon,mDD.epsilon,[T;fixedT]);
    ratioError(i)=max(abs(R-Rmajorana));              % 独立谱在实际观测量上的最大差。
    if i==nN
        [interval,peakT,peakR]=plateau_interval(eND,eDD,cfg.target);
    end
    fprintf('g05: N=%d, mode difference %.3e, R_N difference %.3e\n', ...
        N,max(modeError(i,:)),ratioError(i));
end
rootOptions=nd.numerics; rootOptions.nGrid=gridSize; % 网格数按 plotData.Nlist 排列，其余设置共用。
data.meta=struct('version','v13r6','revision','r1','panels',{{'Fig. 2(c)','Fig. S3(a)','Fig. S3(b)'}}, ...
    'profile',cfg.profile,'units','J=k_B=1','equations',{{'(19),(22); (S66),(S73),(S88)-(S90)'}},'numerics',struct('roots',rootOptions,'integral',integralOptions));
data.plotData=struct('parameters',parameters,'Nlist',Nlist,'displayN',selectedN, ...
    'temperature',T,'fixedTemperature',fixedT,'target',cfg.target, ...
    'finiteRatio',Rfinite,'thermodynamicRatio',Rthermo, ...
    'finiteSizeDifference',correction,'modeEnergyDifference',modeError, ...
    'fivePercentInterval',interval);
data.checkData=struct('modeCountByN',counts,'maxLogResidualByN',residual, ...
    'maxRatioDifferenceByN',ratioError,'largestNPeak',[peakT,peakR]);
end

%% 第一类：文章公式及其必要的数值提升实现
function [interval,peakT,peakR]=plateau_interval(eND,eDD,target)
%PLATEAU_INTERVAL 用 (22)/(S89) 找最大链的峰值及 5% 平台；输入两谱和目标，输出温区与峰值。
ratio=@(t)xx_finite_ratio(eND,eDD,t);
[x,negativePeak]=fminbnd(@(x)-ratio(exp(x)),log(1e-5),log(.8), ...
    optimset('TolX',2e-12,'Display','off'));
peakT=exp(x); peakR=-negativePeak; interval=[];
if peakR>=.95*target && peakR<=1.05*target
    g=@(t)ratio(t)-.95*target;                       % 本组曲线从下方逼近半对数，平台两端均由 95% 决定。
    interval=[fzero(g,[1e-5,peakT]),fzero(g,[peakT,.8])];
end
end
