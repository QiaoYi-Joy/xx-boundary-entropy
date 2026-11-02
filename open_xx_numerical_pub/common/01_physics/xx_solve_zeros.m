function sol=xx_solve_zeros(type,p,options)
%XX_SOLVE_ZEROS 通用普通零求解器；输入边界类型、物理参数 p、可选搜索参数。
% 输出模式、代表根、边界弦、必要诊断及实际求根设置；仅作图时另补全四元组。
% 正文(8)--(11)，SM(S19),(S21),(S24),(S36)；只保留稳定路径，边界弦参加有限 N 迭代。
if nargin<3, options=struct; end
[par,spec]=xx_model_parameters(type,p);
assert(~spec.exceptionalZeroMode, ...
    'cos(delta)=0 has an infinite zero mode: take a limit from generic nonzero cos(delta).');
opt=searchOptions(options,par.N);
if spec.isAnalyticFreeFree
    sol=freeChain(par,spec); sol.numerics=opt; return
end
x=opt.xMin+(opt.xMax-opt.xMin)*linspace(0,1,opt.nGrid)'.^2; % 聚集网格解析相位近端的密集根。
u=x+par.eta/2; [A2,A1,~,X,Delta]=coefficients(u,par,spec);
sd=continuousSqrt(Delta);
records=zeros(0,4);                           % 内部临时列：[x,log残差,二次方程残差,相位残差]。
for branchSign=[1,-1]                        % 二次方程两个解析分支均按整数展开相位扫描。
    branch=struct('Z',[],'sqrt',[]);
    rootFactor=(-A1+branchSign*sd)./(2*A2);
    Z=unwrap(angle(X./rootFactor))/(2*pi);
    branch.Z=pchip(x,real(Z)); branch.sqrt=pchip(x,branchSign*sd);
    integers=ceil(min(Z(isfinite(Z)))):floor(max(Z(isfinite(Z))));
    for integer=integers                    % 找出此相位整数的全部交点，允许非单调分支。
        f=Z-integer; intervals=find(f(1:end-1).*f(2:end)<=0);
        for k=intervals(:)'                 % 每个括区间用原对数方程精修，而非只求插值根。
            bracket=x(k:k+1); seed=fzero(@(v) ppval(branch.Z,v)-integer,bracket);
            root=refineLine(seed,bracket,integer,branch,par,spec,opt);
            [logRes,quadRes,phase]=lineResidual(root,integer,branch,par,spec);
            if root>opt.tolEndpointX && isfinite(logRes) && ...
                    logRes<opt.acceptResidual && abs(phase)<opt.acceptPhaseCycles
                records(end+1,:)=[root,logRes,quadRes,phase]; %#ok<AGROW>
            end
        end
    end
end
records=uniqueRecords(records,opt.tolDuplicateX);
expected=spec.MContinuousExpected;
assert(size(records,1)==expected, ...
    '%s N=%d: expected %d line modes, found %d. Adjust search resolution; extra roots are NOT truncated.', ...
    type,par.N,expected,size(records,1));
strings=refineStrings(spec.boundaryStrings,par,spec);
spec.boundaryStrings=strings;
mu=[records(:,1)+par.eta/2;strings.mu(:)];       % 有限 N 连续根和精修弦合为全体有限模。
diagnostics=struct('status','ok','expectedContinuousModes',expected, ...
    'foundContinuousModes',size(records,1),'expectedFiniteModes',spec.MFiniteExpected, ...
    'foundFiniteModes',numel(mu),'boundaryStringCount',strings.count, ...
    'maxLogEquationResidual',max([0;records(:,2)]), ...
    'maxNormalizedQuadraticResidual',max([0;records(:,3)]), ...
    'maxStringResidual',max([0;strings.normalizedResidual]), ...
    'maxStringNewtonCorrection',max([0;strings.newtonCorrection]));
sol=struct('epsilon',xx_mode_energy(mu,par.J,par.eta),'spec',spec, ...
    'quartetRepresentatives',mu,'numerics',opt, ...
    'boundaryStrings',strings,'diagnostics',diagnostics);
end

%% 第一类：文章公式及其稳定实现
function [A2,A1,A0,X,Delta]=coefficients(u,par,spec)
%COEFFICIENTS 输入零变量 u，输出正文(8)、SM(S19)/(S21)的三系数与 X(u)。
if ~strcmp(spec.left.kind,'ND') && ~strcmp(spec.right.kind,'ND')
    aL=imag(spec.left.alpha); aR=imag(spec.right.alpha); v=cosh(2*u);
    cL=cos(2*aL); cR=cos(2*aR); d=cos(2*(aL-aR));
    % SM(S21)代数约化：原 A1 的三个 exp(8 Re u) 项相消，先消去可保护最低能根。
    A2=(v-cL).*(v-cR)/64; A1=(d*v.^2-cL*cR)/32; A0=(v+cL).*(v+cR)/64;
    Delta=-(v.*sinh(2*u)*sin(2*(aL-aR))).^2/1024; % 同一判别式的因式分解，不改变原零方程。
else
    [B,C,chi]=boundaryFactors(u,par,spec);
    [Bm,Cm]=boundaryFactors(u-par.eta,par,spec);
    A2=B.*Cm; A1=chi.*sinh(2*u).^2.*cosh(2*u).^2-B.*Bm-C.*Cm; A0=C.*Bm;
    if nargout>4, Delta=A1.^2-4*A2.*A0; end
end
if nargout>3, X=(sinh(u+par.eta)./sinh(u)).^(2*par.N); end
end

function [B,C,chi]=boundaryFactors(u,par,spec)
%BOUNDARYFACTORS 输入 u/par/spec，输出左右边界乘积；整数 N 指数相位先解析约化。
[BL,CL]=sideFactors(u,spec.left,par.eta);
[BR,CR]=sideFactors(u,spec.right,par.eta);
B=BL.*BR; C=CL.*CR;
S=spec.left.alpha+spec.right.alpha;
if strcmp(spec.left.kind,'ND'), S=S+spec.left.beta; end
if strcmp(spec.right.kind,'ND'), S=S+spec.right.beta; end
if spec.hasParityProjection
    if mod(par.N,2)==0, term=-sinh(S)^2; else, term=cosh(S)^2; end
    chi=(-1)^par.N/4*(cosh(par.DeltaTheta)^2-term); % exp(-2Nη)=(-1)^N，避免大 N 相位舍入。
else
    chi=exp(2*S)/16;                              % SM(S21)：整数 N 的指数相位精确等于 -1。
end
end

function [B,C]=sideFactors(u,endpoint,eta)
%SIDEFACTORS 输入一端参数，输出 SM(S19)/(S21)的 B_sigma,C_sigma。
a=endpoint.alpha;
if strcmp(endpoint.kind,'ND')
    B=sinh(u-a).*cosh(u-endpoint.beta);
    C=sinh(u+a+eta).*cosh(u+endpoint.beta+eta);
else
    B=.5*exp(-u).*sinh(u-a); C=.5*exp(u+eta).*sinh(u+a+eta); % 严格 D/F 缩放退化。
end
end

function strings=refineStrings(strings,par,spec)
%REFINESTRINGS 输入热力学弦初值，输出有限 N 方程根及少量收敛证据；SM(S19),(S24)。
strings.normalizedResidual=zeros(strings.count,1); strings.newtonCorrection=zeros(strings.count,1);
for j=1:strings.count                             % 每条弦的实部、虚部均参加 Newton 迭代。
    mu=strings.mu(j);
    for iteration=1:40                           % 热力学位置只是初值，不固定为最终弦位置。
        h=1e-6*max(1,abs(mu));
        derivative=(stringEquation(mu+h,par,spec)-stringEquation(mu-h,par,spec))/(2*h);
        assert(isfinite(derivative) && abs(derivative)>0,'Singular boundary-string iteration.');
        step=stringEquation(mu,par,spec)/derivative;
        if abs(step)<=4*eps(max(1,abs(mu))), break; end
        if abs(step)>.05, step=.05*step/abs(step); end % 限步，避免跨到另一条弦。
        mu=mu-step;
    end
    [value,scale]=stringEquation(mu,par,spec);
    h=1e-6*max(1,abs(mu));
    derivative=(stringEquation(mu+h,par,spec)-stringEquation(mu-h,par,spec))/(2*h);
    strings.normalizedResidual(j)=abs(value)/max(scale,realmin);
    strings.newtonCorrection(j)=abs(value/derivative);
    assert(strings.normalizedResidual(j)<1e-12 && strings.newtonCorrection(j)<2e-13, ...
        'Finite-N boundary string did not converge.');
    strings.mu(j)=mu;
end
if strings.count>1
    assert(abs(diff(strings.mu))>1e-11,'Coincident string seeds require a separate degenerate-root limit.');
end
end

function [value,scale]=stringEquation(u,par,spec)
%STRINGEQUATION 返回零方程除以 X^2 的等价式及系数尺度，避免弦附近 X 溢出。
[A2,A1,A0]=coefficients(u,par,spec);
inverseX=(sinh(u)/sinh(u+par.eta))^(2*par.N);
value=A2+A1*inverseX+A0*inverseX^2; scale=abs(A2)+abs(A1)+abs(A0);
end

function [residual,quadratic,phase]=lineResidual(x,integer,branch,par,spec)
%LINERESIDUAL 输入连续根位置，返回由 SM(S36)定义的 Fig.S2(d)对数残差、二次式残差和相位圈差。
u=x+par.eta/2; [A2,A1,A0,X,Delta]=coefficients(u,par,spec);
s=sqrt(Delta); reference=ppval(branch.sqrt,x);
if abs(-s-reference)<abs(s-reference), s=-s; end
R=(-A1+s)/(2*A2); ratio=sinh(u+par.eta)/sinh(u);
logAbs=2*par.N*log(abs(ratio))-log(abs(R)); % 分别计算对数，避免形成巨大 X/R。
angleY=2*par.N*angle(ratio)-angle(R);
principal=atan2(sin(angleY),cos(angleY))/(2*pi);
phase=principal+round(ppval(branch.Z,x)-principal)-integer;
residual=abs(logAbs+2i*pi*phase);
quadratic=abs(A2*X^2+A1*X+A0)/(abs(A2*X^2)+abs(A1*X)+abs(A0)); % 回代原方程的独立检查。
end

function sol=freeChain(par,spec)
%FREECHAIN 自由链闭式解，输入 par/spec，输出与数值求根一致的结构；SM(S82)--(S84)。
k=pi*(1:floor(par.N/2))'/(par.N+1);
x=.5*acosh(1./cos(k)); base=x+par.eta/2;
mu=repelem(base,2); epsilon=repelem(2*par.J*cos(k),2); % 每个有限位置保留二重模简并。
diagnostics=struct('status','ok','expectedContinuousModes',0,'foundContinuousModes',0, ...
    'expectedFiniteModes',spec.MFiniteExpected,'foundFiniteModes',numel(mu), ...
    'boundaryStringCount',0,'maxLogEquationResidual',0,'maxNormalizedQuadraticResidual',0, ...
    'maxStringResidual',0,'maxStringNewtonCorrection',0);
strings=spec.boundaryStrings; strings.normalizedResidual=[]; strings.newtonCorrection=[];
sol=struct('epsilon',epsilon,'spec',spec,'quartetRepresentatives',mu, ...
    'boundaryStrings',strings, ...
    'diagnostics',diagnostics);
end

%% 第二类：辅助数值计算
function opt=searchOptions(user,N)
%SEARCHOPTIONS 输入可选搜索设置，输出唯一一组默认参数；不保留重复旧字段名。
opt=struct('xMin',1e-10,'xMax',8,'nGrid',max(4097,40*N+1), ...
    'tolDuplicateX',1e-11,'tolRefineX',1e-14,'tolEndpointX',1e-9, ...
    'acceptResidual',1e-9,'acceptPhaseCycles',1e-9);
names=fieldnames(user);
for k=1:numel(names)                              % 逐项覆盖显式设置，拼写错误即报错。
    assert(isfield(opt,names{k}),'Unknown root option: %s.',names{k});
    opt.(names{k})=user.(names{k});
end
validateattributes(opt.nGrid,{'numeric'},{'scalar','integer','>=',3});
assert(opt.xMin>=0 && opt.xMax>opt.xMin,'Invalid root interval.');
end

function values=continuousSqrt(discriminant)
%CONTINUOUSSQRT 输入网格判别式，输出连续选支平方根，防止主值支切线造成跳变。
values=sqrt(discriminant(:));
for k=2:numel(values)                             % 每点选与前一点更接近的平方根符号。
    if abs(-values(k)-values(k-1))<abs(values(k)-values(k-1)), values(k)=-values(k); end
end
end

function root=refineLine(seed,bracket,integer,branch,par,spec,opt)
%REFINELINE 输入括区间及初值，输出原方程精修根；相位求根后完整复残差最小化。
root=seed; best=lineResidual(root,integer,branch,par,spec);
settings=optimset('TolX',opt.tolRefineX,'Display','off','MaxIter',400,'MaxFunEvals',800);
try
    trial=fzero(@(x) phaseValue(x,integer,branch,par,spec),bracket,settings);
    value=lineResidual(trial,integer,branch,par,spec);
    if value<best, root=trial; best=value; end
catch
    % 原相位在浮点端点失去异号时，保留已有初值并进行完整残差精修。
end
trial=fminbnd(@(x) lineResidual(x,integer,branch,par,spec)^2,bracket(1),bracket(2),settings);
if lineResidual(trial,integer,branch,par,spec)<best, root=trial; end
end

function phase=phaseValue(x,integer,branch,par,spec)
%PHASEVALUE 输入位置及目标整数，输出相位差供 fzero 使用。
[~,~,phase]=lineResidual(x,integer,branch,par,spec);
end

function records=uniqueRecords(records,tolerance)
%UNIQUERECORDS 输入根记录，输出排序去重记录；重复位置只保留残差较小者。
records=sortrows(records,1); keep=true(size(records,1),1);
for k=1:size(records,1)-1                         % 已排序，离开容差后不再比较远处根。
    if ~keep(k), continue; end
    j=k+1;
    while j<=size(records,1) && records(j,1)-records(k,1)<tolerance
        if records(j,2)<records(k,2), keep(k)=false; break; else, keep(j)=false; end
        j=j+1;
    end
end
records=records(keep,:);
end
