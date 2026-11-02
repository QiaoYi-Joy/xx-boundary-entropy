function [par,spec]=xx_model_parameters(type,p)
%XX_MODEL_PARAMETERS 输入 p.N,J,aLeft,aRight,bLeft,bRight,phiLeft,phiRight，解析九类边界。
% 输出 par（求根参数）和 spec（边界、模数、宇称）；正文(2),(9)--(11)，SM(S24),(S28)。
type=upper(char(type)); parts=strsplit(type,'-');
assert(numel(parts)==2 && all(ismember(parts,{'ND','D','F'})),'Use an ordered ND/D/F pair.');
validateattributes(p.N,{'numeric'},{'scalar','integer','positive'});
J=option(p,'J',1); validateattributes(J,{'numeric'},{'scalar','real','positive','finite'});
left=side(parts{1},p,'Left'); right=side(parts{2},p,'Right');
delta=imag(left.theta-right.theta); hasProjection=strcmp(type,'ND-ND');
exceptional=hasProjection && abs(cos(delta))<32*eps;
kappa=NaN;
if hasProjection
    kappa=sign(cos(delta));                     % 正文(11)：由实际方位角决定，不硬编码为 +1。
    if exceptional, kappa=1; end               % 零模极限两扇区等谱，约定选 +1；从不使用 κ=0。
end
par=struct('N',p.N,'J',J,'eta',1i*pi/2,'DeltaTheta',left.theta-right.theta);
free=strcmp(type,'F-F'); zeroInfinity=double(free && mod(p.N,2)==1);
seeds=[];
for endpoint=[left,right]                      % 逐端检查边界弦；其位置仅作为有限 N 迭代初值。
    a=mod(imag(endpoint.alpha),pi);
    if a<pi/4-1e-10 || a>3*pi/4+1e-10
        seeds(end+1,1)=1i*min(a,pi-a); %#ok<AGROW>
    end
end
strings=struct('mu',seeds,'count',numel(seeds));
M=p.N+double(hasProjection)-zeroInfinity;
continuous=M-numel(seeds); if free, continuous=0; end
spec=struct('type',type,'left',left,'right',right,'MFiniteExpected',M, ...
    'MContinuousExpected',continuous,'isAnalyticFreeFree',free, ...
    'hasParityProjection',hasProjection,'physicalParity',kappa, ...
    'exceptionalZeroMode',exceptional,'zeroModesAtInfinity',zeroInfinity, ...
    'boundaryStrings',strings);
end

%% 第一类：文章公式及其稳定实现
function out=side(kind,p,suffix)
%SIDE 输入端点类别及实参数，输出 alpha=i*a,beta=b,theta=i*phi；正文(2)。
if strcmp(kind,'F'), a=pi/2; else, a=p.(['a' suffix]); end
assert(isreal(a) && isscalar(a) && a>0 && a<pi,'Require 0<a<pi.');
if strcmp(kind,'ND')
    b=p.(['b' suffix]); phi=option(p,['phi' suffix],0);
    assert(isreal(b) && isfinite(b) && isreal(phi) && isfinite(phi),'Require real finite b and phi.');
else
    b=Inf; phi=0;
end
out=struct('kind',kind,'alpha',1i*a,'beta',b,'theta',1i*phi);
end

%% 第二类：辅助数值计算
function value=option(s,name,fallback)
%OPTION 读取唯一标准字段；输入结构、字段名、默认值，输出字段值。
if isfield(s,name), value=s.(name); else, value=fallback; end
end
