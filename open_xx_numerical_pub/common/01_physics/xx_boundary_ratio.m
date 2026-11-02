function [value,check]=xx_boundary_ratio(b,T,opt)
%XX_BOUNDARY_RATIO 输入实数 b、正温 T；输出 v13r6 (19)/(S64) 的 R(b,T)。
% (S66) 作 sinh(2x)=cot(u)/eta 换元；避免 Phi 与边界能相消。
if nargin<3, opt=struct; end
if ~isfield(opt,'absTol'), opt.absTol=2e-12; end
if ~isfield(opt,'relTol'), opt.relTol=2e-11; end
assert(isreal(b) && isscalar(b) && isscalar(T) && T>0,'b 为实标量，T 为正标量。');
eta=xx_sech_stable(2*b); r=T/(2*eta);                  % eta=sech(2b)，能标为 2 eta。
points=unique([0,atan(r*[.1 .3 1 3 10 30]),pi/2]);     % 在热核变化区增加积分分段。
value=xx_piecewise_integral(@(u) kernel(u,r,eta)/pi,points,opt.absTol,opt.relTol);
if nargout>1
    check=struct;
    if isfield(opt,'checkPhi') && opt.checkPhi
        check.phiKernelError=-xx_phi_hat(b,1/T,opt)-xx_boundary_energy(b)/T-value;
    end
end
end

%% 第一类：文章公式及必要的稳定实现
function y=kernel(u,r,eta)
%KERNEL 输入换元坐标 u、r=T/(2 eta)；输出 (S66) 正核的温度部分。
z=sin(u)./(r.*hypot(cos(u),eta.*sin(u)));              % z=epsilon/T，hypot 避免平方溢出。
y=log(2)-log1p(exp(-2*z));                             % 等于 z-log(cosh(z))，无低温相消。
end
