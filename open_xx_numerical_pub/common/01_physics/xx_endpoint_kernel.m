function [value,check]=xx_endpoint_kernel(b,beta,opt)
%XX_ENDPOINT_KERNEL 输入实 b、逆温 beta；输出 (18)/(S44)-(S45) 的移动核积分。
% 可选 check.identityError 检验积分=PhiHat+ln(2)/2，即 (S50)/(S53)。
if nargin<3, opt=struct; end
if ~isfield(opt,'absTol'), opt.absTol=2e-12; end
if ~isfield(opt,'relTol'), opt.relTol=2e-11; end
assert(isreal(b) && isscalar(b) && isscalar(beta) && beta>0);
b=abs(b);                                            % 核关于 b 为偶函数，峰位取 |b|。
fun=@(x) (xx_sech_stable(2*(x+b))+xx_sech_stable(2*(x-b))).* ...
    xx_log2cosh_stable(2*beta*xx_sech_stable(2*x))/pi;  % f(x)=ln[2 cosh(beta epsilon)]。
points=unique([0,max(0,b-8),b,b+8,max(20,b+20),Inf]);   % 跟随 x≈b 的移动峰，避免漏采。
value=xx_piecewise_integral(fun,points,opt.absTol,opt.relTol);
if nargout>1, check=struct('identityError',value-xx_phi_hat(b,beta,opt)-log(2)/2); end
end
