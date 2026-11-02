function value=xx_phi_hat(p,beta,opt)
%XX_PHI_HAT 输入实 p 或纯虚 |Im(p)|<pi/2、逆温 beta>=0；输出 (5)/(S47) 的 PhiHat。
% 只保留稳定前缀求和加解析尾补；精度由截断/尾阶变化及独立核积分检验。
if nargin<3, opt=struct; end
if ~isfield(opt,'nCut'), opt.nCut=2500; end
if ~isfield(opt,'tailOrder'), opt.tailOrder=7; end
if ~isfield(opt,'minTermsPerBeta'), opt.minTermsPerBeta=80; end
assert(isscalar(p) && isfinite(p) && (isreal(p) || (real(p)==0 && abs(imag(p))<pi/2)));
assert(isscalar(beta) && isreal(beta) && isfinite(beta) && beta>=0,'beta 必须是有限非负实标量。');
c=real(cosh(2*p)); assert(1+c>0,'p 太接近角参数奇点，双精度无法分辨 1+cosh(2p)。');
if beta==0, value=0; return; end
M=max(opt.nCut,ceil(opt.minTermsPerBeta*beta));
if c<0                                              % 近虚轴端点时尾部泰勒收敛半径缩小。
    radius=min(1,sqrt(1-c*c));                        % 最近奇点给出 q 的收敛半径 sqrt(1-c^2)。
    M=max(M,ceil(20*beta/(pi*radius)));                % 保证首个尾项 q<=0.1 radius。
end
n=(1:M).'; q=4*beta./(pi*(2*n-1));                    % (5) 中的 Matsubara 参数。
prefix=sum(log1p(q.^2./((sqrt(1+q.^2)+1).*(1+c))));   % 有理化避免小 q 时 sqrt(1+q^2)-1 相消。
a=xx_phi_tail_coefficients(c,opt.tailOrder); tail=0;
for k=1:opt.tailOrder                                % 将 n>M 的每个 q^(2k) 尾项解析求和。
    term=a(k+1)*(2*beta/pi)^(2*k)*xx_hurwitz_zeta_em(2*k,M+.5);
    tail=tail+term;                                   % 奇数分母级数转为 Hurwitz zeta。
end
value=prefix+tail;                                    % 允许的实/纯虚参数均化为实数运算。
end
