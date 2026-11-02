function epsilon = xx_mode_energy(mu,J,eta)
%XX_MODE_ENERGY 输入普通零代表元 mu,J,eta，输出正单模半劈裂；正文(10)、SM(S24)。
if nargin<2 || isempty(J), J=1; end
if nargin<3 || isempty(eta), eta=1i*pi/2; end
mu=mu(:);                                             % 每个零根代表元占一行。
cothLocal=@(z) cosh(z)./sinh(z);
epsilon=J*real(abs(-1i*(cothLocal(mu)+ ...
    cothLocal(-mu-eta))));                            % epsilon(mu)=J|-i[coth(mu)+coth(-mu-eta)]|。
end
