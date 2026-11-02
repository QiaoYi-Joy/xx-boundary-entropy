function coefficients = xx_phi_tail_coefficients(c,order)
%XX_PHI_TAIL_COEFFICIENTS 输入 c=cosh(2p) 与阶数；输出 (5)/(S47) 的 q^(2k) 系数。
% coefficients(k+1) 对应 q^(2k)，常数项为零；用于解析补齐大 n 尾部。
u=zeros(1,order+1); binomialHalf=1;
for k=1:order                                         % 展开 [sqrt(1+q^2)-1]/(1+c)。
    binomialHalf=binomialHalf*(0.5-k+1)/k;
    u(k+1)=binomialHalf/(1+c);
end
coefficients=zeros(1,order+1); power=u;
for r=1:order                                         % 合成 log(1+u) 的幂级数。
    coefficients=coefficients+(-1)^(r+1)*power/r;
    power=conv(power,u); power=power(1:order+1);
end
end
