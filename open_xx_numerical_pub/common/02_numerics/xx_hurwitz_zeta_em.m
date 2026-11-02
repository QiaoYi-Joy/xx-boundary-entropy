function value = xx_hurwitz_zeta_em(s,a)
%XX_HURWITZ_ZETA_EM 输入 s>1、大正数 a；输出 Phi (5)/(S47) 尾补用的 zeta(s,a)。
% Euler--Maclaurin 展开保留 B_2 到 B_16；本包调用时 a>=2500.5。
B=[1/6,-1/30,1/42,-1/30,5/66,-691/2730,7/6,-3617/510];
value=a^(1-s)/(s-1)+a^(-s)/2;                         % 积分项与端点修正。
for k=1:numel(B)                                      % 逐阶加入 B_2 到 B_16 修正。
    n=2*k;
    value=value+B(k)/factorial(n)*prod(s:s+n-2)*a^(-s-n+1);
end
end
