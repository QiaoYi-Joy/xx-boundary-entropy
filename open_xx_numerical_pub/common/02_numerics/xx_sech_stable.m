function y = xx_sech_stable(x)
%XX_SECH_STABLE 输入实数数组 x；输出 sech(x)，用于模式能量及 (S66) 核，避免溢出。
e=exp(-abs(x));
y=2*e./(1+e.^2);                                     % 与 1/cosh(x) 严格等价。
end
