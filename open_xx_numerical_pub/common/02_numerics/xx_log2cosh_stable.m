function y = xx_log2cosh_stable(x)
%XX_LOG2COSH_STABLE 输入实数数组 x；输出 log[2 cosh(x)]，用于 (18) 热核且防溢出。
a=abs(x);
y=a+log1p(exp(-2*a));                                 % 与原式严格等价，但指数不超过 1。
end
