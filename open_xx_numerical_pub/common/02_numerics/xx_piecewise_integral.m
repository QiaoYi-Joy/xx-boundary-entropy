function value = xx_piecewise_integral(fun,points,absTol,relTol)
%XX_PIECEWISE_INTEGRAL 输入函数、递增断点和精度；输出分段积分和，用于 (S66)/(S44)。
points=unique(points(:).','stable');
assert(numel(points)>=2 && all(diff(points)>0),'Breakpoints must increase.');
value=0;
for k=1:numel(points)-1                               % 累加每个已知结构区间的贡献。
    value=value+integral(fun,points(k),points(k+1), ...
        'AbsTol',absTol/(numel(points)-1),'RelTol',relTol,'ArrayValued',true);
end
end
