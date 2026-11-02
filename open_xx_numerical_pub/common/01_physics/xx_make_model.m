function model=xx_make_model(N,J,boundaryType,parameters)
%XX_MAKE_MODEL 输入 N,J,边界类型及实参数，输出左右物理场；正文(1),(2)。
assert(isscalar(N) && isreal(N) && isfinite(N) && N==round(N) && N>=2, ...
    'N must be an integer >=2.');
assert(isscalar(J) && isreal(J) && isfinite(J) && J>0,'J must be positive.');
assert(isstruct(parameters),'parameters must be a structure.');

parameters.N=N; parameters.J=J;
[~,spec]=xx_model_parameters(boundaryType,parameters); % 复用唯一的边界类型和参数解析入口。
hLeft=physicalField(spec.left,+1,J);                   % 左端场采用正文式(2)的正方向号。
hRight=physicalField(spec.right,-1,J);                 % 右端场采用正文式(2)的负方向号。
model=struct('N',double(N),'J',double(J),'boundaryType',spec.type, ...
    'hLeft',hLeft,'hRight',hRight);
end

function field=physicalField(parsed,orientation,J)
%PHYSICALFIELD 将解析参数还原为 a、b、phi，并计算一个端点的三分量边界场。
a=imag(parsed.alpha);
if strcmp(parsed.kind,'ND'), b=parsed.beta; phi=imag(parsed.theta); else, b=Inf; phi=0; end
field=xx_boundary_field(parsed.kind,orientation,J,a,b,phi); % 正文式(2)。
end
