function [value,check]=xx_finite_ratio(epsND,epsDD,T)
%XX_FINITE_RATIO 输入同 N 的 ND-D、D-D 非负模式能量和温度；输出列向量 R_N。
% v13r6 (22)/(S73)；用 (S89) 先解析消去各自 E0，避免低温大数相减。
epsND=epsND(:); epsDD=epsDD(:); T=T(:);
assert(numel(epsND)==numel(epsDD) && all(epsND>=0) && all(epsDD>=0));
assert(all(T>0),'温度必须为正。');
value=zeros(size(T));
for k=1:numel(T)                                     % 每个温度分别累加两套能谱的热激发贡献。
    value(k)=sum(log1p(exp(-2*epsDD/T(k))))-sum(log1p(exp(-2*epsND/T(k))));
end
if nargout>1, check=struct('deltaE0',sum(epsDD)-sum(epsND)); end
end
