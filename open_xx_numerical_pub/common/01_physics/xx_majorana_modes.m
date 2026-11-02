function out=xx_majorana_modes(model)
%XX_MAJORANA_MODES 独立 Majorana 能谱；输入 xx_make_model 模型，输出正半能级及模数。
% v13r6 (S90)：A=[0,C;-C',0]，故 epsilon_j=svd(C)/2；不形成大矩阵 A。
% 直接求 C 的奇异值，避免 C'*C 将条件数平方、损失大 N 软模的精度。
assert(any(strcmp(model.boundaryType,{'ND-D','D-D'})), ...
    '此 Majorana 接口仅用于 ND-D 和 D-D。');
N=model.N; J=model.J;
rows=[2;(3:N+1).';(2:N).';N+1];
cols=[1;(1:N-1).';(2:N).';N];
values=[2*model.hLeft(3);-2*J*ones(N-1,1); ...
    -2*J*ones(N-1,1);2*model.hRight(3)]; % 纵场在两端，交换项在相邻格点。
if strcmp(model.boundaryType,'ND-D')
    rows=[1;rows]; cols=[1;cols];
    values=[2*hypot(model.hLeft(1),model.hLeft(2));values]; % 横场方位可由全局绕 z 轴旋转消去。
end
C=sparse(rows,cols,values,N+1,N);
if strcmp(model.boundaryType,'D-D'), C=C(2:end,:); end % 对角链不需辅助 Majorana。
out.epsilon=sort(svd(full(C),'econ'))/2;              % 仅求奇异值，不求奇异向量。
out.diagnostics=struct('modeCount',numel(out.epsilon),'minimumEnergy',out.epsilon(1));
end
