function out=xx_spin_spectrum(model)
%XX_SPIN_SPECTRUM 独立自旋对角化；输入 xx_make_model 输出，返回 spectrum/E0/diagnostics。
% 正文(1),(2)：只保留稀疏组装一路，小 N 时求完整 2^N 能谱。
N=model.N; assert(N<=12,'Full spin diagonalization is intended for N<=12.');
sx=sparse([0 1;1 0]); sy=sparse([0 -1i;1i 0]); sz=sparse([1 0;0 -1]);
H=sparse(2^N,2^N); bond=model.J*(kron(sx,sx)+kron(sy,sy));
for j=1:N-1                                  % 加入 N-1 条相邻 XX 交换项。
    H=H+kron(kron(speye(2^(j-1)),bond),speye(2^(N-j-1)));
end
pauli={sx,sy,sz};
for a=1:3                                    % 将左右端三个场分量分别嵌入物理端点。
    H=H+kron(model.hLeft(a)*pauli{a},speye(2^(N-1))) ...
        +kron(speye(2^(N-1)),model.hRight(a)*pauli{a});
end
raw=eig(full(H));                             % 不作能量平移：比较绝对谱。
out=struct('spectrum',sort(real(raw)),'E0',min(real(raw)), ...
    'diagnostics',struct('hermiticityError',norm(H-H','fro')/max(1,norm(H,'fro')), ...
    'maxImagSpectrum',max(abs(imag(raw)))));
end
