function spectrum=xx_reconstruct_spectrum(epsilon,spec)
%XX_RECONSTRUCT_SPECTRUM 输入单模半劈裂及边界分类，输出排序多体能量；正文(10)--(13)、SM(S28)。
epsilon=real(epsilon(:)); M=numel(epsilon);
assert(M<=20,'Use full spectrum reconstruction only for small M.');
if M==0, signs=ones(1,0); else, signs=2*(dec2bin(0:2^M-1,M)-'0')-1; end % 枚举 s_j=±1。
if spec.hasParityProjection
    signs=signs(prod(signs,2)==spec.physicalParity,:); % 只保留物理扇区 prod(s)=κ。
end
spectrum=sort(signs*epsilon);
if spec.zeroModesAtInfinity>0
    spectrum=sort(repelem(spectrum,2^spec.zeroModesAtInfinity)); % 奇数自由链零模保留简并度。
end
end
