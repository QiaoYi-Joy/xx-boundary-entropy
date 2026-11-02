function roots = xx_complete_quartets(representatives,eta,tolerance)
%XX_COMPLETE_QUARTETS 输入代表元、eta、去重容差，输出所有允许位置；正文(9)。
if nargin<2 || isempty(eta), eta=1i*pi/2; end
if nargin<3 || isempty(tolerance), tolerance=1e-11; end
z=representatives(:);
raw=[z;z+eta;-z-eta;-z];                              % Q(z)={z,z+eta,-z-eta,-z}。
raw=real(raw)+1i*(mod(imag(raw)+pi/2,pi)-pi/2);       % 折回 Im(z) 属于 [-pi/2,pi/2)。
roots=zeros(0,1);
for z=raw.'                                           % 逐个保留尚未出现的周期等价类。
    dz=z-roots;
    dz=real(dz)+1i*(mod(imag(dz)+pi/2,pi)-pi/2);
    if isempty(roots) || all(abs(dz)>=tolerance), roots(end+1,1)=z; end %#ok<AGROW>
end
[~,order]=sortrows([real(roots),imag(roots)],[1 2]);   % 固定作图时的根排列。
roots=roots(order);
end
