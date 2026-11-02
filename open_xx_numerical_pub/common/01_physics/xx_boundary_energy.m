function value=xx_boundary_energy(p)
%XX_BOUNDARY_ENERGY 输入边界参数 p，输出 v13r6 (16) 的 e(p)，支持数组。
value=-4*p./(pi*sinh(2*p));                            % e(p)=-4p/[pi sinh(2p)]。
value(p==0)=-2/pi;                                    % 可去奇点 p=0 的解析极限。
end
