function field = xx_boundary_field(kind,orientation,J,a,b,phi)
%XX_BOUNDARY_FIELD 输入端点类别、方向号、J,a,b,phi；按正文(2)输出 [h_x,h_y,h_z]。
switch upper(kind)
    case 'ND'
        hPerp=J/(sin(a)*cosh(b));                     % 非对角边界横场强度 h_perp。
        field=[hPerp*cos(phi),-hPerp*sin(phi), ...
            orientation*J*cot(a)*tanh(b)];            % 左/右端由 orientation=+1/-1 区分。
    case 'D'
        field=[0,0,orientation*J*cot(a)];             % 严格对角边界：h_perp=0。
    case 'F'
        field=[0,0,0];                                % 自由边界：边界场为零。
    otherwise
        error('Unknown boundary kind "%s".',kind);
end
end
