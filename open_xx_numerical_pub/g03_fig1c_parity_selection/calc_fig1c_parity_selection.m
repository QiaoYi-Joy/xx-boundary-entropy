function data=calc_fig1c_parity_selection(cfg)
%CALC_FIG1C_PARITY_SELECTION 输入模型和任意组翻转位置；输出宇称示意数据，正文(10),(11)、SM(S28)。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'common','01_physics'));
[~,spec]=xx_model_parameters('ND-ND',cfg.parameters); M=cfg.parameters.N+1;
signs=ones(numel(cfg.exampleFlips),M); labels=cell(numel(cfg.exampleFlips),1);
for row=1:numel(cfg.exampleFlips)                    % 每行选择一组上槽，其余列取下槽；标签与实际选择同步。
    flips=cfg.exampleFlips{row};
    assert(all(flips>=1 & flips<=M & flips==fix(flips)) && numel(unique(flips))==numel(flips));
    signs(row,flips)=-1;                            % Im z=+π/4 对应负能量符号。
    if isempty(flips)
        labels{row}='All lower slots';
    else
        names={'One flip','Two flips'}; prefix='Flips';
        if numel(flips)<=2, prefix=names{numel(flips)}; end
        labels{row}=sprintf('%s: j=%s',prefix,strjoin(cellstr(string(flips(:).')),','));
    end
end
parity=prod(signs,2); physical=parity==spec.physicalParity;
status=repmat({'Excluded'},size(parity)); status(physical)={'Allowed physical state'};
data.meta=struct('version','v13r6','revision','r1','panels',{{'Fig. 1(c)'}}, ...
    'profile',cfg.profile,'units','dimensionless','equations',{{'(10)-(13); (S27)-(S28)'}},'numerics',struct);
data.plotData=struct('parameters',cfg.parameters,'M',M,'kappa',spec.physicalParity, ...
    'selectedSigns',signs,'parity',parity,'isPhysical',physical,'rowLabels',{labels},'status',{status});
data.checkData=struct('selectedPairCount',sum(abs(signs)==1,2));
end
