function data=calc_spectrum_scan(cfg)
%CALC_SPECTRUM_SCAN 输入尺寸和九类边界，输出 Fig.S1(a,b)、S2(d)共用的三层数据。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'common','01_physics'));
nN=numel(cfg.Nlist); nC=numel(cfg.classes); error=zeros(nN,nC); residual=error;
rootCount=error; stateCount=error; stringResidual=error; stringCorrection=error;
panel=struct; gridSize=zeros(nN,1);
for i=1:nN                                         % 逐链长比较完整物理能谱。
    p=cfg.parameters; p.N=cfg.Nlist(i);
    for j=1:nC                                     % 九个有序边界组合使用同一个公共求根器。
        type=cfg.classes{j}; z=xx_solve_zeros(type,p,cfg.rootOptions); gridSize(i)=z.numerics.nGrid; % 正文(8)、SM(S19)/(S21)。
        zeroSpectrum=xx_reconstruct_spectrum(z.epsilon,z.spec);    % 正文(10)--(13)：含实际 κ 投影。
        spin=xx_spin_spectrum(xx_make_model(p.N,p.J,type,p));       % 正文(1),(2)：独立自旋 ED。
        assert(numel(zeroSpectrum)==2^p.N,'Physical state count mismatch.');
        error(i,j)=max(abs(zeroSpectrum-spin.spectrum))/p.J;      % 绝对谱差，不平移、不拟合。
        residual(i,j)=z.diagnostics.maxLogEquationResidual;
        rootCount(i,j)=numel(z.epsilon); stateCount(i,j)=numel(zeroSpectrum);
        stringResidual(i,j)=z.diagnostics.maxStringResidual;
        stringCorrection(i,j)=z.diagnostics.maxStringNewtonCorrection;
        if i==nN && strcmp(type,'ND-ND')
            panel=struct('N',p.N,'zeroSpectrum',zeroSpectrum/p.J, ...
                'spinSpectrum',spin.spectrum/p.J);                 % S1(a)保留最大尺寸的一个代表谱。
        end
    end
end
rootOptions=z.numerics; rootOptions.nGrid=gridSize; % 与 Nlist 一一对应，不保存网格数组。
data.meta=struct('version','v13r6','revision','r1','panels',{{'Fig. S1(a,b)','Fig. S2(d)'}}, ...
    'profile',cfg.profile,'units','energies divided by J','equations',{{'(1)-(2),(8)-(13); (S19),(S21),(S24),(S28),(S36),(S86); Fig.S2(d) caption'}},'numerics',rootOptions);
data.plotData=struct('parameters',cfg.parameters,'Nlist',cfg.Nlist,'classes',{cfg.classes}, ...
    'spectrumError',error,'lineRootResidual',residual,'panel',panel);
data.checkData=struct('modeCount',rootCount,'stateCount',stateCount, ...
    'maxStringResidual',max(stringResidual,[],'all'), ...
    'maxStringNewtonCorrection',max(stringCorrection,[],'all'));
end
