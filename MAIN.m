%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
val = 0.2;    % control radius to plot
LW = 2;        % line width
factor = 3;    % fatigue life scatter factor
PolDeg = 1;    % degree of the polynomial
CP_name = 'FS';    % Cp factor to be used (SWT or FS are implemented at the moment)
ControlRadius = 'CONST'; % Control radius to be used ('CONST' constant or 'VAR' variable with the number of cycles to failure)
kFS  = 0.4;      % Fatemi Socie
Sy = 606.2;      % Fatemi Socie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  PRE-ALLOCATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CalibFiles = dir(['CALIBRATION', '\*.txt']);
CorrFiles = dir(['CORRECTOR', '\*.txt']);
Calibration = cell(9,length(CalibFiles)/2);
% Contral radius variable definition
rc = linspace(0,0.5,500);
NumCycles = logspace(4,7,1000);
CPrcNf = zeros(length(NumCycles),length(rc),length(CalibFiles)/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  CALIBRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = summer(4);         % colormap
bb = copper(length(CalibFiles)/2);        % colormap
index = length(CalibFiles)/2 + 1;
j = 0;
q = 1;
for i = 1 : length(CalibFiles)/2
    Calibration{1,i} = importdata([CalibFiles(i).folder,'\', CalibFiles(i).name]);
    Calibration{2,i} = importdata([CalibFiles(index).folder,'\', CalibFiles(index).name]);
    LoadRatio = extractAfter(CalibFiles(i).name,"R_");
    Calibration{3,i} = str2double(LoadRatio(1:end-4));
    [Row1, Columns1] = size(Calibration{1,i});
    [Row2, Columns2] = size(Calibration{2,i});
    Calibration{4,i} = Row2;
    Calibration{5,i} = Columns2;
    switch CP_name
        case 'SWT'
            if Row1 == 1
                Calibration{6,i} = polyval(Calibration{1,i},rc);
                Calibration{7,i} = Calibration{2,i}(:,3).^2.*Calibration{6,i}(1,:);
            else
                for kk = 1 : Calibration{4,i}
                    CALIBfitCP = polyval(Calibration{1,i}(kk,:),rc);
                    Calibration{7,i}(kk,:) = CALIBfitCP;
                end
            end

        case 'FS'
            if Row1 == 2 % Can be critic if I have a multiaxial dataset with only two data points
                Calibration{6,i}(1,:) = polyval(Calibration{1,i}(1,:),rc); % Delta gamma
                Calibration{6,i}(2,:) = polyval(Calibration{1,i}(2,:),rc); % Max normal stress
                Calibration{7,i} = 0.5*(Calibration{2,i}(:,3).*Calibration{6,i}(1,:) + Calibration{2,i}(:,3).^2.*(kFS/Sy).*Calibration{6,i}(1,:).*Calibration{6,i}(2,:)); % FS
            else
                for kk = 1 : Calibration{4,i}
                    CALIBfitCP = polyval(Calibration{1,i}(kk,:),rc);
                    Calibration{7,i}(kk,:) = 0.5*(CALIBfitCP);
                end
            end
            [a,b] = size(CorrFiles);
            if a == 0
            else
                if isfile([CorrFiles(1).folder,'\','3_Corr_', extractAfter(CalibFiles(i).name,"1_")]) % Add the correction for calculating Delta Gamma
                    Corrector = importdata([CorrFiles(1).folder,'\','3_Corr_', extractAfter(CalibFiles(i).name,"1_")]);
                    Calibration{7,i} = Corrector.*Calibration{7,i};
                else
                    % File does not exist.
                end
            end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%  PRE-ALLOCATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Linear regression
    j = j + Calibration{4,i};
    LogNf(q:j,2) = log(Calibration{2,i}(:,1));
    LogNf(q:j,1) = 1;
    LogCP(q:j,:) = log(Calibration{7,i});

    SOL = pinv(LogNf(q:j,:))*LogCP(q:j,:);
    Calibration{10,i} = exp(SOL(1,:)); % c
    Calibration{11,i} = SOL(2,:); % m
    Calibration{12,i} = Calibration{10,i}.*NumCycles'.^Calibration{11,i}; % P50% function of rc
    % Fatigue life prediction
    CP(q:j,:) = Calibration{7,i};
    Nf(q:j,1) = Calibration{2,i}(:,1);
    q = q + Calibration{4,i};
    index = index + 1;
    % CP(NumCycles, rc) matrix
    CPrcNf(:,:,i) = Calibration{12,i};

end
%% CP(NumCycles, rc)
Diff = std(CPrcNf, 0, 3);
[VAL, Indexrc] = min(Diff,[],2);
rvsNf = rc(Indexrc);
% Fit power law %
f = fit(NumCycles',rvsNf','C*x^B');

% Fatigue life prediction
SOL = pinv(LogNf)*LogCP;
c = exp(SOL(1,:));
m = SOL(2,:);

Nexpected = (CP./c).^(1./m);
Diff = abs(Nf - Nexpected);
% Fatigue life prediction
j = 0;
q = 1;
for i = 1 : length(CalibFiles)/2
    j = j + Calibration{4,i};
    Calibration{8,i} = (Calibration{7,i}./c).^(1./m);
    Residual(q:j,:) = log(Calibration{2,i}(:,1)) - log(Calibration{8,i});
    q = q + Calibration{4,i};
end
% Residual and standard deviation as a fucntion of rc
S = std(Residual);
d_10 = 1.282;
d_90 = -1.282;
c_10 = c./exp(S*d_10.*m);
c_90 = c./exp(S*d_90.*m);

%Nnum = linspace(0.01*min(Nf), 5*max(Nf), 1000);
Nnum = NumCycles;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  PRE-ALLOCATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PrevFiles = dir(['PREVISION', '\*.txt']);
Prevision = cell(9,length(PrevFiles)/2);
bbb = jet(length(PrevFiles)/2);          % colormap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  PREVISION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index = length(PrevFiles)/2 + 1;
for i = 1 : length(PrevFiles)/2
    Prevision{1,i} = importdata([PrevFiles(i).folder,'\', PrevFiles(i).name]);
    Prevision{2,i} = importdata([PrevFiles(index).folder,'\', PrevFiles(index).name]);
    LoadRatio = extractAfter(PrevFiles(i).name,"R_");
    Prevision{3,i} = str2double(LoadRatio(1:end-4));
    [Row1, Columns1] = size(Prevision{1,i});
    [Row2, Columns2] = size(Prevision{2,i});
    Prevision{4,i} = Row2;
    Prevision{5,i} = Columns2;
    switch CP_name
        case 'SWT'
            if Row1 == 1
                Prevision{6,i} = polyval(Prevision{1,i},rc);
                Prevision{7,i} = Prevision{2,i}(:,3).^2.*Prevision{6,i}(1,:);
            else
                for kk = 1 : Prevision{4,i}
                    PREVfitCP = polyval(Prevision{1,i}(kk,:),rc);
                    Prevision{7,i}(kk,:) = PREVfitCP;
                end
            end
        case 'FS'
            if Row1 == 2 % Can be critic if I have a multiaxial dataset with only two data points
                Prevision{6,i}(1,:) = polyval(Prevision{1,i}(1,:),rc); % Delta gamma
                Prevision{6,i}(2,:) = polyval(Prevision{1,i}(2,:),rc); % Max normal stress
                Prevision{7,i} = 0.5*(Prevision{2,i}(:,3).*Prevision{6,i}(1,:) + Prevision{2,i}(:,3).^2.*(kFS/Sy).*Prevision{6,i}(1,:).*Prevision{6,i}(2,:)); % FS
            else
                for kk = 1 : Prevision{4,i}
                    PREVfitCP = polyval(Prevision{1,i}(kk,:),rc);
                    Prevision{7,i}(kk,:) = 0.5*(PREVfitCP);
                end
            end
            [a,b] = size(CorrFiles);
            if a == 0
            else
                if isfile([CorrFiles(1).folder,'\','3_Corr_', extractAfter(PrevFiles(i).name,"1_")]) % Add the correction for calculating Delta Gamma
                    Corrector = importdata([CorrFiles(1).folder,'\','3_Corr_', extractAfter(PrevFiles(i).name,"1_")]);
                    Prevision{7,i} = Corrector.*Prevision{7,i};
                else
                    % File does not exist.
                end
            end
    end
    Prevision{8,i} = (Prevision{7,i}./c).^(1./m);
    index = index + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% CORRECTION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FIGURES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ControlRadius
    case 'CONST'
        subplot(2,2,1);
        [d, ix] = min(abs(rc-val));
        CPnum = c(ix).*Nnum.^m(ix);
        CPnum_10 = c_10(ix).*Nnum.^m(ix);
        CPnum_90 = c_90(ix).*Nnum.^m(ix);
        Nnum_10 = exp((log(CPnum) - log(c_10(ix)))./(m(ix)));
        factor_10 = max(Nnum_10./Nnum);
        % Graph SN curve with data
        hold on
        for i = 1 : length(CalibFiles)/2
            plot(Calibration{2,i}(:,1), Calibration{7,i}(:,ix), 'o','MarkerFaceColor',bb(i,:),'MarkerEdgeColor',bb(i,:), 'DisplayName', CalibFiles(i).name)
        end

        Tcalib = c_10(ix)/c_90(ix);
        CalibScatter = table(Tcalib, c_10(ix), c(ix), c_90(ix), m(ix));
        writetable(CalibScatter, ['RESULTS/CalibScatter_' CP_name], 'Delimiter', '\t','WriteVariableNames',0)

        plot(Nnum, CPnum, '-black','LineWidth',LW, 'DisplayName', 'P_{50}')
        plot(Nnum, CPnum_10, '--black','LineWidth',LW, 'DisplayName', 'P_{10} - P_{90}')
        plot(Nnum, CPnum_90, '--black','LineWidth',LW, HandleVisibility='off' )
        ylim([0.5*min(CP(:,ix)) 5*max(CP(:,ix))])
        xlim([1E03 1E08])
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        grid on
        xlabel('N_f')
        ylabel([CP_name '_{avg}'])
        legend show
        % 1 - create smaller axes in bottom left, and plot on it
        p = get(gca, 'Position');
        h = axes('Parent', gcf, 'Position', [p(1)+.023 p(2)+.030 p(3)-.24 p(4)-.24]);
        hold on
        for i = 1 : length(CalibFiles)/2
            plot(h, rc, Calibration{7,i}(1,:), 'Color',bb(i,:), 'LineWidth', LW, 'DisplayName', CalibFiles(i).name)
            % Take care that the plot is given just for the first row (i.e. the first load)

        end
        box on
        xlabel('r_c')
        ylabel(CP_name)
        hold off
        % 2 - create smaller axes in top left, and plot on it
        f = axes('Parent', gcf, 'Position', [p(1)+.22 p(2)+.22 p(3)-.24 p(4)-.24]);
        plot(f, rc, S, 'Color', [0 0 0], 'LineWidth', LW)
        box on
        xlabel('r_c')
        ylabel('Residuals (Nf)')
        [minS, IndexS] = min(S);
        Optrc = rc(IndexS);
        text(0.2*max(rc),0.8*max(S),'r_c = ' + string(Optrc))
        hold off
    case 'VAR'
        subplot(2,2,1);
        [d, ix] = min(abs(rc-val));
        CPnum = c(Indexrc).*Nnum.^m(Indexrc);
        CPnum_10 = c_10(Indexrc).*Nnum.^m(Indexrc);
        CPnum_90 = c_90(Indexrc).*Nnum.^m(Indexrc);
        Nnum_10 = exp((log(CPnum) - log(c_10(Indexrc)))./(m(Indexrc)));
        factor_10 = Nnum_10./Nnum;
        CP_VAR = [];
        Nf_VAR = [];
        % Graph SN curve with data
        hold on 
        for i = 1 : length(CalibFiles)/2
            [~, IndexNf] = min(abs(Calibration{2,i}(:,1) - repmat(Nnum, length(Calibration{2,i}(:,1)), 1)),[],2);
            A = rvsNf(IndexNf);
            [~, Calibration{13,i}] = min( abs(A' - repmat(rc, length(A), 1)), [], 2);
            plot(Calibration{2,i}(:,1), diag(Calibration{7,i}(:,Calibration{13,i})), 'o','MarkerFaceColor',bb(i,:),'MarkerEdgeColor',bb(i,:), 'DisplayName', CalibFiles(i).name)
            CP_VAR = [CP_VAR; diag(Calibration{7,i}(:,Calibration{13,i}))];
            Nf_VAR = [Nf_VAR;  Calibration{2,i}(:,1)];
        end
        % Calculate P_10% P_50% and P_90% of the 'VAR' data
        [T_ECP_VAR, c_10_VAR, c_VAR, c_90_VAR, m_VAR] = FatigueScatter(CP_VAR, Nf_VAR, d_10, d_90)


        plot(Nnum, CPnum, '-black','LineWidth',LW, 'DisplayName', 'P_{50}')
        plot(Nnum, CPnum_10, '--black','LineWidth',LW, 'DisplayName', 'P_{10} - P_{90}')
        plot(Nnum, CPnum_90, '--black','LineWidth',LW, HandleVisibility='off' )
        ylim([0.5*min(CP(:,ix)) 5*max(CP(:,ix))])
        xlim([1E03 1E08])
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        grid on
        xlabel('N_f')
        ylabel([CP_name '_{avg}'])
        legend show
        % 1 - create smaller axes in bottom left, and plot on it
        p = get(gca, 'Position');
        h = axes('Parent', gcf, 'Position', [p(1)+.023 p(2)+.030 p(3)-.24 p(4)-.24]);
        hold on
        plot(h,NumCycles, rvsNf, 'LineWidth', LW)
        plot(f,'k')
        box on
        xlabel('N_f')
        ylabel('r_c')
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        hold off

        % Write results to file %
        radiusvsNf = table(rvsNf', NumCycles');
        writetable(radiusvsNf, ['RESULTS/radiusvsNf_' CP_name], 'Delimiter', '\t','WriteVariableNames',0)

end
subplot(2,2,2);
switch ControlRadius
    case 'CONST'
        % Graph SN curve with data
        hold on
        for i = 1 : length(CalibFiles)/2
            plot(Calibration{2,i}(:,1), Calibration{8,i}(:,ix), 'o','MarkerFaceColor',bb(i,:),'MarkerEdgeColor',bb(i,:), 'DisplayName', CalibFiles(i).name)
        end
        plot(Nnum, Nnum, '-black','LineWidth',LW)
        plot(Nnum, factor.*Nnum, '--','Color', [.7 .7 .7],'LineWidth',LW)
        plot(Nnum, (1/factor).*Nnum, '--','Color', [.7 .7 .7],'LineWidth',LW)
        plot(Nnum, factor_10.*Nnum, '--black','LineWidth',LW)
        plot(Nnum, (1/factor_10).*Nnum, '--black','LineWidth',LW)
        grid on
        axis equal
        xlim([1E03 1E08])
        ylim([1E03 1E08])
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        xlabel('N_f')
        ylabel('N_{expected}')
    case 'VAR'
        % Graph SN curve with data
        hold on
        for i = 1 : length(CalibFiles)/2
            plot(Calibration{2,i}(:,1), diag(Calibration{8,i}(:,Calibration{13,i})), 'o','MarkerFaceColor',bb(i,:),'MarkerEdgeColor',bb(i,:), 'DisplayName', CalibFiles(i).name)
        end
        plot(Nnum, Nnum, '-black','LineWidth',LW)
        plot(Nnum, factor.*Nnum, '--','Color', [.7 .7 .7],'LineWidth',LW)
        plot(Nnum, (1./factor).*Nnum, '--','Color', [.7 .7 .7],'LineWidth',LW)
        plot(Nnum, factor_10.*Nnum, '--black','LineWidth',LW)
        plot(Nnum, (1./factor_10).*Nnum, '--black','LineWidth',LW)
        grid on
        axis equal
        xlim([1E03 1E08])
        ylim([1E03 1E08])
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        xlabel('N_f')
        ylabel('N_{expected}')
end
subplot(2,2,3);
switch ControlRadius
    case 'CONST'
        % Graph SN curve with prevision data
        hold on
        for i = 1 : length(PrevFiles)/2
            plot(Prevision{2,i}(:,1), Prevision{7,i}(:,ix), 'o','MarkerFaceColor',bbb(i,:),'MarkerEdgeColor',bbb(i,:), 'DisplayName', PrevFiles(i).name)
        end
        plot(Nnum, CPnum, '-black','LineWidth',LW, 'DisplayName', 'P_{50}')
        plot(Nnum, CPnum_10, '--black','LineWidth',LW, 'DisplayName', 'P_{10} - P_{90}')
        plot(Nnum, CPnum_90, '--black','LineWidth',LW, HandleVisibility='off' )
        ylim([0.5*min(CP(:,ix)) 5*max(CP(:,ix))])
        xlim([1E03 1E08])
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        grid on
        xlabel('N_f')
        ylabel([CP_name '_{avg}'])
        legend show
        % create smaller axes in top right, and plot on it
        p = get(gca, 'Position');
        h = axes('Parent', gcf, 'Position', [p(1)+.023 p(2)+.030 p(3)-.24 p(4)-.24]);
        hold on
        for i = 1 : length(PrevFiles)/2
            plot(h, rc, Prevision{7,i}(1,:), 'Color',bbb(i,:), 'LineWidth', LW, 'DisplayName', PrevFiles(i).name)
            % Take care that the plot is given just for the first row (i.e. the first load)
        end
        box on
        xlabel('r_c')
        ylabel(CP_name)
        hold off
    case 'VAR'
        % Graph SN curve with prevision data
        hold on
        for i = 1 : length(PrevFiles)/2
            [~, IndexNf] = min(abs(Prevision{2,i}(:,1) - repmat(Nnum, length(Prevision{2,i}(:,1)), 1)),[],2);
            A = rvsNf(IndexNf);
            [~, Prevision{13,i}] = min( abs(A' - repmat(rc, length(A), 1)), [], 2);

            plot(Prevision{2,i}(:,1), diag(Prevision{7,i}(:,Prevision{13,i})), 'o','MarkerFaceColor',bbb(i,:),'MarkerEdgeColor',bbb(i,:), 'DisplayName', PrevFiles(i).name)
        end
        plot(Nnum, CPnum, '-black','LineWidth',LW, 'DisplayName', 'P_{50}')
        plot(Nnum, CPnum_10, '--black','LineWidth',LW, 'DisplayName', 'P_{10} - P_{90}')
        plot(Nnum, CPnum_90, '--black','LineWidth',LW, HandleVisibility='off' )
        ylim([0.5*min(CP(:,ix)) 5*max(CP(:,ix))])
        xlim([1E03 1E08])
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        grid on
        xlabel('N_f')
        ylabel([CP_name '_{avg}'])
        legend show
        % create smaller axes in top right, and plot on it
        p = get(gca, 'Position');
        h = axes('Parent', gcf, 'Position', [p(1)+.023 p(2)+.030 p(3)-.24 p(4)-.24]);
        hold on
        for i = 1 : length(PrevFiles)/2
            plot(h, rc, Prevision{7,i}(1,:), 'Color',bbb(i,:), 'LineWidth', LW, 'DisplayName', PrevFiles(i).name)
            % Take care that the plot is given just for the first row (i.e. the first load)
        end
        box on
        xlabel('r_c')
        ylabel(CP_name)
        hold off
end
subplot(2,2,4);
switch ControlRadius
    case 'CONST'
        % Graph SN curve with prevision data
        hold on
        j = 0;
        q = 1;
        for i = 1 : length(PrevFiles)/2
            plot(Prevision{2,i}(:,1), Prevision{8,i}(:,ix), 'o','MarkerFaceColor',bbb(i,:),'MarkerEdgeColor',bbb(i,:), 'DisplayName', PrevFiles(i).name)
            j = j + Prevision{4,i};
            ErrNf(q:j,1) = log10(Prevision{2,i}(:,1)./Prevision{8,i}(:,ix));
            q = q + Prevision{4,i};
        end
        STD = std(ErrNf);
        AVG = mean(ErrNf);
        plot(Nnum, Nnum, '-black','LineWidth',LW)
        plot(Nnum, factor.*Nnum, '--','Color', [.7 .7 .7],'LineWidth',LW)
        plot(Nnum, (1/factor).*Nnum, '--','Color', [.7 .7 .7],'LineWidth',LW)
        plot(Nnum, factor_10.*Nnum, '--black','LineWidth',LW)
        plot(Nnum, (1/factor_10).*Nnum, '--black','LineWidth',LW)
        text(3000,50000000,'Std = ' + string(STD))
        text(3000,20000000,'Avg = ' + string(AVG))
        grid on
        axis equal
        xlim([1E03 1E08])
        ylim([1E03 1E08])
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        xlabel('N_f')
        ylabel('N_{expected}')
    case 'VAR'
        % Graph SN curve with prevision data
        hold on
        j = 0;
        q = 1;
        for i = 1 : length(PrevFiles)/2
            plot(Prevision{2,i}(:,1), diag(Prevision{8,i}(:,Prevision{13,i})), 'o','MarkerFaceColor',bbb(i,:),'MarkerEdgeColor',bbb(i,:), 'DisplayName', PrevFiles(i).name)
            j = j + Prevision{4,i};
            ErrNf(q:j,1) = log10(Prevision{2,i}(:,1)./diag(Prevision{8,i}(:,Prevision{13,i})));
            q = q + Prevision{4,i};
        end
        STD = std(ErrNf);
        AVG = mean(ErrNf);
        plot(Nnum, Nnum, '-black','LineWidth',LW)
        plot(Nnum, factor.*Nnum, '--','Color', [.7 .7 .7],'LineWidth',LW)
        plot(Nnum, (1./factor).*Nnum, '--','Color', [.7 .7 .7],'LineWidth',LW)
        %plot(Nnum, factor_10.*Nnum, '--black','LineWidth',LW)
        %plot(Nnum, (1./factor_10).*Nnum, '--black','LineWidth',LW)
        text(3000,50000000,'Std = ' + string(STD))
        text(3000,20000000,'Avg = ' + string(AVG))
        grid on
        axis equal
        xlim([1E03 1E08])
        ylim([1E03 1E08])
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        xlabel('N_f')
        ylabel('N_{expected}')
end
%% Calculate fatigue scatter for each dataset
switch ControlRadius
    case 'CONST'
        for i = 1 : length(PrevFiles)/2
            [Tsigma, c_10, c, c_90, m] = FatigueScatter(Prevision{7,i}(:,ix), Prevision{2,i}(:,1), d_10, d_90);
            TsigmaTable = table(Tsigma, c_10, c, c_90, m);
            writetable(TsigmaTable, ['RESULTS/FatigueScatter_' CP_name '_' PrevFiles(i).name], 'Delimiter', '\t','WriteVariableNames',0)

        end
        for i = 1 : length(CalibFiles)/2
            [Tsigma, c_10, c, c_90, m] = FatigueScatter(Calibration{7,i}(:,ix), Calibration{2,i}(:,1), d_10, d_90);
            TsigmaTable = table(Tsigma, c_10, c, c_90, m);
            writetable(TsigmaTable, ['RESULTS/FatigueScatter_' CP_name '_' CalibFiles(i).name], 'Delimiter', '\t','WriteVariableNames',0)

        end
    case 'VAR'
        for i = 1 : length(PrevFiles)/2
            [Tsigma, c_10, c, c_90, m] = FatigueScatter(diag(Prevision{7,i}(:,Prevision{13,i})), Prevision{2,i}(:,1), d_10, d_90);
            TsigmaTable = table(Tsigma, c_10, c, c_90, m);
            writetable(TsigmaTable, ['RESULTS/FatigueScatter_' CP_name '_' PrevFiles(i).name], 'Delimiter', '\t','WriteVariableNames',0)
        end
        for i = 1 : length(CalibFiles)/2
            [Tsigma, c_10, c, c_90, m] = FatigueScatter(diag(Calibration{7,i}(:,Calibration{13,i})), Calibration{2,i}(:,1), d_10, d_90);
            TsigmaTable = table(Tsigma, c_10, c, c_90, m);
            writetable(TsigmaTable, ['RESULTS/FatigueScatter_' CP_name '_' CalibFiles(i).name], 'Delimiter', '\t','WriteVariableNames',0)
        end
end
%% Calculate fatigue scatter for the entire dataset without the Calbration Data
CP = [];
Nf = [];
switch ControlRadius
    case 'CONST'
        for i = 1 : length(PrevFiles)/2
            CP = [CP; Prevision{7,i}(:,ix)];
            Nf = [Nf; Prevision{2,i}(:,1)];
        end
        % for i = 1 : length(CalibFiles)/2
        %     CP = [CP; Calibration{7,i}(:,ix)];
        %     Nf = [Nf; Calibration{2,i}(:,1)];
        % end
    case 'VAR'
        for i = 1 : length(PrevFiles)/2
            CP = [CP; diag(Prevision{7,i}(:,Prevision{13,i}))];
            Nf = [Nf; Prevision{2,i}(:,1)];
        end
        % for i = 1 : length(CalibFiles)/2
        %     CP = [CP; diag(Calibration{7,i}(:,Calibration{13,i}))];
        %     Nf = [Nf; Calibration{2,i}(:,1)];
        % end
end
[Tsigma, c_10, c, c_90, m] = FatigueScatter(CP, Nf, d_10, d_90);
TsigmaTable = table(Tsigma, c_10, c, c_90, m);
writetable(TsigmaTable,  ['RESULTS/TotalFatigueScatter_' CP_name], 'Delimiter', '\t','WriteVariableNames',0)
%% Create a table with the data and variable names
switch ControlRadius
    case 'CONST'
        for i = 1 : length(CalibFiles)/2
            Calibration{9,i} = table(Calibration{2,i}(:,1), Calibration{7,i}(:,ix), Calibration{8,i}(:,ix));
            writetable(Calibration{9,i}, ['RESULTS/' CP_name '_' CalibFiles(i).name], 'Delimiter', '\t','WriteVariableNames',0)
        end
        for i = 1 : length(PrevFiles)/2
            Prevision{9,i} = table(Prevision{2,i}(:,1), Prevision{7,i}(:,ix), Prevision{8,i}(:,ix));
            writetable(Prevision{9,i}, ['RESULTS/' CP_name '_' PrevFiles(i).name], 'Delimiter', '\t','WriteVariableNames',0)
        end
        residualsTable = table(rc', S');
        writetable(residualsTable, ['RESULTS/' CP_name '_Residuals'], 'Delimiter', '\t','WriteVariableNames',0)
    case 'VAR'
        for i = 1 : length(CalibFiles)/2
            Calibration{9,i} = table(Calibration{2,i}(:,1), diag(Calibration{7,i}(:,Calibration{13,i})), diag(Calibration{8,i}(:,Calibration{13,i})));
            writetable(Calibration{9,i}, ['RESULTS/' CP_name '_' CalibFiles(i).name], 'Delimiter', '\t','WriteVariableNames',0)
        end
        for i = 1 : length(PrevFiles)/2
            Prevision{9,i} = table(Prevision{2,i}(:,1), diag(Prevision{7,i}(:,Prevision{13,i})), diag(Prevision{8,i}(:,Prevision{13,i})));
            writetable(Prevision{9,i}, ['RESULTS/' CP_name '_' PrevFiles(i).name], 'Delimiter', '\t','WriteVariableNames',0)
        end
        residualsTable = table(rc', S');
        writetable(residualsTable, ['RESULTS/' CP_name '_Residuals'], 'Delimiter', '\t','WriteVariableNames',0)
end