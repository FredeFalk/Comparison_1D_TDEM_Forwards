%% ANALYSIS FOR "Comparison of 1D TDEM Forward Functions"

% This MATLAB code produces the results for Falk et al. (Submitted to Applied Computing & Geosciences 2024)
% The code relies on proper installation of the three forward algorithms. We have
% segmented this code such that you can run the tests that are possible for
% you. Below you'll find some points with respect to the installation of
% each algorith (text taken from the manuscript).

%% INSTALLING SIMPEG

% SimPEG is a framework for simulation and parameter estimation in geophysics, implemented as a Python library.
% The installation of SimPEG is as easy as typing "pip install simpeg" in the Python console.
% We interface to SimPEG using Python scripts which we call from MATLAB using the \verb'pyrunfile' function,
% hich is native to MATLAB versions R2021b or newer. In order to use \verb'pyrunfile' in MATLAB the user must
% first set up the Python environment in MATLAB which requires a valid version of the Python software installed.

%% INSTALLING AARHUSINV

% The AarhusInv 1D forward code works by calling a binary executable file, called "AarhusInv64.exe", which reads
% the layered model from an ASCII text file with a ".mod" extension, in which the ".tem" file to be read is also
% specified. The ".tem" file is another ASCII text file which contains the TEM system specifications. For each
% forward response one single ASCII file is written by AarhusInv to the same directory as the ".mod" and ".tem"
% files, with the extension ".fwr". forward modelling settings may be set in one last ASCII file with the ".con"
% extension. We refer to the AarhusInv users manual for more information \citep{AarhusInvManual}.
% The executable file "AarhusInv64.exe" may be called from a command line, and thus also from within MATLAB,
% given that the user has the appropriate license to do so, which should be activated using another executable
% file called "AarhusInvLic.exe"

% We have a working implementation of AarhusInv, such that a new user of this code would only need to register
% their license in the AarhusInvLic.exe file.

%% INSTALLING GA-AEM

% The GA-AEM forward algorithm is implemented in C++ and comes with interfaces for MATLAB and Python, through the
% use of a library. However, in order to use the code one must install a series of ".dll" files in the windows
% search path. Furthermore, if the user would like to use the MATLAB implementation it is also necessary to set
% the MATLAB search path to the folder containing the GA-AEM ".m" files.

%% CODE BEGINS HERE

close all
clear all
clc

rng(1)

%% Which analyses should be performed? (0 = no, 1 = yes)
prec1 = 0;
prec2 = 0;
time1 = 1;
time2 = 0;

%% Which forward functions should be included? (0 = no, 1 = yes)
UseAarhusInv = 1;
UseSimPEG = 1;
UseGAAEM = 1;

%% SET PARAMETERS %%

%Number of models in precision analysis 1 (halfspace and analytic comparison)
Nmod_prec1 = 100;

%Number of models in precision analysis 2 (forwards relative to each other)
Nmod_prec2 = 100;

%All computation time analyses
%Set number of data points
NDL = 3;

%Computation time analysis 2a (fixed number of models, variable number of layers)
%Set number of models
NmodTest = 10;

%Computation time analysis 2b (1 model at a time, variable number of layers)
%Set number of models
NSingleModTest = 1;

%Computation time analysis 1 (variable number of models, fixed number of layers)
%Set number of layers
NLT = 5;

%Set the maximum number of models for time analysis 1
MaxMods = power(10,2);
MaxM = log10(MaxMods);

%Discretize the interval 1 to MaxM using NDL points for time analysis 2
Nmodz = round(logspace(0,MaxM,NDL));
if Nmodz(end) > 99999
    Nmodz(end) = 99999;
end

%SimPEG Npoints
Nmodz_SimPEG = Nmodz;

Nlays = round(logspace(0,2,NDL));

%% ESTIMATE TIME REQUIRED
prec1_time_est = prec1*0.5*Nmod_prec1;
prec2_time_est = prec2*(110+12*Nmod_prec2);
time1_time_est = time1*(sum(Nmodz_SimPEG)*UseSimPEG*(5+5+10+10+40+40)/1000+sum(Nmodz)*(12*25)*UseAarhusInv/1000+UseGAAEM*20*sum(Nmodz)/1000);
time2_time_est = time2*((2000*12+200+3*2500+800+100+1)/1000*(NSingleModTest*NDL)+NmodTest*NDL*(UseAarhusInv*12*90+UseSimPEG*(900+800+90+80+10+5)+UseGAAEM*100)/1000);

t_est = prec1_time_est+prec2_time_est+time1_time_est+time2_time_est;

ndays = t_est/(3600*24);
nhours = t_est/(3600);
nminutes = t_est/60;
nseconds = t_est;

if ndays < 1
    s = [num2str(nhours), ' hours'];
    if nhours < 1
        s = [num2str(nminutes), ' minutes'];
        if nminutes < 1
            s = [num2str(nseconds), ' seconds'];
        end
    end
else
    s = [num2str(ndays), ' days'];
end

disp(['Poor estimate of time required: ',s])
T0 = tic();
pause(2)

%% NAME THE RUN
STR1 = ['___P1_',num2str(Nmod_prec1)];
STR2 = ['___P2_',num2str(Nmod_prec2)];
STR3 = ['___E1_NMod_',num2str(NmodTest)];
STR4 = ['___E2_MaxMods_',num2str(MaxMods)];

Name = date;
%% SETUP TDEM SYSTEMS 1
%Ideal System
NoFilterGex = 'NoFilterGex\IdealSystem.gex';

%Near Surface (tTEM system)
NoFilterGex2 = 'NoFilterGex\TEM_nofilter.gex';

%Airborne System
NoFilterGex3 = 'NoFilterGex\ATEM_nofilter.gex';

%% SETUP TDEM SYSTEMS 2
% First, use a step-off waveform and circular loop, the relevant .gex file is "CircStepOff.gex"
GexFiles{1} = [pwd,'\VariousGexFiles\DipoleStepOff.gex'];

% Then, use a step-off waveform and rectangular loop, the relevant .gex file is "RectStepOff.gex"
GexFiles{2} = [pwd,'\VariousGexFiles\CircStepOff.gex'];

% Then, use a user-defined waveform and rectangular loop, the relevant .gex file is "RectUserDef.gex"
GexFiles{3} = [pwd,'\VariousGexFiles\RectStepOff.gex'];

% Then, use a user-defined waveform and rectangular loop, the relevant .gex file is "RectUserDef.gex"
GexFiles{4} = [pwd,'\VariousGexFiles\DipoleUserDef.gex'];

% Then, use the user-defined tTEM waveform and circular loop, the relevant .gex file is "CircUserDef.gex"
GexFiles{5} = [pwd,'\VariousGexFiles\CircUserDef.gex'];

% Then, use a user-defined waveform and rectangular loop, the relevant .gex file is "RectUserDef.gex"
GexFiles{6} = [pwd,'\VariousGexFiles\RectUserDef.gex'];

NGex = numel(GexFiles);

if prec1 == 1
    a = tic();

    resies = power(10,linspace(0,3,Nmod_prec1));
    model = cell(1,Nmod_prec1);

    for i = 1:Nmod_prec1
        depths = 0;
        ress = resies(i);
        model{i} = [depths;ress];
    end

    if UseAarhusInv == 1
        [Aarhus, AarhusInvTime] = calculate_forward_1D(model,'AarhusInv',0,NoFilterGex);
        AarhusUseGates = find(Aarhus.HM.UseGates);
    end

    [Analytic, AnalyticTime] = calculate_forward_1D(model,'Analytic',0,NoFilterGex);

    if UseSimPEG == 1
        [SimPEG, SimPEGTime] = calculate_forward_1D(model,'SimPEG',0,NoFilterGex);
        SimPEGUseGates = find(SimPEG.HM.UseGates);
    end

    if UseGAAEM == 1
        [GAAEM, GAAEMTime] = calculate_forward_1D(model,'GAAEM',0,NoFilterGex);
        GAAEMUseGates = find(GAAEM.HM.UseGates);

        GAAEM.HM.dBdt = GAAEM.HM.dBdt/2;
        GAAEM.LM.dBdt = GAAEM.LM.dBdt/2;
    end
    %%
    c = parula;
    col = @(res) c(round(255/500*res)+1,:);
    gates = Analytic.HM.Gates;
    ngates = numel(gates);




    if UseAarhusInv == 1

        AarhusInv_residuals = NaN(ngates,Nmod_prec1);
        AarhusInv_relerrors = NaN(ngates,Nmod_prec1);

        for i = 1:Nmod_prec1
            cur_analytic = Analytic.HM.dBdt(i,:);
            cur_AarhusInv = Aarhus.HM.dBdt(i,:);

            cur_AarhusInv_residual = cur_AarhusInv-cur_analytic(AarhusUseGates);
            cur_AarhusInv_relerror = cur_AarhusInv_residual./cur_analytic(AarhusUseGates);

            AarhusInv_residuals(AarhusUseGates,i) = cur_AarhusInv_residual;
            AarhusInv_relerrors(AarhusUseGates,i) = cur_AarhusInv_relerror;

        end
    end

    if UseGAAEM == 1

        GAAEM_residuals = NaN(ngates,Nmod_prec1);
        GAAEM_relerrors = NaN(ngates,Nmod_prec1);

        for i = 1:Nmod_prec1
            cur_analytic = Analytic.HM.dBdt(i,:);
            cur_GAAEM = GAAEM.HM.dBdt(i,:);

            cur_GAAEM_residual = cur_GAAEM-cur_analytic(GAAEMUseGates);
            cur_GAAEM_relerror = cur_GAAEM_residual./cur_analytic(GAAEMUseGates);

            GAAEM_residuals(GAAEMUseGates,i) = cur_GAAEM_residual;
            GAAEM_relerrors(GAAEMUseGates,i) = cur_GAAEM_relerror;

        end
    end

    if UseSimPEG == 1

        SimPEG_residuals = NaN(ngates,Nmod_prec1);
        SimPEG_relerrors = NaN(ngates,Nmod_prec1);

        for i = 1:Nmod_prec1
            cur_analytic = Analytic.HM.dBdt(i,:);
            cur_SimPEG = SimPEG.HM.dBdt(i,:);

            cur_SimPEG_residual = cur_SimPEG-cur_analytic(SimPEGUseGates);
            cur_SimPEG_relerror = cur_SimPEG_residual./cur_analytic(SimPEGUseGates);

            SimPEG_residuals(SimPEGUseGates,i) = cur_SimPEG_residual;
            SimPEG_relerrors(SimPEGUseGates,i) = cur_SimPEG_relerror;
        end
    end

    %% Calculate a sliding mean on approximate diffusion distance:
    Nwindows = 200;

    D_approx1=1260*sqrt(resies.*gates);
    D_approx = D_approx1(:);
    [~,inds] = sort(D_approx);
    sorted_d = D_approx(inds);

    window_indices = round(linspace(1,max(inds),Nwindows+1));
    mean_d = zeros(1,Nwindows);

    if UseGAAEM == 1
        GR=GAAEM_relerrors(:);
        sorted_g = GR(inds);
        mean_g = zeros(1,Nwindows);

        for i = 1:Nwindows
            start_ind = window_indices(i);
            end_ind = window_indices(i+1)-1;

            cur_ds = sorted_d(start_ind:end_ind);
            cur_gs = sorted_g(start_ind:end_ind);

            cur_mean_d = mean(cur_ds);
            cur_mean_g = mean(cur_gs);

            mean_d(i) = cur_mean_d;
            mean_g(i) = cur_mean_g;
        end
    end

    if UseAarhusInv == 1
        AR=AarhusInv_relerrors(:);
        sorted_a = AR(inds);
        mean_a = zeros(1,Nwindows);

        for i = 1:Nwindows
            start_ind = window_indices(i);
            end_ind = window_indices(i+1)-1;

            cur_ds = sorted_d(start_ind:end_ind);
            cur_as = sorted_a(start_ind:end_ind);

            cur_mean_d = mean(cur_ds);
            cur_mean_a = mean(cur_as);

            mean_d(i) = cur_mean_d;
            mean_a(i) = cur_mean_a;
        end
    end

    if UseSimPEG == 1
        SR=SimPEG_relerrors(:);
        sorted_s = SR(inds);
        mean_s = zeros(1,Nwindows);


        for i = 1:Nwindows
            start_ind = window_indices(i);
            end_ind = window_indices(i+1)-1;

            cur_ds = sorted_d(start_ind:end_ind);
            cur_ss = sorted_s(start_ind:end_ind);

            cur_mean_d = mean(cur_ds);
            cur_mean_s = mean(cur_ss);

            mean_d(i) = cur_mean_d;
            mean_s(i) = cur_mean_s;
        end
    end

    %% PRODUCE FIGURE 1
    nrowperpanel = 4;
    skip_nrowperpanel = 1;

    ncol = 3;

    nrow = (nrowperpanel+skip_nrowperpanel)*3+1;

    figure();
    p = @(pl) 1+(pl-1)*(nrowperpanel+skip_nrowperpanel)*ncol+ncol*((1:nrowperpanel)-1);

    % AarhusInv Subplots
    subplot(nrow,ncol,p(1))
    if UseAarhusInv == 1
        semilogx(gates,100*(AarhusInv_relerrors(:,2:end)),'-r','handlevisibility','off')
        hold on
        L(1)=semilogx(gates,100*(AarhusInv_relerrors(:,1)),'-r','displayname','Errors');
        L(4)=semilogx(gates,mean(100*(AarhusInv_relerrors),2),'--k','LineWidth',2,'DisplayName','Mean Error');
    end
    xlim([gates(1),gates(end)])
    ylim([-4,4])

    ylabel('Error [%]')
    title('AarhusInv')
    xlabel('Gate Time [s]')
    grid on
    L(2)=plot([1E-6 1E-3],[3 3],'-k','LineWidth',2,'DisplayName','3% Relative Error');
    hold on
    plot([1E-6 1E-3],[-3 -3],'-k','LineWidth',2,'DisplayName','3% Relative Error','HandleVisibility','off')
    L(3)=plot([1E-6 1E-3],[0 0],'-k','LineWidth',1,'DisplayName','0% Relative Error');
    box on
    set(gca,'xtick',[1E-6,1E-5,1E-4,1E-3])

    subplot(nrow,ncol,p(2))
    if UseAarhusInv == 1
        semilogx(resies,(AarhusInv_relerrors*100),'-r','HandleVisibility','off')
        hold on
        loglog(resies,mean(100*(AarhusInv_relerrors)),'--k','LineWidth',2,'DisplayName','Mean Error','HandleVisibility','off')
    end
    plot([1 1E3],[3 3],'-k','LineWidth',2,'DisplayName','3% Relative Error','HandleVisibility','off')
    plot([1 1E3],[-3 -3],'-k','LineWidth',2,'DisplayName','3% Relative Error','HandleVisibility','off')
    plot([1 1E3],[0 0],'-k','LineWidth',1,'DisplayName','0% Relative Error','HandleVisibility','off')
    grid on
    ylim([-4,4])
    ylabel('Error [%]')
    xlabel('Resitivity [\Omega m]')
    box on

    subplot(nrow,ncol,p(3))
    if UseAarhusInv == 1
        semilogx(1260*sqrt(resies.*gates),100*(AarhusInv_relerrors),'-r','HandleVisibility','off')
        hold on
        grid on
        loglog(mean_d,100*mean_a,'--k','LineWidth',2)
    end
    xlabel('Approximate Diffusion Distance [m]')
    plot([1 1E3],[3 3],'-k','LineWidth',2,'DisplayName','3% Relative Error','HandleVisibility','off')
    plot([1 1E3],[-3 -3],'-k','LineWidth',2,'DisplayName','3% Relative Error','HandleVisibility','off')
    plot([1 1E3],[0 0],'-k','LineWidth',1,'DisplayName','3% Relative Error','HandleVisibility','off')
    ylim([-4,4])
    ylabel('Error [%]')
    box on

    %GA-AEM Subplots
    subplot(nrow,ncol,p(1)+1)
    if UseGAAEM == 1
        semilogx(gates,100*(GAAEM_relerrors(:,2:end)),'-r','handlevisibility','off')
        hold on
        L(1)=semilogx(gates,100*(GAAEM_relerrors(:,1)),'-r','displayname','Errors');
        L(4)=semilogx(gates,mean(100*(GAAEM_relerrors),2),'--k','LineWidth',2,'DisplayName','Mean Error');
    end
    xlim([gates(1),gates(end)])
    ylim([-4,4])
    title('GAAEM')
    grid on
    xlabel('Gate Time [s]')
    plot([1E-6 1E-3],[3 3],'-k','LineWidth',2,'DisplayName','3% Relative Error','HandleVisibility','off')
    plot([1E-6 1E-3],[-3 -3],'-k','LineWidth',2,'DisplayName','3% Relative Error','HandleVisibility','off')
    plot([1E-6 1E-3],[0 0],'-k','LineWidth',1,'DisplayName','0% Relative Error','HandleVisibility','off')
    box on
    set(gca,'xtick',[1E-6,1E-5,1E-4,1E-3],'ytick',[])

    subplot(nrow,ncol,p(2)+1)
    if UseGAAEM == 1
        semilogx(resies,(GAAEM_relerrors*100),'-r','HandleVisibility','off')
        hold on
        loglog(resies,mean(100*(GAAEM_relerrors)),'--k','LineWidth',2,'DisplayName','Mean Error','HandleVisibility','off')
    end
    plot([1 1E3],[3 3],'-k','LineWidth',2,'DisplayName','3% Relative Error','HandleVisibility','off')
    plot([1 1E3],[-3 -3],'-k','LineWidth',2,'DisplayName','3% Relative Error','HandleVisibility','off')
    plot([1 1E3],[0 0],'-k','LineWidth',1,'DisplayName','3% Relative Error','HandleVisibility','off')
    grid on
    ylim([-4,4])
    xlabel('Resitivity [\Omega m]')
    box on

    subplot(nrow,ncol,p(3)+1)
    if UseGAAEM == 1
        semilogx(1260*sqrt(resies.*gates),100*(GAAEM_relerrors),'-r','HandleVisibility','off')
        hold on
        grid on
        loglog(mean_d,100*mean_g,'--k','LineWidth',2)
    end
    xlabel('Approximate Diffusion Distance [m]')
    plot([1 1E3],[3 3],'-k','LineWidth',2,'DisplayName','3% Relative Error','HandleVisibility','off')
    plot([1 1E3],[-3 -3],'-k','LineWidth',2,'DisplayName','3% Relative Error','HandleVisibility','off')
    plot([1 1E3],[0 0],'-k','LineWidth',1,'DisplayName','3% Relative Error','HandleVisibility','off')
    ylim([-4,4])
    box on

    %SimPEG subplots
    subplot(nrow,ncol,p(1)+2)
    if UseSimPEG == 1
        semilogx(gates,100*(SimPEG_relerrors(:,2:end)),'-r','handlevisibility','off')
        hold on
        L(1)=semilogx(gates,100*(SimPEG_relerrors(:,1)),'-r','displayname','Errors');
        L(4)=semilogx(gates,mean(100*(SimPEG_relerrors),2),'--k','LineWidth',2,'DisplayName','Mean Error');
    end
    xlim([gates(1),gates(end)])
    ylim([-4,4])
    grid on
    xlabel('Gate Time [s]')
    title('SimPEG')
    plot([1E-6 1E-3],[3 3],'-k','LineWidth',2,'DisplayName','3% Relative Error','HandleVisibility','off')
    plot([1E-6 1E-3],[-3 -3],'-k','LineWidth',2,'DisplayName','3% Relative Error','HandleVisibility','off')
    plot([1E-6 1E-3],[0 0],'-k','LineWidth',1,'DisplayName','0% Relative Error','HandleVisibility','off')
    box on
    set(gca,'xtick',[1E-6,1E-5,1E-4,1E-3],'ytick',[])

    subplot(nrow,ncol,p(2)+2)
    if UseSimPEG == 1
        semilogx(resies,(SimPEG_relerrors*100),'-r','HandleVisibility','off')
        hold on
        loglog(resies,mean(100*(SimPEG_relerrors)),'--k','LineWidth',2,'DisplayName','Mean Error','HandleVisibility','off')
    end
    plot([1 1E3],[3 3],'-k','LineWidth',2,'DisplayName','3% Relative Error','HandleVisibility','off')
    plot([1 1E3],[-3 -3],'-k','LineWidth',2,'DisplayName','3% Relative Error','HandleVisibility','off')
    plot([1 1E3],[0 0],'-k','LineWidth',1,'DisplayName','3% Relative Error','HandleVisibility','off')
    grid on
    xlabel('Resitivity [\Omega m]')
    ylim([-4,4])
    box on

    subplot(nrow,ncol,p(3)+2)
    if UseSimPEG == 1
        semilogx(1260*sqrt(resies.*gates),100*(SimPEG_relerrors),'-r','HandleVisibility','off')
        hold on
        grid on
        loglog(mean_d,100*mean_s,'--k','LineWidth',2)
    end
    xlabel('Approximate Diffusion Distance [m]')
    plot([1 1E3],[3 3],'-k','LineWidth',2,'DisplayName','3% Relative Error','HandleVisibility','off')
    plot([1 1E3],[-3 -3],'-k','LineWidth',2,'DisplayName','3% Relative Error','HandleVisibility','off')
    plot([1 1E3],[0 0],'-k','LineWidth',1,'DisplayName','3% Relative Error','HandleVisibility','off')
    ylim([-4,4])
    box on

    % Create a tile on the right column to get its position
    ax = subplot(15,1,15,'Visible','off');
    axPos = ax.Position;
    delete(ax)

    % Construct a Legend with the data from the sub-plots
    hL = legend(L,'orientation','horizontal');
    % Move the legend to the position of the extra axes
    hL.Position = axPos;

    %% SAVE RESULTS
    Name = [Name,STR1];
    save(Name)
    T1 = toc(a);
    clear a

end

if prec2 == 1
    a = tic();
    %% All Relerrors between AarhusInv, GA-AEM and SimPEG, as a function of Nlayers

    %In order to do this analysis we need at least 2 forward functions!
    filter = (UseAarhusInv+UseSimPEG+UseGAAEM)>1;

    if filter
        Ntest = 3;
        nlayz = logspace(0,2,Ntest);
        model_rel = cell(Ntest,Nmod_prec2);

        for j = 1:Ntest
            curNlay = nlayz(j);

            for i = 1:Nmod_prec2
                depths = [0,linspace(1,100,curNlay-1)];
                ress = 5+0.999*power(10,3*rand(1,curNlay));

                model_rel{j,i} = [depths;ress];
            end
        end

        clear Aar Sim GAA rAG rGA rAS rSA rGS rSG
        clear Aar2 Sim2 GAA2 rAG2 rGA2 rAS2 rSA2 rGS2 rSG2

        if UseAarhusInv == 1
            for j = 1:Ntest
                [Aarhus2, ~] = calculate_forward_1D(model_rel(j,:),'AarhusInv',0,NoFilterGex2);
                Aar{j} = Aarhus2;
                gates = Aarhus2.HM.Gates;
                AarhusUseGates = find(Aarhus2.HM.UseGates);
            end
        end

        if UseGAAEM == 1
            for j = 1:Ntest
                [GAAEM2, ~] = calculate_forward_1D(model_rel(j,:),'GAAEM',0,NoFilterGex2);
                GAA{j} = GAAEM2;
                gates = GAAEM2.HM.Gates;
                GAAEMUseGates = find(GAAEM2.HM.UseGates);
            end
        end

        if UseSimPEG == 1
            for j = 1:Ntest
                [SimPEG2, ~] = calculate_forward_1D(model_rel(j,:),'SimPEG',0,NoFilterGex2);
                Sim{j} = SimPEG2;
                gates = SimPEG2.HM.Gates;
            end
        end

        gates1 = gates;

        % Plot Relative Errors between GAAEM and AarhusInv
        for j = 1:Ntest
            if UseAarhusInv == 1
                A = Aar{j};
            end

            if UseGAAEM == 1
                G = GAA{j};
            end

            if UseSimPEG == 1
                S = Sim{j};
            end

            if and(UseAarhusInv==1,UseGAAEM==1)
                %AarhusInv Relative to GAAEM
                cur = A;
                ref = G;
                [rel_AG,~,~] = get_relerror(ref.HM.dBdt,cur.HM.dBdt,ref.HM.UseGates,cur.HM.UseGates);
                rAG(:,:,j) = rel_AG;
            end

            if and(UseAarhusInv==1,UseSimPEG==1)
                %AarhusInv Relative to SimPEG
                cur = A;
                ref = S;
                [rel_AS,~,~] = get_relerror(ref.HM.dBdt,cur.HM.dBdt,ref.HM.UseGates,cur.HM.UseGates);
                rAS(  :,:,j) = rel_AS;
            end

            if and(UseAarhusInv==1,UseGAAEM==1)
                %GAAEM Relative to AarhusInv
                cur = G;
                ref = A;
                [rel_GA,~,~] = get_relerror(ref.HM.dBdt,cur.HM.dBdt,ref.HM.UseGates,cur.HM.UseGates);
                rGA(:,:,j) = rel_GA;
            end

            if and(UseSimPEG==1,UseGAAEM==1)
                %GAAEM Relative to SimPEG
                cur = G;
                ref = S;
                [rel_GS,~,~] = get_relerror(ref.HM.dBdt,cur.HM.dBdt,ref.HM.UseGates,cur.HM.UseGates);
                rGS(:,:,j) = rel_GS;
            end

            if and(UseAarhusInv==1,UseSimPEG==1)
                %SimPEG Relative to AarhusInv
                cur = S;
                ref = A;
                [rel_SA,~,~] = get_relerror(ref.HM.dBdt,cur.HM.dBdt,ref.HM.UseGates,cur.HM.UseGates);
                rSA(:,:,j) = rel_SA;
            end

            if and(UseSimPEG==1,UseGAAEM==1)
                %SimPEG Relative to GAAEM
                cur = S;
                ref = G;
                [rel_SG,~,~] = get_relerror(ref.HM.dBdt,cur.HM.dBdt,ref.HM.UseGates,cur.HM.UseGates);
                rSG(:,:,j) = rel_SG;
            end
        end

        for j = 1:Ntest
            if UseAarhusInv == 1
                [Aarhus2, ~] = calculate_forward_1D(model_rel(j,:),'AarhusInv',0,NoFilterGex3);
                Aar2{j} = Aarhus2;
                gates = Aarhus2.HM.Gates;
                AarhusUseGates = find(Aarhus2.HM.UseGates);
            end

            if UseSimPEG == 1
                [SimPEG2, ~] = calculate_forward_1D(model_rel(j,:),'SimPEG',0,NoFilterGex3);
                Sim2{j} = SimPEG2;
                gates = SimPEG2.HM.Gates;
                SimPEGUseGates = find(SimPEG2.HM.UseGates);
            end

            if UseGAAEM == 1
                [GAAEM2, ~] = calculate_forward_1D(model_rel(j,:),'GAAEM',0,NoFilterGex3);
                GAA2{j} = GAAEM2;
                gates = GAAEM2.HM.Gates;
                GAAEMUseGates = find(GAAEM2.HM.UseGates);
            end
        end

        gates2 = gates;
        ngates = numel(gates);

        for j = 1:Ntest
            if UseAarhusInv == 1
                A = Aar2{j};
            end

            if UseGAAEM == 1
                G = GAA2{j};
            end

            if UseSimPEG == 1
                S = Sim2{j};
            end

            if and(UseAarhusInv==1,UseGAAEM==1)
                %AarhusInv Relative to GAAEM
                cur = A;
                ref = G;
                [rel_AG,~,~] = get_relerror(ref.HM.dBdt,cur.HM.dBdt,ref.HM.UseGates,cur.HM.UseGates);
                rAG2(:,:,j) = rel_AG;
            end

            if and(UseAarhusInv==1,UseSimPEG==1)
                %AarhusInv Relative to SimPEG
                cur = A;
                ref = S;
                [rel_AS,~,~] = get_relerror(ref.HM.dBdt,cur.HM.dBdt,ref.HM.UseGates,cur.HM.UseGates);
                rAS2(  :,:,j) = rel_AS;
            end

            if and(UseAarhusInv==1,UseGAAEM==1)
                %GAAEM Relative to AarhusInv
                cur = G;
                ref = A;
                [rel_GA,~,~] = get_relerror(ref.HM.dBdt,cur.HM.dBdt,ref.HM.UseGates,cur.HM.UseGates);
                rGA2(:,:,j) = rel_GA;
            end

            if and(UseSimPEG==1,UseGAAEM==1)
                %GAAEM Relative to SimPEG
                cur = G;
                ref = S;
                [rel_GS,~,~] = get_relerror(ref.HM.dBdt,cur.HM.dBdt,ref.HM.UseGates,cur.HM.UseGates);
                rGS2(:,:,j) = rel_GS;
            end

            if and(UseAarhusInv==1,UseSimPEG==1)
                %SimPEG Relative to AarhusInv
                cur = S;
                ref = A;
                [rel_SA,~,~] = get_relerror(ref.HM.dBdt,cur.HM.dBdt,ref.HM.UseGates,cur.HM.UseGates);
                rSA2(:,:,j) = rel_SA;
            end

            if and(UseSimPEG==1,UseGAAEM==1)
                %SimPEG Relative to GAAEM
                cur = S;
                ref = G;
                [rel_SG,~,~] = get_relerror(ref.HM.dBdt,cur.HM.dBdt,ref.HM.UseGates,cur.HM.UseGates);
                rSG2(:,:,j) = rel_SG;
            end
        end


        %% PLOTTING FIGURE 2
        Nrows = 2;
        Ncols = 3;

        if UseGAAEM == 1
            in = logical(GAAEMUseGates);
        elseif UseSimPEG == 1
            in = logical(SimPEGUseGates);
        elseif UseAarhusInv == 1
            in = logical(AarhusUseGates);
        end



        figure()
        clear lin lin2 lin3 linestyle linestyl
        W = 100;
        LI = 3;

        linestyl{1} = '-r';
        linestyl{2} = '-b';
        linestyl{3} = '-g';

        linestyle{1} = '-r';
        linestyle{2} = '-b';
        linestyle{3} = '-g';

        H1=subplot(14,3,[1,4,7,10,13]);
        for j = 1:Ntest

            if and(UseGAAEM==1,UseAarhusInv==1)
                curplot = 100*squeeze((abs((rAG(:,:,j)))))';


                %Calculate 0.025,0.5 and 0.975 percentiles
                P025 = prctile(curplot,2.5);
                P500 = prctile(curplot,50.0);
                P975 = prctile(curplot,97.5);

                if j == 1
                    lin(2) = loglog(gates1,P025,linestyl{j},'displayname','95% Confidence Bound');
                    grid on
                    hold on
                    lin(1) = loglog(gates1,P500,linestyle{j},'LineWidth',2,'displayname',['Median RAD, ',num2str(nlayz(j)),' layer(s)']);

                elseif j == 2
                    lin2(2) = loglog(gates1,P025,linestyl{j},'displayname','95% Confidence Bound');
                    grid on
                    hold on
                    lin2(1) = loglog(gates1,P500,linestyle{j},'LineWidth',2,'displayname',['Median RAD, ',num2str(nlayz(j)),' layer(s)']);
                else
                    lin3(2) = loglog(gates1,P025,linestyl{j},'displayname','95% Confidence Bound');
                    grid on
                    hold on
                    lin3(1) = loglog(gates1,P500,linestyle{j},'LineWidth',2,'displayname',['Median RAD, ',num2str(nlayz(j)),' layer(s)']);
                end

                loglog(gates1,P975,linestyl{j})

                if j == 1
                    title('AarhusInv/GA-AEM')
                    lin0=plot([2.1E-6,1E-2],[LI,LI],'-k','DisplayName','3% Relative Absolute Difference','LineWidth',2);
                    xlim([2.1E-6,1E-2])
                    ylim([1E-2,W])
                    xlabel('Gate Time [s]')
                    set(gca,'xtick',[1E-5,1E-4,1E-3,1E-2],'yticklabel',{})
                end
            end
        end

        H2=subplot(14,3,[1,4,7,10,13]+1);
        for j = 1:Ntest
            if and(UseSimPEG==1,UseAarhusInv==1)
                curplot = 100*squeeze((abs((rAS(:,:,j)))))';

                %Calculate 0.025,0.5 and 0.975 percentiles
                P025 = prctile(curplot,2.5);
                P500 = prctile(curplot,50.0);
                P975 = prctile(curplot,97.5);

                if j == 1
                    lin(2) = loglog(gates1,P025,linestyl{j},'displayname','95% Confidence Bound');
                    grid on
                    hold on
                    lin(1) = loglog(gates1,P500,linestyle{j},'LineWidth',2,'displayname',['Median RAD, ',num2str(nlayz(j)),' layer(s)']);

                elseif j == 2
                    lin2(2) = loglog(gates1,P025,linestyl{j},'displayname','95% Confidence Bound');
                    grid on
                    hold on
                    lin2(1) = loglog(gates1,P500,linestyle{j},'LineWidth',2,'displayname',['Median RAD, ',num2str(nlayz(j)),' layer(s)']);
                else
                    lin3(2) = loglog(gates1,P025,linestyl{j},'displayname','95% Confidence Bound');
                    grid on
                    hold on
                    lin3(1) = loglog(gates1,P500,linestyle{j},'LineWidth',2,'displayname',['Median RAD, ',num2str(nlayz(j)),' layer(s)']);
                end

                loglog(gates1,P975,linestyl{j})
                if j == 1
                    title('AarhusInv/SimPEG')
                    lin0=plot([2.1E-6,1E-2],[LI,LI],'-k','DisplayName','3% Relative Absolute Difference','LineWidth',2);
                    xlim([2.1E-6,1E-2])
                    ylim([1E-2,W])
                    xlabel('Gate Time [s]')
                    set(gca,'xtick',[1E-5,1E-4,1E-3,1E-2],'yticklabel',{})
                end


            end
        end

        H3=subplot(14,3,[1,4,7,10,13]+2);

        for j = 1:Ntest
            if and(UseGAAEM==1,UseSimPEG==1)
                curplot = 100*squeeze((abs((rSG(:,:,j)))))';

                %Calculate 0.025,0.5 and 0.975 percentiles
                P025 = prctile(curplot,2.5);
                P500 = prctile(curplot,50.0);
                P975 = prctile(curplot,97.5);

                if j == 1
                    lin(2) = loglog(gates1,P025,linestyl{j},'displayname','95% Confidence Bound');
                    grid on
                    hold on
                    lin(1) = loglog(gates1,P500,linestyle{j},'LineWidth',2,'displayname',['Median RAD, ',num2str(nlayz(j)),' layer(s)']);

                elseif j == 2
                    lin2(2) = loglog(gates1,P025,linestyl{j},'displayname','95% Confidence Bound');
                    grid on
                    hold on
                    lin2(1) = loglog(gates1,P500,linestyle{j},'LineWidth',2,'displayname',['Median RAD, ',num2str(nlayz(j)),' layer(s)']);
                else
                    lin3(2) = loglog(gates1,P025,linestyl{j},'displayname','95% Confidence Bound');
                    grid on
                    hold on
                    lin3(1) = loglog(gates1,P500,linestyle{j},'LineWidth',2,'displayname',['Median RAD, ',num2str(nlayz(j)),' layer(s)']);
                end

                loglog(gates1,P975,linestyl{j})

                if j == 1
                    title('SimPEG/GA-AEM')
                    lin0=plot([2.1E-6,1E-2],[LI,LI],'-k','DisplayName','3% Relative Absolute Difference','LineWidth',2);
                    xlim([2.1E-6,1E-2])
                    ylim([1E-2,W])
                    xlabel('Gate Time [s]')
                    set(gca,'xtick',[1E-5,1E-4,1E-3,1E-2],'yticklabel',{})
                end
            end
        end

        H4=subplot(14,3,[1,4,7,10,13]+21);

        for j = 1:Ntest
            if and(UseGAAEM==1,UseAarhusInv==1)
                curplot = 100*squeeze((abs((rAG2(:,:,j)))))';

                %Calculate 0.025,0.5 and 0.975 percentiles
                P025 = prctile(curplot,2.5);
                P500 = prctile(curplot,50.0);
                P975 = prctile(curplot,97.5);

                if j == 1
                    lin(2) = loglog(gates2,P025,linestyl{j},'displayname','95% Confidence Bound');
                    grid on
                    hold on
                    lin(1) = loglog(gates2,P500,linestyle{j},'LineWidth',2,'displayname',['Median RAD, ',num2str(nlayz(j)),' layer(s)']);

                elseif j == 2
                    lin2(2) = loglog(gates2,P025,linestyl{j},'displayname','95% Confidence Bound');
                    grid on
                    hold on
                    lin2(1) = loglog(gates2,P500,linestyle{j},'LineWidth',2,'displayname',['Median RAD, ',num2str(nlayz(j)),' layer(s)']);
                else
                    lin3(2) = loglog(gates2,P025,linestyl{j},'displayname','95% Confidence Bound');
                    grid on
                    hold on
                    lin3(1) = loglog(gates2,P500,linestyle{j},'LineWidth',2,'displayname',['Median RAD, ',num2str(nlayz(j)),' layer(s)']);
                end

                loglog(gates2,P975,linestyl{j})

                if j == 1
                    lin0=plot([2.1E-6,1E-2],[LI,LI],'-k','DisplayName','3% Relative Absolute Difference','LineWidth',2);

                    xlim([2.1E-6,1E-2])
                    ylim([1E-2,W])
                    xlabel('Gate Time [s]')
                    set(gca,'xtick',[1E-5,1E-4,1E-3,1E-2],'yticklabel',{})
                end
            end
        end

        H5=subplot(14,3,[1,4,7,10,13]+22);
        for j = 1:Ntest
            if and(UseAarhusInv==1,UseSimPEG==1)
                curplot = 100*squeeze((abs((rAS2(:,:,j)))))';

                %Calculate 0.025,0.5 and 0.975 percentiles
                P025 = prctile(curplot,2.5);
                P500 = prctile(curplot,50.0);
                P975 = prctile(curplot,97.5);

                if j == 1
                    lin(2) = loglog(gates2,P025,linestyl{j},'displayname','95% Confidence Bound');
                    grid on
                    hold on
                    lin(1) = loglog(gates2,P500,linestyle{j},'LineWidth',2,'displayname',['Median RAD, ',num2str(nlayz(j)),' layer(s)']);

                elseif j == 2
                    lin2(2) = loglog(gates2,P025,linestyl{j},'displayname','95% Confidence Bound');
                    grid on
                    hold on
                    lin2(1) = loglog(gates2,P500,linestyle{j},'LineWidth',2,'displayname',['Median RAD, ',num2str(nlayz(j)),' layer(s)']);
                else
                    lin3(2) = loglog(gates2,P025,linestyl{j},'displayname','95% Confidence Bound');
                    grid on
                    hold on
                    lin3(1) = loglog(gates2,P500,linestyle{j},'LineWidth',2,'displayname',['Median RAD, ',num2str(nlayz(j)),' layer(s)']);
                end

                loglog(gates2,P975,linestyl{j})

                if j == 1
                    lin0=plot([2.1E-6,1E-2],[LI,LI],'-k','DisplayName','3% Relative Absolute Difference','LineWidth',2);

                    xlim([2.1E-6,1E-2])
                    ylim([1E-2,W])
                    xlabel('Gate Time [s]')
                    set(gca,'xtick',[1E-5,1E-4,1E-3,1E-2],'yticklabel',{})
                end
            end
        end

        H6=subplot(14,3,[1,4,7,10,13]+23);

        for j = 1:Ntest

            if and(UseSimPEG==1,UseGAAEM==1)
                curplot = 100*squeeze((abs((rSG2(:,:,j)))))';

                %Calculate 0.025,0.5 and 0.975 percentiles
                P025 = prctile(curplot,2.5);
                P500 = prctile(curplot,50.0);
                P975 = prctile(curplot,97.5);

                if j == 1
                    lin(2) = loglog(gates2,P025,linestyl{j},'displayname','95% Confidence Bound');
                    grid on
                    hold on
                    lin(1) = loglog(gates2,P500,linestyle{j},'LineWidth',2,'displayname',['Median RAD, ',num2str(nlayz(j)),' layer(s)']);

                elseif j == 2
                    lin2(2) = loglog(gates2,P025,linestyl{j},'displayname','95% Confidence Bound');
                    grid on
                    hold on
                    lin2(1) = loglog(gates2,P500,linestyle{j},'LineWidth',2,'displayname',['Median RAD, ',num2str(nlayz(j)),' layer(s)']);
                else
                    lin3(2) = loglog(gates2,P025,linestyl{j},'displayname','95% Confidence Bound');
                    grid on
                    hold on
                    lin3(1) = loglog(gates2,P500,linestyle{j},'LineWidth',2,'displayname',['Median RAD, ',num2str(nlayz(j)),' layer(s)']);
                end

                loglog(gates2,P975,linestyl{j})

                if j == 1
                    lin0=plot([2.1E-6,1E-2],[LI,LI],'-k','DisplayName','3% Relative Absolute Difference','LineWidth',2);

                    xlim([2.1E-6,1E-2])
                    ylim([1E-2,W])
                    xlabel('Gate Time [s]')
                    set(gca,'xtick',[1E-5,1E-4,1E-3,1E-2],'yticklabel',{})
                end
            end
        end

        % Create a tile on the right column to get its position
        ax1 = subplot(20,1,20,'Visible','off');
        axPos1 = ax1.Position;
        delete(ax1)

        % Construct a Legend with the data from the sub-plots
        hL = legend([lin lin0 lin2 lin0 lin3 lin0],'orientation','vertical');
        % Move the legend to the position of the extra axes
        hL.Position = axPos1;
        hL.NumColumns=3;

        %Position Subplots
        set(H1,'Position',[.12, .65, .22, .30])
        set(H4,'Position',[.12, .30, .22, .30])
        set(H3,'Position',[.66, .65, .22, .30])
        set(H6,'Position',[.66, .30, .22, .30])
        set(H2,'Position',[.39, .65, .22, .30])
        set(H5,'Position',[.39, .30, .22, .30])

        %Place labels
        t1 = annotation("textbox");
        t1.FontSize = 8;
        t1.String = "a)";
        t1.Position = [.135 .92 0.03 0.03];
        t1.EdgeColor = "None";
        t1.FontWeight = "bold";

        t2 = annotation("textbox");
        t2.FontSize = 8;
        t2.String = "b)";
        t2.Position = [.405 .92 0.03 0.03];
        t2.EdgeColor = "None";
        t2.FontWeight = "bold";

        t3 = annotation("textbox");
        t3.FontSize = 8;
        t3.String = "c)";
        t3.Position = [.675 .92 0.03 0.03];
        t3.EdgeColor = "None";
        t3.FontWeight = "bold";

        t4 = annotation("textbox");
        t4.FontSize = 8;
        t4.String = "d)";
        t4.Position = [.135 .57 0.03 0.03];
        t4.EdgeColor = "None";
        t4.FontWeight = "bold";

        t5 = annotation("textbox");
        t5.FontSize = 8;
        t5.String = "e)";
        t5.Position = [.405 .57 0.03 0.03];
        t5.EdgeColor = "None";
        t5.FontWeight = "bold";

        t6 = annotation("textbox");
        t6.FontSize = 8;
        t6.String = "f)";
        t6.Position = [.675 .57 0.03 0.03];
        t6.EdgeColor = "None";
        t6.FontWeight = "bold";
    
    %% SAVE RESULTS
    Name = [Name,STR2];
    save(Name)
    else
        disp('Skipping Precision Analysis (Figure 2) since at least 2 forward algorithms are needed.')
    end
    T2 = toc(a);
    clear a
end


%% Now do computation time analysis
NmodzRep = 1;
clear AT GT ST AT3 GT3 ST3
model2 = cell(NmodTest,NDL);
model3 = cell(1,NDL);
model4 = cell(1,NDL);

for j = 1:NDL
    curNlays = Nlays(j);
    for i = 1:NmodTest

        if curNlays < 1
            depths = 0;
        end
        depths = [0, sort(99*rand(1,curNlays-1))];
        extra_depths = linspace(0,(curNlays-1)*0.1,curNlays);
        depths = depths+extra_depths;

        ress = 0.1+0.999*power(10,3*rand(1,curNlays));

        model2{i,j} = [depths;ress];
    end
end

for j = 1:NDL
    curNmod = Nmodz(j);
    tempmod = cell(1,curNmod);

    for i = 1:curNmod

        curNlays = NLT;

        if curNlays < 1
            depths = 0;
        else
            depths = linspace(0,100,curNlays);
        end

        ress = 0.1+0.999*power(10,3*rand(1,curNlays));
        tempmod{i} = [depths;ress];
    end

    model3{j} = tempmod;
end

for j = 1:NDL
    curNmod = Nmodz_SimPEG(j);
    tempmod = cell(1,curNmod);

    for i = 1:curNmod

        curNlays = NLT;

        if curNlays < 1
            depths = 0;
        else
            depths = linspace(0,100,curNlays);
        end

        ress = 0.1+0.999*power(10,3*rand(1,curNlays));
        tempmod{i} = [depths;ress];
    end

    model4{j} = tempmod;
end

if time2 == 1 %THIS IS HOW FAR I MADE IT APRIL 24TH
    %% Run the test 4 times, varying Loop Geometry between circular and rectangular geometries and varying the waveform between a step-off waveform and the tTEM waveform
    %% Calculate TIME PER 1 MOD/LAY (1 at a time!!)
    a=tic();
    t = waitbar(0,'');
    AT3c = zeros(NDL,2,NGex);
    GT3c = zeros(NDL,2,NGex);
    ST3c = zeros(NDL,2,NGex);
    NmodTest2 = min(NmodTest,NSingleModTest);

    for k = 1:NGex
        curgex = GexFiles{k};

        for i = 1:NDL
            for j = 1:NmodTest2
                cur_progress = (k-1)/NGex;

                if UseAarhusInv == 1
                    waitbar(cur_progress,t,['Nlays, 1 Mod, ','[',num2str(cur_progress*100),'%] AarhusInv 2 (',num2str(i),'/',num2str(NDL),')'])
                    [A3, AarhusInvTimeI,AarhusInvTimeW] = calculate_forward_1D(model2(j,i)','AarhusInv',0,curgex);
                    AT3c(i,1,k) = AT3c(i,1,k)+AarhusInvTimeI/NmodTest2;
                    AT3c(i,2,k) = AT3c(i,2,k)+AarhusInvTimeW/(2*NmodTest2);
                end

                if UseSimPEG == 1
                    waitbar(cur_progress,t,['Nlays, 1 Mod, ','[',num2str(cur_progress*100),'%] SimPEGTime 2 (',num2str(i),'/',num2str(NDL),')'])
                    [S3, SimPEGTimeI, SimPEGTimeW] = calculate_forward_1D(model2(j,i)','SimPEG',0,curgex);
                    ST3c(i,1,k) = ST3c(i,1,k)+SimPEGTimeI/NmodTest2;
                    ST3c(i,2,k) = ST3c(i,1,k)+SimPEGTimeW/NmodTest2;
                end
                if UseGAAEM == 1
                    if k == 4
                        waitbar(cur_progress,t,['Nlays, 1 Mod, ','[',num2str(cur_progress*100),'%] GAAEM 2 (',num2str(i),'/',num2str(NDL),')'])
                        [G3, GAAEMTimeI, GAAEMTimeW] = calculate_forward_1D(model2(j,i)','GAAEM',0,curgex);
                        GT3c(i,1,k) = GT3c(i,1,k)+GAAEMTimeI/NmodTest2;
                        GT3c(i,2,k) = GT3c(i,1,k)+GAAEMTimeW/(2*NmodTest2);
                    end
                end
            end
        end
    end
    T3a = toc(a);
    clear a

    close(t)

    disp(['Time Spent: ',num2str(T3a/3600),' hours'])

    %% Calculate TIME PER MOD/LAY
    a=tic();
    t = waitbar(0,'');
    AT3 = zeros(NDL,2,NGex);
    GT3 = zeros(NDL,2,NGex);
    ST3 = zeros(NDL,2,NGex);

    for k = 1:NGex
        curgex = GexFiles{k};

        for j = 1:NmodzRep
            for i = 1:NDL
                cur_progress = (k-1)/NGex;

                if UseAarhusInv == 1
                    waitbar(cur_progress,t,['Nlays, Many Mods, ','[',num2str(cur_progress*100),'%] AarhusInv 2 (',num2str(i),'/',num2str(NDL),')'])
                    [A3, AarhusInvTimeI,AarhusInvTimeW] = calculate_forward_1D(model2(:,i)','AarhusInv',0,curgex);
                    AT3(i,1,k) = AarhusInvTimeI;
                    AT3(i,2,k) = AarhusInvTimeW/2;
                end

                if UseSimPEG == 1
                    waitbar(cur_progress,t,['Nlays, Many Mods, ','[',num2str(cur_progress*100),'%] SimPEGTime 2 (',num2str(i),'/',num2str(NDL),')'])
                    [S3, SimPEGTimeI, SimPEGTimeW] = calculate_forward_1D(model2(:,i)','SimPEG',0,curgex);
                    ST3(i,1,k) = SimPEGTimeI;
                    ST3(i,2,k) = SimPEGTimeW;
                end

                if UseGAAEM == 1
                    if k == 4
                        waitbar(cur_progress,t,['Nlays, Many Mods, ','[',num2str(cur_progress*100),'%] GAAEM 2 (',num2str(i),'/',num2str(NDL),')'])
                        [G3, GAAEMTimeI, GAAEMTimeW] = calculate_forward_1D(model2(:,i)','GAAEM',0,curgex);
                        GT3(i,1,k) = GAAEMTimeI;
                        GT3(i,2,k) = GAAEMTimeW/2;
                    end
                end
            end
        end
    end
    T3b = toc(a);
    clear a

    close(t)

    disp(['Time Spent: ',num2str(T3b/3600),' hours'])

    %% NLAYS: TIME PER MOD/LAY
    figure()
    inum = 0;
    for k = 1:NGex
        %k loops over systems
        switch k
            case 1
                linestyle = '.-g';
                linestyle2 = 'og';
                linestyle3 = '.--g';
                linestyle4 = 'xg';
                dpn = ['S4'];
            case 2
                linestyle = '.-k';
                linestyle2 = 'ok';
                linestyle3 = '.--k';
                linestyle4 = 'xk';
                dpn = ['S5'];
            case 3
                linestyle = '.-r';
                linestyle2 = 'or';
                linestyle3 = '.--r';
                linestyle4 = 'xr';
                dpn = ['S6'];
            case 4
                linestyle = '.-b';
                linestyle2 = 'ob';
                linestyle3 = '.--b';
                linestyle4 = 'xb';
                dpn = ['S7'];
            case 5
                linestyle = '.-c';
                linestyle2 = 'oc';
                linestyle3 = '.--c';
                linestyle4 = 'xc';
                dpn = ['S8'];
            case 6
                linestyle = '.-m';
                linestyle2 = 'om';
                linestyle3 = '.--m';
                linestyle4 = 'xm';
                dpn = ['S9'];
        end

        dpn2 = ['Implemented: ',dpn];
        dpn1 = ['Raw: ',dpn];

        dpn2 = [dpn];
        dpn1 = [dpn];

        dpn3 = dpn1;
        dpn4 = dpn2;

        H1=subplot(14,3,[1,4,7,10,13]);
        inum = inum+1;
        lines(inum)=loglog(Nlays,1E3*AT3(:,1,k)/NmodTest,[linestyle],'DisplayName',dpn1,'LineWidth',1,'MarkerSize',7);
        hold on
        grid on
        inum = inum+1;
        lines(inum)=loglog(Nlays,1E3*AT3(:,2,k)/NmodTest,[linestyle2],'DisplayName',dpn2,'LineWidth',1,'MarkerSize',5);
        title('AarhusInv')
        ylabel('Time per response [ms]')
        ylim([0.5,20000])
        set(gca,'Xtick',[1E0,1E1,1E2],'XtickLabel',{})
        set(gca,'Ytick',[1E0,1E1,1E2,1E3,1E4])
        set(gca,'XLim',[0.8,125])

        H2=subplot(14,3,[1,4,7,10,13]+21);
        loglog(Nlays,1E3*AT3c(:,1,k),[linestyle],'DisplayName',dpn3,'LineWidth',1,'MarkerSize',7);
        hold on
        grid on
        loglog(Nlays,1E3*AT3c(:,2,k),[linestyle2],'DisplayName',dpn4,'LineWidth',1,'MarkerSize',5);
        xlabel('Number of layers')
        ylabel('Time per response [ms]')
        ylim([0.5,20000])
        set(gca,'Xtick',[1,1E1,1E2])
        set(gca,'Ytick',[1E0,1E1,1E2,1E3,1E4])
        set(gca,'XLim',[0.8,125])

        H3=subplot(14,3,[1,4,7,10,13]+2);
        loglog(Nlays,1E3*ST3(:,1,k)/NmodTest,[linestyle],'DisplayName',dpn1,'LineWidth',1,'MarkerSize',7)
        hold on
        grid on
        loglog(Nlays',1E3*ST3(:,2,k)/NmodTest,[linestyle2],'DisplayName',dpn2,'LineWidth',1,'MarkerSize',5)
        title('SimPEG')
        ylim([0.5,20000])
        set(gca,'Xtick',[1,1E1,1E2])
        set(gca,'Ytick',[1E0,1E1,1E2,1E3,1E4])
        set(gca,'Yticklabel',{},'XtickLabel',{})
        set(gca,'XLim',[0.8,125])

        H4=subplot(14,3,[1,4,7,10,13]+23);
        loglog(Nlays,1E3*ST3c(:,1,k),[linestyle],'DisplayName',dpn3,'LineWidth',1,'MarkerSize',7);
        hold on
        grid on
        loglog(Nlays,1E3*ST3c(:,2,k),[linestyle2],'DisplayName',dpn4,'LineWidth',1,'MarkerSize',5);
        ylim([0.5,20000])
        xlabel('Number of layers')
        set(gca,'Xtick',[1,1E1,1E2])
        set(gca,'Ytick',[1E0,1E1,1E2,1E3,1E4])
        set(gca,'Yticklabel',{})
        set(gca,'XLim',[0.8,125])

        if k == 4
            H5=subplot(14,3,[1,4,7,10,13]+1);
            loglog(Nlays,1E3*GT3(:,1,k)/NmodTest,[linestyle],'DisplayName',dpn1,'LineWidth',1,'MarkerSize',7)
            hold on
            grid on
            loglog(Nlays,1E3*GT3(:,2,k)./NmodTest,[linestyle2],'DisplayName',dpn2,'LineWidth',1,'MarkerSize',5)
            title('GA-AEM')
            ylim([0.5,20000])
            set(gca,'Xtick',[1,1E1,1E2])
            set(gca,'Ytick',[1E0,1E1,1E2,1E3,1E4])
            set(gca,'Yticklabel',{},'XtickLabel',{})
            set(gca,'XLim',[0.8,125])

            H6=subplot(14,3,[1,4,7,10,13]+22);
            loglog(Nlays,1E3*GT3c(:,1,k),[linestyle],'DisplayName',dpn3,'LineWidth',1,'MarkerSize',5);
            hold on
            grid on
            loglog(Nlays,1E3*GT3c(:,2,k),[linestyle2],'DisplayName',dpn4,'LineWidth',1,'MarkerSize',5);
            xlabel('Number of layers')
            ylim([0.5,20000])
            set(gca,'Xtick',[1,1E1,1E2])
            set(gca,'Ytick',[1E0,1E1,1E2,1E3,1E4])
            set(gca,'Yticklabel',{})
            set(gca,'XLim',[0.8,125])
        end
    end

    % Create a tile on the right column to get its position
    ax = subplot(20,1,20,'Visible','off');
    axPos = ax.Position;
    axPos(3) = 0.74;
    axPos(2) = 0.05;
    delete(ax)

    % Construct a Legend with the data from the sub-plots
    hL = legend(lines(1:2:end),'orientation','vertical');
    % Move the legend to the position of the extra axes
    hL.Position = axPos;
    hL.NumColumns=6;

    %Position Subplots
    set(H1,'Position',[.12, .60, .22, .35])
    set(H2,'Position',[.12, .20, .22, .35])
    set(H3,'Position',[.66, .60, .22, .35])
    set(H4,'Position',[.66, .20, .22, .35])
    set(H5,'Position',[.39, .60, .22, .35])
    set(H6,'Position',[.39, .20, .22, .35])

    %Place labels
    t1 = annotation("textbox");
    t1.FontSize = 8;
    t1.String = "a)";
    t1.Position = [.12 .92 0.03 0.03];
    t1.EdgeColor = "None";
    t1.FontWeight = "bold";

    t2 = annotation("textbox");
    t2.FontSize = 8;
    t2.String = "b)";
    t2.Position = [.39 .92 0.03 0.03];
    t2.EdgeColor = "None";
    t2.FontWeight = "bold";

    t3 = annotation("textbox");
    t3.FontSize = 8;
    t3.String = "c)";
    t3.Position = [.66 .92 0.03 0.03];
    t3.EdgeColor = "None";
    t3.FontWeight = "bold";

    t4 = annotation("textbox");
    t4.FontSize = 8;
    t4.String = "d)";
    t4.Position = [.12 .52 0.03 0.03];
    t4.EdgeColor = "None";
    t4.FontWeight = "bold";

    t5 = annotation("textbox");
    t5.FontSize = 8;
    t5.String = "e)";
    t5.Position = [.39 .52 0.03 0.03];
    t5.EdgeColor = "None";
    t5.FontWeight = "bold";

    t6 = annotation("textbox");
    t6.FontSize = 8;
    t6.String = "f)";
    t6.Position = [.66 .52 0.03 0.03];
    t6.EdgeColor = "None";
    t6.FontWeight = "bold";
    %%
    % SAVE RESULTS
    Name = [Name,STR3];
    save(Name)

end

%% NMODS: CALCULATE TIME PER MOD
if time1 == 1
    a=tic();
    AT3b = zeros(NDL,2,NGex);
    GT3b = zeros(NDL,2,NGex);
    ST3b = zeros(NDL,2,NGex);
    t = waitbar(0,'');
    for k = 1:NGex
        curgex = GexFiles{k};

        for j = 1:NmodzRep
            for i = 1:NDL
                cur_progress = (k-1)/NGex;

                if UseAarhusInv == 1
                    waitbar(cur_progress,t,['Nmods, ','[',num2str(cur_progress*100),'%] AarhusInv 2 (',num2str(i),'/',num2str(NDL),')'])
                    [A3, AarhusInvTimeI, AarhusInvTimeW] = calculate_forward_1D(model3{i},'AarhusInv',0,curgex);
                    AT3b(i,1,k) = AarhusInvTimeI;
                    AT3b(i,2,k) = AarhusInvTimeW/2;
                end

                waitbar(cur_progress,t,['Nmods, ','[',num2str(cur_progress*100),'%] SimPEGTime 2 (',num2str(i),'/',num2str(NDL),')'])

                z=1;

                if UseSimPEG == 1
                    [S3, SimPEGTimeI,SimPEGTimeW] = calculate_forward_1D(model4{i},'SimPEG',0,curgex);
                    ST3b(i,1,k) = SimPEGTimeI;
                    ST3b(i,2,k) = SimPEGTimeW;
                end

                if UseGAAEM == 1
                    if k == 4
                        waitbar(cur_progress,t,['Nmods, ','[',num2str(cur_progress*100),'%] GAAEM 2 (',num2str(i),'/',num2str(NDL),')'])
                        [G3, GAAEMTimeI, GAAEMTimeW] = calculate_forward_1D(model3{i},'GAAEM',0,curgex);
                        GT3b(i,1,k) = GAAEMTimeI;
                        GT3b(i,2,k) = GAAEMTimeW/2;
                    end
                end
            end
        end

    end
    b3 = toc(a);
    close(t)

    disp(['Time Spent: ',num2str(b3/3600),' hours'])

    %% NMODS: TIME PER MOD / NMOD
    figure()

    for k = 1:NGex
        %k loops over systems
        switch k
            case 1
                linestyle = '.-g';
                linestyle2 = 'og';
                dpn = ['S4'];
            case 2
                linestyle = '.-k';
                linestyle2 = 'ok';
                dpn = ['S5'];
            case 3
                linestyle = '.-r';
                linestyle2 = 'or';
                dpn = ['S6'];
            case 4
                linestyle = '.-b';
                linestyle2 = 'ob';
                dpn = ['S7'];
            case 5
                linestyle = '.-c';
                linestyle2 = 'oc';
                dpn = ['S8'];
            case 6
                linestyle = '.-m';
                linestyle2 = 'om';
                dpn = ['S9'];
        end

        dpn1 = ['Raw: ',dpn];
        dpn2 = ['Implemented: ',dpn];

        dpn1 = dpn;
        dpn2 = dpn;

        H1=subplot(7,3,[1,4,7,10,13]);
        lines(2*(k-1)+1)=loglog(Nmodz',1E3*AT3b(:,1,k)./Nmodz',[linestyle],'DisplayName',dpn1,'LineWidth',1,'MarkerSize',7);
        hold on
        grid on
        lines(2*k)=loglog(Nmodz',1E3*AT3b(:,2,k)./Nmodz',[linestyle2],'DisplayName',dpn2,'LineWidth',1,'MarkerSize',5);
        title('AarhusInv')
        xlabel('Number of responses')
        ylabel('Average time per response [ms]')
        ylim([1,10000])
        set(gca,'Xtick',[1E1,1E2,1E3,1E4,1E5])
        set(gca,'Ytick',[1E-1,1E0,1E1,1E2,1E3,1E4])
        set(gca,'XLim',[0.5,max(Nmodz)*2])

        H2=subplot(7,3,[1,4,7,10,13]+2);
        loglog(Nmodz_SimPEG',1E3*ST3b(:,1,k)./Nmodz_SimPEG',[linestyle],'DisplayName',dpn1,'LineWidth',1,'MarkerSize',7)
        hold on
        grid on
        loglog(Nmodz_SimPEG',1E3*ST3b(:,2,k)./Nmodz_SimPEG',[linestyle2],'DisplayName',dpn2,'LineWidth',1,'MarkerSize',5)
        title('SimPEG')
        xlabel('Number of responses')
        ylim([1,10000])
        set(gca,'Xtick',[1E1,1E2,1E3,1E4,1E5])
        set(gca,'XLim',[0.5,max(Nmodz_SimPEG)*2])
        set(gca,'Yticklabel',{})

        if k == 4
            H3=subplot(7,3,[1,4,7,10,13]+1);
            loglog(Nmodz',1E3*GT3b(:,1,k)./Nmodz',[linestyle],'DisplayName',dpn,'LineWidth',1,'MarkerSize',7)
            hold on
            grid on
            loglog(Nmodz',1E3*GT3b(:,2,k)./Nmodz',[linestyle2],'DisplayName',dpn,'LineWidth',1,'MarkerSize',5)
            title('GA-AEM')
            xlabel('Number of responses')
            ylim([1,10000])
            set(gca,'Xtick',[1E1,1E2,1E3,1E4,1E5])
            set(gca,'Yticklabel',{})
            set(gca,'XLim',[0.5,max(Nmodz)*2])
        end
    end

    % Create a tile on the right column to get its position
    ax = subplot(20,1,20,'Visible','off');
    axPos = ax.Position;
    axPos(3) = 0.7675;
    axPos(2) = 0.1;
    delete(ax)

    %Position Subplots
    set(H1,'Position',[.12, .25, .22, .67])
    set(H2,'Position',[.66, .25, .22, .67])
    set(H3,'Position',[.39, .25, .22, .67])

    % Construct a Legend with the data from the sub-plots
    hL = legend(lines(1:2:end),'orientation','vertical');
    % Move the legend to the position of the extra axes
    hL.Position = axPos;
    hL.NumColumns=6;

    t1 = annotation("textbox");
    t1.FontSize = 8;
    t1.String = "a)";
    t1.Position = [.13 .88 0.03 0.03];
    t1.EdgeColor = "None";
    t1.FontWeight = "bold";

    t2 = annotation("textbox");
    t2.FontSize = 8;
    t2.String = "b)";
    t2.Position = [.40 .88 0.03 0.03];
    t2.EdgeColor = "None";
    t2.FontWeight = "bold";

    t3 = annotation("textbox");
    t3.FontSize = 8;
    t3.String = "c)";
    t3.Position = [.67 .88 0.03 0.03];
    t3.EdgeColor = "None";
    t3.FontWeight = "bold";
    %%
    % SAVE RESULTS
    Name = [Name,STR4];
    save(Name)

end

T = toc(T0);

ndays = T/(3600*24);
nhours = T/(3600);
nminutes = T/60;
nseconds = T;

if ndays < 1
    s2 = [num2str(nhours), ' hours'];
    if nhours < 1
        s2 = [num2str(nminutes), ' minutes'];
        if nminutes < 1
            s2 = [num2str(nseconds), ' seconds'];
        end
    end
else
    s2 = [num2str(ndays), ' days'];
end
disp(['Predicted Time: ',s])
disp(['Actual Time Consumed: ',s2])
disp('')
