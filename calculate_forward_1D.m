function [out,calc_time_int,calc_time] = calculate_forward_1D(models,type,do_log_transform,gexname)
%Function returns the TDEM forward response for a 1D layered earth.

if nargin < 3
    do_log_transform = 0;
end

%Read system file, which should be the only .gex file in the local
%directory
if nargin < 4
    S = read_gex();
else
    S = read_gex(gexname);
end

LoopArea = S.General.TxLoopArea;

System_Name = S.General.Description{1};

switch type
    case 'GAAEM'
        

        %GAAEM is set to high resolution
        Nhank1 = 280;
        Nfreq1 = 12;
        Ndig1 = 10.35E6;
        
        switch System_Name
            case 'IdealSystem'
   
                %GAAEM is set to high resolution
                Nhank1 = 3400;
                Nfreq1 = 30;
                Ndig1 = 10.35E6;

                %Shrink the Gates for idealized dBdt values
                ActualGates = S.General.GateArray;
                ShrinkFactor = 0.125+0.45*power(linspace(0,1,max(size(ActualGates))),0.05);
                
                ShrinkFactor(1) = 0.125;
                ShrinkFactor(2) = 0.525;
                ShrinkFactor(3) = 0.625;
                ShrinkFactor(4) = 0.55;
                ShrinkFactor(5:14) = 0.6;
                ShrinkFactor(15:25) = 0.565;
                % ShrinkFactor(30) = 0.7;

                ActualGateWidths = diff(ActualGates(:,2:3)');
                NewGateWidths = ShrinkFactor.*ActualGateWidths;
                
                S.General.TxCoilPosition1(3) = -0.01;
                S.General.TxCoilPosition1(3) = -0.01;

                NewGates = repmat(ActualGates(:,1),1,3);
                NewGates(:,2) = NewGates(:,2)-NewGateWidths'/2;
                NewGates(:,3) = NewGates(:,3)+NewGateWidths'/2;

                S.General.GateArray = NewGates;
                
        end

        WriteSTMfiles(S,System_Name,Nhank1,Nfreq1,Ndig1);
        B = tic;
        [fwr,gates,calc_time_int] = tTEM_forward_GAAEM(models,S,System_Name);
        calc_time = toc(B);
        
        LM = struct;
        LM.dBdt = fwr.LM';
        LM.Gates = gates{1};

        HM = struct;
        HM.dBdt = fwr.HM';
        HM.Gates = gates{2};

        NGatesLM = S.Channel1.NoGates;
        NGatesHM = S.Channel2.NoGates;
        NGates = max(NGatesLM,NGatesHM);
        
        SkipGatesLM = S.Channel1.RemoveInitialGates;
        SkipGatesHM = S.Channel2.RemoveInitialGates;
        
        UseGatesLM = zeros(1,NGates);
        UseGatesHM = zeros(1,NGates);

        UseGatesLM(1+SkipGatesLM:NGatesLM) = 1;
        UseGatesHM(1+SkipGatesHM:NGatesHM) = 1;

        LM.UseGates = UseGatesLM;
        HM.UseGates = UseGatesHM;

        out.HM = HM;
        out.LM = LM;
        

    case 'AarhusInv'
        D = tic;

        WriteTEMfiles(S,System_Name)
        [out,calc_time_int]=tTEM_forward_AarhusInv(models,LoopArea,System_Name);
        calc_time = toc(D);

        NGatesLM = S.Channel1.NoGates;
        NGatesHM = S.Channel2.NoGates;
        
        NGates = max(NGatesLM,NGatesHM);
        
        SkipGatesLM = S.Channel1.RemoveInitialGates;
        SkipGatesHM = S.Channel2.RemoveInitialGates;
        
        UseGatesLM = zeros(1,NGates);
        UseGatesHM = zeros(1,NGates);

        UseGatesLM(1+SkipGatesLM:NGatesLM) = 1;
        UseGatesHM(1+SkipGatesHM:NGatesHM) = 1;

        out.LM.UseGates = UseGatesLM;
        out.HM.UseGates = UseGatesHM;

        

    case 'SimPEG'

        C = tic;
        [out,calc_time_int,calc_time_imp] = SimPEG_forward(models,S,LoopArea);
        calc_time = toc(C);
        calc_time = calc_time_imp;

    case 'Analytic'
        C = tic;
        out = get_analytic_tdem_circLoop(models,S,LoopArea);
        calc_time_int = toc(C);
        calc_time = calc_time_int;
end

%% LOG TRANSFORM %%
if do_log_transform == 1
    
    out.LM.dBdt = real(log10((out.LM.dBdt)));
    out.HM.dBdt = real(log10((out.HM.dBdt)));

end
end