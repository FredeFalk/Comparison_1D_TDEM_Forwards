function [out,time_calc_int,time_calc_imp] = SimPEG_forward(models,S,LoopArea)
%Requires pyrunfile to be set up.

%NOT COMPLETE

waveform_LM = S.General.WaveformLM;
waveform_HM = S.General.WaveformHM;

timeshift_LM = S.Channel1.GateTimeShift;
timeshift_HM = S.Channel2.GateTimeShift;

%Gate Times with Timeshift Correction
gate_times = S.General.GateArray;
LM_gatetimes = gate_times+timeshift_LM;
HM_gatetimes = gate_times+timeshift_HM;

%Receiver
Receiver = S.General.RxCoilPosition1;
Receiver(3) = abs(Receiver(3));
Receiver(1) = Receiver(1);

%Transmitter Geometry
TxArr = S.General.TxArray;
TxArray = zeros(size(TxArr,1),size(TxArr,2)+1)+abs(S.General.TxCoilPosition1);
TxArray(:,1:2) = TxArr;

TxArray = [TxArray;TxArray(1,:)];

Radius = sqrt(LoopArea/pi);

%Prepare variables
gates_LM = LM_gatetimes(:,1);
gates_HM = HM_gatetimes(:,1);

NGates_LM = numel(gates_LM);
NGates_HM = numel(gates_HM);
NModels = numel(models);

fwr_LM = zeros(NModels,NGates_LM);
fwr_HM = zeros(NModels,NGates_HM);

condmods = zeros(NModels,2);
thickmods = zeros(NModels,2);
nlays = zeros(1,NModels);

for i = 1:NModels
    cur_model = models{i};
    curNlay = size(cur_model,2);
    nlays(i) = curNlay;
    
    curdepths = abs(cur_model(1,:));

    if numel(curdepths)==1
        curthickz = curdepths;
    else
        curthickz = diff(curdepths);
    end

    thickmods(i,1:numel(curthickz)) = curthickz;
    condmods(i,1:curNlay) = 1./cur_model(2,:);
end

Description2 = S.General.Description2{1};
Description3 = S.General.Description3{1};

switch Description2
    case 'IdealSystem'
        pyfile = "SimPEG\tTEM_forward_SimPEG_halfspace.py";
    case 'Circ'
        pyfile = "SimPEG\tTEM_forward_SimPEG_Circ.py";
        
        TxArray = mean(TxArray(1:end-1,:)); %set xyz to [0,0,0]
        Receiver = TxArray;
    case 'Dipole'
        pyfile = "SimPEG\tTEM_forward_SimPEG_Dipole.py";
    case 'tTEM42'
        pyfile = "SimPEG\tTEM_forward_SimPEG.py";
end

switch Description3
    case 'stepoff'
        UseStep = 1;
    case 'userdefined'
        UseStep = 0;
end
NModels
    disp(pyfile)
    [SimPEG_response_HM,time_calc_int,time_calc_imp] = pyrunfile(pyfile,["dpred","time_calc_int","time_calc_imp"],...
                            waveform_times = waveform_HM(:,1),...
                            waveform_current = waveform_HM(:,2),...
                            times = gates_HM,...
                            source = TxArray,...
                            receiver = Receiver,...
                            thicknesses = thickmods,...
                            conductivities = condmods,...
                            Nlayers = nlays,...
                            R = Radius,...
                            usestep = UseStep,...
                            nmodz=NModels);

    % SimPEG_response_LM = pyrunfile("SimPEG\tTEM_forward_SimPEG_halfspace.py","dpred",...
    %                         waveform_times = waveform_LM(:,1),...
    %                         waveform_current = waveform_LM(:,2),...
    %                         times = gates_LM,...
    %                         source = TxArray,...
    %                         receiver = Receiver,...
    %                         thicknesses = curthicks',...
    %                         conductivities = curconds',...
    %                         R = Radius);

    SimPEG_response_LM = SimPEG_response_HM;

fwr_LM = double(-SimPEG_response_LM)/LoopArea;
fwr_HM = double(-SimPEG_response_HM)/LoopArea;

%Out should be a struct with two sub-structs named HM and LM, respectively.
%Each sub-struct has two variables, a [n_fwr,n_gate] double called "dBdt"
%and a [n_gate,1] double with gatetimes called "Gates".
out = struct;
HM = struct;
LM = struct;

HM.dBdt = fwr_HM;
HM.Gates = gates_HM;

LM.dBdt = fwr_LM;
LM.Gates = gates_LM;

LM.UseGates = gates_LM*0+1;
HM.UseGates = gates_HM*0+1;

out.HM = HM;
out.LM = LM;

end