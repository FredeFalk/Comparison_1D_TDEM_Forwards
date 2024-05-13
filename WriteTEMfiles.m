function WriteTEMfiles(S,system_name)

if nargin < 2
    system_name = S.General.Description{1};
end

format long

LoopArray = S.General.TxArray;
NSegments = numel(LoopArray(:))/2;

WaveFormLM = S.General.WaveformLM;
WaveFormHM = S.General.WaveformHM;

NSubWavesLM = numel(WaveFormLM)/2-1;
NSubWavesHM = numel(WaveFormHM)/2-1;

WFLM = zeros(1,NSubWavesLM*4);
WFHM = zeros(1,NSubWavesHM*4);

%Rearrange waveform to piecewise subwave-format
for i = 1:NSubWavesLM
    CurTimes = WaveFormLM(i:i+1,1);
    CurCurrs = WaveFormLM(i:i+1,2);

    WFLM(((i-1)*4+1):(i*4)) = [CurTimes(:),CurCurrs(:)];
end

for i = 1:NSubWavesHM
    CurTimes = WaveFormHM(i:i+1,1);
    CurCurrs = WaveFormHM(i:i+1,2);

    WFHM(((i-1)*4+1):(i*4)) = [CurTimes(:),CurCurrs(:)];
end

windows = S.General.GateArray;

LastWin_LM = S.Channel1.NoGates;
LastWin_HM = S.Channel2.NoGates;

SkipWin_LM = S.Channel1.RemoveInitialGates;
SkipWin_HM = S.Channel2.RemoveInitialGates;

%% 
TimeShiftLM = S.Channel1.GateTimeShift;
TimeShiftHM = S.Channel2.GateTimeShift;

FrontGateLM = S.Channel1.FrontGateTime+TimeShiftLM+S.General.FrontGateDelay;
FrontGateHM = S.Channel2.FrontGateTime+TimeShiftHM+S.General.FrontGateDelay;

windows_LM = windows(SkipWin_LM+1:LastWin_LM,:)+TimeShiftLM;
windows_HM = windows(SkipWin_HM+1:LastWin_HM,:)+TimeShiftHM;

NWin_LM = size(windows_LM,1);
NWin_HM = size(windows_HM,1);

GateArrayLM = zeros(NWin_LM,7);
GateArrayHM = zeros(NWin_HM,7);

GateArrayLM(:,1) = windows_LM(:,1);
GateArrayLM(:,6:7) = windows_LM(:,2:3);
GateArrayLM(:,2) = 1234e-4;
GateArrayLM(:,3) = S.Channel1.UniformDataSTD;
GateArrayLM(:,4:5) = 1;

GateArrayHM(:,1) = windows_HM(:,1);
GateArrayHM(:,6:7) = windows_HM(:,2:3);
GateArrayHM(:,2) = 1234e-4;
GateArrayHM(:,3) = S.Channel2.UniformDataSTD;
GateArrayHM(:,4:5) = 1;

LoopType = S.General.LoopType;

Description = S.General.Description{1};
Description2 = S.General.Description2{1};
Description3 = S.General.Description3{1};

usestepoff = 0;
switch Description3
    case 'userdefined'
        usestepoff = 0;
    case 'stepoff'
        usestepoff = 1;
end

switch Description2
    case 'Circ'
        LoopType = 12;
end

GateFormat = '\t %10d %5d %4d %2d %2d %10d %10d \n';

switch Description
    case 'IdealSystem'
        GateArrayHM = GateArrayHM(:,1:5);
        GateArrayLM = GateArrayLM(:,1:5);
        GateFormat = '\t %10d %5d %4d %2d %2d \n';
end


LM_name = [pwd,'\AarhusInv\ForwardModellingFiles\',system_name,'_LM.tem'];
HM_name = [pwd,'\AarhusInv\ForwardModellingFiles\',system_name,'_HM.tem'];
%% START WRITING FILE
%make new file
fID_LM = fopen(LM_name,'w');
fID_HM = fopen(HM_name,'w');

%% Write LM file
%Line 1: text label
fprintf(fID_LM,'%s\n','LM tem file produced using "WriteTEMfiles.m" by Frederik Falk');

%Line 2: 2 integers, source type and receiver polarization
fprintf(fID_LM,'%s \n',[num2str(LoopType),' 3']);

if LoopType < 7
    %Dipole
    switch LoopType
        case 1
            disp('horinztal dipoles not supported yet!')
            %Horizontal Dipole (x-oriented)

            %Line 3: 6 reals, Tx and Rx
            fprintf(fID_LM,'%s \n',[num2str(S.General.TxCoilPosition1),' ',num2str(S.General.RxCoilPosition1)]);

        case 2
        case 3
            %Vertical Dipole

            %Line 3: 6 reals, Tx and Rx
            fprintf(fID_LM,'%s \n',[num2str(S.General.TxCoilPosition1),' ',num2str(S.General.RxCoilPosition1)]);

            %Line 4: Tx Loop Dimensions (Radius of Loop) (DUMMY)
            fprintf(fID_LM,'\n');

        case 4
        case 5

    end

elseif LoopType == 7
    %Square Loop
    disp('Square Loop not supported yet!')
    
elseif and(LoopType>=12,LoopType<=14)
    %Circular Loop

    %Line 3: 6 reals, Tx and Rx
    fprintf(fID_LM,'%s \n',[num2str(S.General.TxCoilPosition1),' ',num2str(S.General.RxCoilPosition1)]);
    
    %Line 4: Tx Loop Dimensions (Radius of Loop)
    fprintf(fID_LM,'%f \n',[sqrt(S.General.TxLoopArea/pi)]);

elseif or(LoopType==72,LoopType==73)
    %Segmented Loop

    %Line 3: 6 reals, Tx and Rx
    fprintf(fID_LM,'%s \n',[num2str(S.General.TxCoilPosition1),' ',num2str(S.General.RxCoilPosition1)]);

    %Line 4: Tx Loop Dimensions (Number og segments, area of Loop)
    fprintf(fID_LM,'\t %f %f \n',[NSegments,S.General.TxLoopArea]);
    fprintf(fID_LM,'\t %4e %4e \n',LoopArray');
end

%Line 5: Data Transforms (3 = dBdt, see em1dinv_manual.pdf for details)
fprintf(fID_LM,'%s \n','3 3 3');

%Line 6: Transmitter Waveform Type (0 = Step Function, 1 = impulse (Step Function time Derivative) 3 = user-defined)
%and number of waveforms (1).
if usestepoff == 1
    fprintf(fID_LM,'%s \n','1 1');

    %Line 10: number of filters, modelling of front gate and damping of primary
    %field
    fprintf(fID_LM,'%s \n','1 0 0');

    %Line 11-16: Filters:
    fprintf(fID_LM,'%s \n',['1 ',num2str(S.General.RxCoilLPFilter1)]);
    fprintf(fID_LM,'%s \n','0');
    fprintf(fID_LM,GateFormat,GateArrayLM');
else
    fprintf(fID_LM,'%s \n','3 1');

    %Line 7-9: Tx Waveform Defintions
    fprintf(fID_LM,'%s \n',[num2str(NSubWavesLM),' ',num2str(WFLM)]);

    %Line 10: number of filters, modelling of front gate and damping of primary
    %field
    fprintf(fID_LM,'%s \n','1 1 1');
    
    %Line 11-16: Filters:
    fprintf(fID_LM,'%s \n',['1 ',num2str(S.General.RxCoilLPFilter1)]);
    fprintf(fID_LM,'%s \n','0');
    fprintf(fID_LM,'%s \n',FrontGateLM);
    fprintf(fID_LM,'%s \n',['1 ',num2str(S.Channel1.TiBLowPassFilter)]);
    fprintf(fID_LM,'%s \n','0');
    fprintf(fID_LM,GateFormat,GateArrayLM');
end

%% Write HM file
%Line 1: text label
fprintf(fID_HM,'%s\n','HM tem file produced using "WriteTEMfiles.m" by Frederik Falk');

%Line 2: 2 integers, source type and receiver polarization
fprintf(fID_HM,'%s \n',[num2str(LoopType),' 3']);

if LoopType < 7

    %Dipole
    switch LoopType
        case 1
            disp('horinztal dipoles not supported yet!')
            %Horizontal Dipole (x-oriented)

            %Line 3: 6 reals, Tx and Rx
            fprintf(fID_HM,'%s \n',[num2str(S.General.TxCoilPosition1),' ',num2str(S.General.RxCoilPosition1)]);
            
            %Line 4: Tx Loop Dimensions (Radius of Loop) (DUMMY)
            fprintf(fID_HM,'%f \n',[sqrt(S.General.TxLoopArea/pi)]);
        case 2
        case 3
            %Vertical Dipole

            %Line 3: 6 reals, Tx and Rx
            fprintf(fID_HM,'%s \n',[num2str(S.General.TxCoilPosition1),' ',num2str(S.General.RxCoilPosition1)]);

            %Line 4: Tx Loop Dimensions (Radius of Loop) (DUMMY)
            fprintf(fID_HM,' \n');

        case 4
        case 5

    end

elseif LoopType == 7
    %Square Loop
    disp('Square Loop not supported yet!')
    
elseif and(LoopType>=12,LoopType<=14)
    %Circular Loop

    %Line 3: 6 reals, Tx and Rx
    fprintf(fID_HM,'%s \n',[num2str(S.General.TxCoilPosition1),' ',num2str(S.General.RxCoilPosition1)]);
    
    %Line 4: Tx Loop Dimensions (Radius of Loop)
    fprintf(fID_HM,'%f \n',[sqrt(S.General.TxLoopArea/pi)]);

elseif or(LoopType==72,LoopType==73)
    %Segmented Loop

    %Line 3: 6 reals, Tx and Rx
    fprintf(fID_HM,'%s \n',[num2str(S.General.TxCoilPosition1),' ',num2str(S.General.RxCoilPosition1)]);

    %Line 4: Tx Loop Dimensions (Number og segments, area of Loop)
    fprintf(fID_HM,'\t %f %f \n',[NSegments,S.General.TxLoopArea]);
    fprintf(fID_HM,'\t %4e %4e \n',LoopArray');
end

%Line 5: Data Transforms (3 = dBdt, see em1dinv_manual.pdf for details)
fprintf(fID_HM,'%s \n','3 3 3');

%Line 6: Transmitter Waveform Type (0 = Step Function, 3 = user-defined)
%and number of waveforms (1).
if usestepoff == 1
    fprintf(fID_HM,'%s \n','1 1');

    %Line 10: number of filters, modelling of front gate and damping of primary
    %field
    fprintf(fID_HM,'%s \n','1 0 0');

    %Line 11-16: Filters:
    fprintf(fID_HM,'%s \n',['1 ',num2str(S.General.RxCoilLPFilter1)]);
    fprintf(fID_HM,'%s \n','0');
    fprintf(fID_HM,GateFormat,GateArrayHM');
else
    fprintf(fID_HM,'%s \n','3 1');

    %Line 7-9: Tx Waveform Defintions
    fprintf(fID_HM,'%s \n',[num2str(NSubWavesHM),' ',num2str(WFHM)]);

    %Line 10: number of filters, modelling of front gate and damping of primary
    %field
    fprintf(fID_HM,'%s \n','1 1 1');

    %Line 11-16: Filters:
    fprintf(fID_HM,'%s \n',['1 ',num2str(S.General.RxCoilLPFilter1)]);
    fprintf(fID_HM,'%s \n','0');
    fprintf(fID_HM,'%s \n',FrontGateHM);
    fprintf(fID_HM,'%s \n',['1 ',num2str(S.Channel2.TiBLowPassFilter)]);
    fprintf(fID_HM,'%s \n','0');
    fprintf(fID_HM,GateFormat,GateArrayHM');
end


fclose(fID_HM);
fclose(fID_LM);
end