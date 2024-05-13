function WriteSTMfiles(S,system_name,Nhank,Nfreq,Ndig)
format long

if nargin < 2

    NumAbsHM = 250;
    NumAbsLM = 150;
    NumFreqHM = 10;
    NumFreqLM = 6;
    DigitFreq = 4E6;
    system_name = S.General.Description{1};

elseif nargin < 3

    NumAbsHM = 250;
    NumAbsLM = 150;
    NumFreqHM = 10;
    NumFreqLM = 6;
    DigitFreq = 4E6;


elseif nargin < 4

    NumFreqHM = 10;
    NumFreqLM = 6;
    DigitFreq = 4E6;

elseif nargin < 5

    DigitFreq = 4E6;
    
else

    NumAbsHM = Nhank;
    NumAbsLM = Nhank;
    NumFreqHM = Nfreq;
    NumFreqLM = Nfreq;
    DigitFreq = Ndig;

end

windows = S.General.GateArray;

LastWin_LM = S.Channel1.NoGates;
LastWin_HM = S.Channel2.NoGates;

SkipWin_LM = S.Channel1.RemoveInitialGates;
SkipWin_HM = S.Channel2.RemoveInitialGates;

TimeShiftLM = S.Channel1.GateTimeShift;
TimeShiftHM = S.Channel2.GateTimeShift;

windows_LM = windows(SkipWin_LM+1:LastWin_LM,2:3)+TimeShiftLM;
windows_HM = windows(SkipWin_HM+1:LastWin_HM,2:3)+TimeShiftHM;

NWin_LM = size(windows_LM,1);
NWin_HM = size(windows_HM,1);

%% PREPARE WAVEFORMS
LMWF = S.General.WaveformLM;
HMWF = S.General.WaveformHM;

LMWFTime1 = LMWF(1,1);
LMWFTime2 = LMWF(end,1);

HMWFTime1 = HMWF(1,1);
HMWFTime2 = HMWF(end,1);

LMWF_Period = 1./S.Channel1.RepFreq;
HMWF_Period = 1./S.Channel2.RepFreq;

%Check if at least half waveform is defined
LMWF_ishalf = (LMWFTime2-LMWFTime1) >= LMWF_Period/2;
HMWF_ishalf = (HMWFTime2-HMWFTime1) >= HMWF_Period/2;

if ~LMWF_ishalf
    LMWF = [LMWF;LMWFTime1+LMWF_Period/2,LMWF(end)];
end

if ~HMWF_ishalf
    HMWF = [HMWF;HMWFTime1+HMWF_Period/2,HMWF(end)];
end
%%

stm_dir = [pwd,'\GAAEM\stmfiles\'];

%Make sure the output folder exists
if ~isfolder(stm_dir)
    mkdir(stm_dir);
end

LM_name = [stm_dir,system_name,'_LM.stm'];
HM_name = [stm_dir,system_name,'_HM.stm'];

%make new file
fID_LM = fopen(LM_name,'w');
fID_HM = fopen(HM_name,'w');

fprintf(fID_LM,'%s\n','System Begin');
fprintf(fID_LM,'\t %s \n',['Name = ',S.General.Description{1}]);
fprintf(fID_LM,'\t %s \n','Type = Time Domain');

fprintf(fID_LM,'\n');

fprintf(fID_LM,'\t %s \n','Transmitter Begin');
fprintf(fID_LM,'\t \t %s \n','NumberOfTurns = 1');
fprintf(fID_LM,'\t \t %s \n','PeakCurrent = 1');
fprintf(fID_LM,'\t \t %s \n','LoopArea = 1');
fprintf(fID_LM,'\t \t %s \n',['BaseFrequency = ', num2str(S.Channel1.RepFreq)]);
fprintf(fID_LM,'\t \t %s \n',['WaveformDigitisingFrequency = ',num2str(DigitFreq)]);
fprintf(fID_LM,'\t \t %s \n','WaveFormCurrent Begin');
fprintf(fID_LM,'\t \t \t %10e %10e \n',LMWF');
fprintf(fID_LM,'\t \t %s \n','WaveFormCurrent End');

fprintf(fID_LM,'\n');

fprintf(fID_LM,'\t %s \n','Transmitter End');

fprintf(fID_LM,'\n');

fprintf(fID_LM,'\t %s \n','Receiver Begin');

fprintf(fID_LM,'\t \t %s \n',['NumberOfWindows = ',num2str(NWin_LM)]);
fprintf(fID_LM,'\t \t %s \n','WindowWeightingScheme = AreaUnderCurve');
fprintf(fID_LM,'\t \t %s \n','WindowTimes Begin');
fprintf(fID_LM,'\t \t \t %10e %10e \n',windows_LM');
fprintf(fID_LM,'\t \t %s \n','WindowTimes End');

fprintf(fID_LM,'\n');

TiBFilt = S.Channel1.TiBLowPassFilter;

fprintf(fID_LM,'\t \t %s \n','LowPassFilter Begin');
fprintf(fID_LM,'\t \t \t %s \n',['CutOffFrequency = ',num2str(TiBFilt(2))]);
fprintf(fID_LM,'\t \t \t %s \n',['Order = ',num2str(TiBFilt(1))]);
fprintf(fID_LM,'\t \t %s \n','LowPassFilter End');

fprintf(fID_LM,'\n');

fprintf(fID_LM,'\t %s \n','Receiver End');

fprintf(fID_LM,'\n');

fprintf(fID_LM,'\t %s \n','ForwardModelling Begin');

fprintf(fID_LM,'\n');

fprintf(fID_LM,'\t \t %s \n',['ModellingLoopRadius = ', num2str(sqrt(S.General.TxLoopArea/pi))]);
fprintf(fID_LM,'\t \t %s \n','OutputType = dB/dt');
fprintf(fID_LM,'\t \t %s \n','SaveDiagnosticFiles = no');
fprintf(fID_LM,'\t \t %s \n','XOutputScaling = 0');
fprintf(fID_LM,'\t \t %s \n','YOutputScaling = 0');
fprintf(fID_LM,'\t \t %s \n','ZOutputScaling = 1');
fprintf(fID_LM,'\t \t %s \n','SecondaryFieldNormalisation = none');
fprintf(fID_LM,'\t \t %s \n',['FrequenciesPerDecade = ',num2str(NumFreqLM)]);
fprintf(fID_LM,'\t \t %s \n',['NumberOfAbsiccaInHankelTransformEvaluation = ',num2str(NumAbsLM)]);
fprintf(fID_LM,'\t %s \n','ForwardModelling End');

fprintf(fID_LM,'\n');

fprintf(fID_LM,'%s \n','System End');



fprintf(fID_HM,'%s\n','System Begin');
fprintf(fID_HM,'\t %s \n',['Name = ',S.General.Description{1}]);
fprintf(fID_HM,'\t %s \n','Type = Time Domain');

fprintf(fID_HM,'\n');

fprintf(fID_HM,'\t %s \n','Transmitter Begin');
fprintf(fID_HM,'\t \t %s \n','NumberOfTurns = 1');
fprintf(fID_HM,'\t \t %s \n','PeakCurrent = 1');
fprintf(fID_HM,'\t \t %s \n','LoopArea = 1');
fprintf(fID_HM,'\t \t %s \n',['BaseFrequency = ', num2str(S.Channel2.RepFreq)]);
fprintf(fID_HM,'\t \t %s \n',['WaveformDigitisingFrequency = ',num2str(DigitFreq)]);
fprintf(fID_HM,'\t \t %s \n','WaveFormCurrent Begin');
fprintf(fID_HM,'\t \t \t %10e %10e \n',HMWF');
fprintf(fID_HM,'\t \t %s \n','WaveFormCurrent End');

fprintf(fID_HM,'\n');

fprintf(fID_HM,'\t %s \n','Transmitter End');

fprintf(fID_HM,'\n');

fprintf(fID_HM,'\t %s \n','Receiver Begin');

fprintf(fID_HM,'\t \t %s \n',['NumberOfWindows = ',num2str(NWin_HM)]);
fprintf(fID_HM,'\t \t %s \n','WindowWeightingScheme = AreaUnderCurve');
fprintf(fID_HM,'\t \t %s \n','WindowTimes Begin');
fprintf(fID_HM,'\t \t \t %10e %10e \n',windows_HM');
fprintf(fID_HM,'\t \t %s \n','WindowTimes End');

fprintf(fID_HM,'\n');

TiBFilt = S.Channel2.TiBLowPassFilter;

fprintf(fID_HM,'\t \t %s \n','LowPassFilter Begin');
fprintf(fID_HM,'\t \t \t %s \n',['CutOffFrequency = ',num2str(TiBFilt(2))]);
fprintf(fID_HM,'\t \t \t %s \n',['Order = ',num2str(TiBFilt(1))]);
fprintf(fID_HM,'\t \t %s \n','LowPassFilter End');

fprintf(fID_HM,'\n');

fprintf(fID_HM,'\t %s \n','Receiver End');

fprintf(fID_HM,'\n');

fprintf(fID_HM,'\t %s \n','ForwardModelling Begin');

fprintf(fID_HM,'\n');

fprintf(fID_HM,'\t \t %s \n',['ModellingLoopRadius = ', num2str(sqrt(S.General.TxLoopArea/pi))]);
fprintf(fID_HM,'\t \t %s \n','OutputType = dB/dt');
fprintf(fID_HM,'\t \t %s \n','SaveDiagnosticFiles = no');
fprintf(fID_HM,'\t \t %s \n','XOutputScaling = 0');
fprintf(fID_HM,'\t \t %s \n','YOutputScaling = 0');
fprintf(fID_HM,'\t \t %s \n','ZOutputScaling = 1');
fprintf(fID_HM,'\t \t %s \n','SecondaryFieldNormalisation = none');
fprintf(fID_HM,'\t \t %s \n',['FrequenciesPerDecade = ',num2str(NumFreqHM)]);
fprintf(fID_HM,'\t \t %s \n',['NumberOfAbsiccaInHankelTransformEvaluation = ',num2str(NumAbsHM)]);
fprintf(fID_HM,'\t %s \n','ForwardModelling End');

fprintf(fID_HM,'\n');

fprintf(fID_HM,'%s \n','System End');


fclose(fID_HM);
fclose(fID_LM);

end