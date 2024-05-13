function [fwd,TM,calc_time_int] = tTEM_forward_GAAEM(models,S,System_Name)

%resistivities is a cell array of size (1,N) where N is the number of
%models. Each cell contains an array of size (2,m) where the number of
%model layers is m. row 1 contains depths and row 2 contains resistivities.

%depths has the same structure as resistivites but contains depth values.
%The depth is assigned to the top of the respective layer.


%Load the shared library
gatdaem1d_loadlibrary();

%Create a LM system object, get its handle, and some basic info
LM.stmfile = [pwd,'\GAAEM\stmfiles\',System_Name,'_LM.stm'];
LM.hS  = gatdaem1d_getsystemhandle(LM.stmfile);                       
LM.nw  = gatdaem1d_nwindows(LM.hS);
LM.wt  = gatdaem1d_windowtimes(LM.hS);
LM.wfm = gatdaem1d_waveform(LM.hS);

%Create a HM system object, get its handle, and some basic info
HM.stmfile = [pwd,'\GAAEM\stmfiles\',System_Name,'_HM.stm'];
HM.hS  = gatdaem1d_getsystemhandle(HM.stmfile);
HM.nw  = gatdaem1d_nwindows(HM.hS);
HM.wt  = gatdaem1d_windowtimes(HM.hS);
HM.wfm = gatdaem1d_waveform(HM.hS);

rx_x = S.General.RxCoilPosition1(1);
rx_y = S.General.RxCoilPosition1(2);
rx_z = -S.General.RxCoilPosition1(3);

tx_x = S.General.TxCoilPosition1(1);
tx_y = S.General.TxCoilPosition1(2);
tx_z = -S.General.TxCoilPosition1(3);

txrx_dz = rx_z-tx_z;
txrx_dy = rx_y-tx_y;
txrx_dx = rx_x-tx_x;

    %Setup geometry
    G.tx_height = tx_z;
    G.tx_roll   = 0;         G.tx_pitch  = 0;       G.tx_yaw    = 0;
    G.txrx_dx   = txrx_dx;   G.txrx_dy   = txrx_dy; G.txrx_dz   = txrx_dz;
    G.rx_roll   = 0;         G.rx_pitch  = 0;        G.rx_yaw    = 0;  

NMOD = numel(models);
L = cell(1,NMOD);
H = cell(1,NMOD);

%Compute responses (put in loop change G and E as required)
wb1 = waitbar(0,'Calculating Forward');
c = tic;
for k=1:NMOD

    if mod(k,1000) == 0
        cur_rate = toc(c)/k;
        timeleft = ((NMOD-k)*cur_rate)/60;
        waitbar(k/NMOD,wb1,['Calculating Forwards, estimated ',num2str(timeleft),' minutes left'])
    end

    curthick = abs(diff(models{k}(1,:)));
    curres = models{k}(2,:);

    %Setup earth
    E = [];
    E.conductivity = (1./curres);        
    E.thickness    = curthick;        
            
    OutL = gatdaem1d_fm_dlogc(LM.hS,G,E);
    L{k} = OutL;
end
calc_time_int = toc(c);
for k=1:NMOD

    if mod(k,1000) == 0
        cur_rate = toc(c)/k;
        timeleft = ((NMOD-k)*cur_rate)/60;
        waitbar(k/NMOD,wb1,['Calculating Forwards, estimated ',num2str(timeleft),' minutes left'])
    end

    curthick = abs(diff(models{k}(1,:)));
    curres = models{k}(2,:);

    %Setup earth
    E = [];
    E.conductivity = (1./curres);        
    E.thickness    = curthick;        
            
    OutH = gatdaem1d_fm_dlogc(HM.hS,G,E);
    H{k} = OutH;
end

for k = 1:NMOD
    if k == 1
        NLM = numel(L{k}.FM.SZ);
        NHM = numel(H{k}.FM.SZ);
        LMz = zeros(NLM,NMOD);
        HMz = zeros(NHM,NMOD);
    end

    LMz(:,k) = L{k}.FM.SZ;
    HMz(:,k) = H{k}.FM.SZ;    
end

fwd.LM = -LMz;
fwd.HM = -HMz;

TM{1}=LM.wt.centre;
TM{2}=HM.wt.centre;

%Finished modelling so delete the system objects and unload the dll
gatdaem1d_freesystemhandle(LM.hS);
gatdaem1d_freesystemhandle(HM.hS);
gatdaem1d_unloadlibrary();
close(wb1)
end