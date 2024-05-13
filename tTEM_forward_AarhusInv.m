function [out,calc_time] = tTEM_forward_AarhusInv(models,LoopArea,system_name)
%Chunksize is for splitting the process into parts
Nmod = numel(models);
ChunkSize = Nmod;
Nmod = numel(models);

NChunk = ceil(Nmod/ChunkSize);

if nargin < 3
    system_name = 'tTEM42';
end

%RES is a M-by-N 2D double array with resistivites
% N is the number of soundings
% M is the number of resistivities/layers in a single sounding

%THICK is a (M-1)-by-N 2D double array with thicknesses
% N is the number of soundings
% M is the number of layers in a single sounding. The last layer is
% infinitely thick, which is why only M-1 are required.

% f = waitbar(0,'AarhusInv');

N = Nmod;

dBdt_HM = zeros(N,100);
dBdt_LM = zeros(N,100);

for i = 1:NChunk

ChunkStart = (i-1)*ChunkSize+1;

if i == NChunk
    curmodels = models(ChunkStart:end);
else
    ChunkEnd = i*ChunkSize;
    curmodels = models(ChunkStart:ChunkEnd);
end

CurChunkSize = numel(curmodels);

%Start by writing the model file
write_modfile(curmodels,'LM',system_name);
write_modfile(curmodels,'HM',system_name);


command_LM = string([pwd,'\AarhusInv\AarhusInv64.exe ' pwd,'\AarhusInv\ForwardModellingFiles\',system_name,'_LM.mod ' pwd,'\AarhusInv\Aarhusinv.con']);
command_HM = string([pwd,'\AarhusInv\AarhusInv64.exe ' pwd,'\AarhusInv\ForwardModellingFiles\',system_name,'_HM.mod ' pwd,'\AarhusInv\Aarhusinv.con']);

%Run Forward Response on the model file
a = tic(); 
system(command_LM)
calc_time = toc(a); %Take the internal calculation time
system(command_HM)


for j = 1:CurChunkSize
    %Read FWR file
    % ex is the name's leading zeros

    ex = '00000';
    ex = ex(1:end-length(num2str(j)));

    %Define j'th .fwr filename'
    filename_LM = [pwd,'\AarhusInv\ForwardModellingFiles\',system_name,'_LM',ex,num2str(j),'.fwr'];
    filename_HM = [pwd,'\AarhusInv\ForwardModellingFiles\',system_name,'_HM',ex,num2str(j),'.fwr'];

    %Read all lines
    data_LM_char = readlines(filename_LM);
    data_HM_char = readlines(filename_HM);

    %Find data start line
    if j == 1
        NumCharLastLine_LM = numel(cell2mat(data_LM_char(end-1)'));
        NumCharLastLine_HM = numel(cell2mat(data_HM_char(end-1)'));

        curline = 0;
        linefound = 0;
        while linefound == 0
            curline = curline+1;
                NumCharCurLine = numel(cell2mat(data_HM_char(curline)'));
                is_same_num_char = NumCharCurLine == NumCharLastLine_HM;
                
                if is_same_num_char
                    linefound = 1;
                    Data_Startline_HM = curline;
                end
                
        end

        curline = 0;
        linefound = 0;
        while linefound == 0
            curline = curline+1;
                NumCharCurLine = numel(cell2mat(data_LM_char(curline)'));
                is_same_num_char = NumCharCurLine == NumCharLastLine_LM;

                if is_same_num_char
                    linefound = 1;
                    Data_Startline_LM = curline;
                end
                
        end
    end

    %Take out data
    data_LM_mat = cell2mat(data_LM_char(Data_Startline_LM:end)');
    data_HM_mat = cell2mat(data_HM_char(Data_Startline_HM:end)');

    %Data contains 7 columns
    switch system_name
        case 'tTEM42'
            data_LM = reshape(str2double(strsplit(data_LM_mat(2:end),' ')),7,[])';
            data_HM = reshape(str2double(strsplit(data_HM_mat(2:end),' ')),7,[])';
        case 'IdealSystem'
            data_LM = reshape(str2double(strsplit(data_LM_mat(2:end),' ')),5,[])';
            data_HM = reshape(str2double(strsplit(data_HM_mat(2:end),' ')),5,[])';
    end
    
    if j == 1
        NG_LM = size(data_LM,1);
        NG_HM = size(data_HM,1);

        gates_LM = data_LM(:,1);
        gates_HM = data_HM(:,1);
    end

    dBdt_LM((i-1)*ChunkSize+j,1:NG_LM) = data_LM(:,2);
    dBdt_HM((i-1)*ChunkSize+j,1:NG_HM) = data_HM(:,2);
    
    delete(filename_LM)
    delete(filename_HM)
end


%OUTPUT: dBdt array has dimensions G-by-N where G is the number of gates
dBdt_LM = dBdt_LM(:,1:NG_LM);
dBdt_HM = dBdt_HM(:,1:NG_HM);

out.HM.dBdt = dBdt_HM/LoopArea;
out.HM.Gates = gates_HM;
out.LM.dBdt = dBdt_LM/LoopArea;
out.LM.Gates = gates_LM;

end