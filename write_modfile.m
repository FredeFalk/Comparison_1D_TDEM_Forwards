function write_modfile(models,type,system_name)

%RES is a M-by-N 2D double array with resistivites
% N is the number of soundings
% M is the number of resistivities/layers in a single sounding

%THICK is a (M-1)-by-N 2D double array with thicknesses
% N is the number of soundings
% M is the number of layers in a single sounding. The last layer is
% infinitely thick, which is why only M-1 are required.

% Set up modfile folder
if ~isfolder([pwd,'\ModFiles'])
    mkdir([pwd,'\ModFiles'])
end

N = numel(models);

% %Number of layer interfaces
switch type
    case 'HM'
        fileID = fopen([pwd,'\AarhusInv\ForwardModellingFiles\',system_name,'_HM.mod'],'w');

        %Setting up file head
        fprintf(fileID,'%6s\n','Model filer, Forward tTEM by Frederik Falk');
        fprintf(fileID,'%3.0f %3.0f\n',N,0);
        
        for i = 1:N
            fprintf(fileID,'%3.0f %3.0f %18s\n',i,1,[system_name,'_HM.tem']);
        end

    case 'LM'
        fileID = fopen([pwd,'\AarhusInv\ForwardModellingFiles\',system_name,'_LM.mod'],'w');

        %Setting up file head
        fprintf(fileID,'%6s\n','Model filer, Forward tTEM by Frederik Falk');
        fprintf(fileID,'%3.0f %3.0f\n',N,0);
        
        for i = 1:N
            fprintf(fileID,'%3.0f %3.0f %18s\n',i,1,[system_name,'_LM.tem']);
        end
end

fprintf(fileID,'%0.0f\n',-1);

%add N models
    for i = 1:N
        curmod = models{i};

        curdepths = abs(curmod(1,:));
        curthicks = diff(curdepths);

        curres = curmod(2,:);
        
        %models have M resistivity values followed by M-1 thickness values
        %and M-1 depth values
        
        M = numel(curres);
        N_interface = M-1;

        fprintf(fileID,'%0.0f\n',M);

        %Write all M resistivities
        for j = 1:M
           fprintf(fileID,'%9.2f %4.0f\n',curres(j),-1);
        end
        
        %Write all M-1 thicknesses
        for j = 1:N_interface
            fprintf(fileID,'%9.2f %4.0f\n',curthicks(j),-1);
        end
        
        %Write all M-1 depths
        for j = 1:N_interface
            fprintf(fileID,'%9.2f %4.0f\n',curdepths(j+1),-1);
        end
                
    end
    
    fclose(fileID);
end
