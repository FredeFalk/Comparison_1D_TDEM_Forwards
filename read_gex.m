function S = read_gex(name)
%This function gets information from a given gex file

%Find local system file if "name" is not provided
if nargin < 1
    local_dir = dir;
    Nfiles = numel(local_dir);
    
    %Search for a .gex file in local_dir
    gex_found = 0;

    for i = 1:Nfiles
        cur_file = local_dir(i);
        cur_name = cur_file.name;

        cur_name_length = numel(cur_name);
        
        if cur_name_length >= 3
            cur_file_extension = cur_file.name(end-3:end);
        
            gex_Check = cumprod(cur_file_extension == '.gex');
        
            if gex_Check == 1

                gex_found = 1;
                name = cur_file.name;

            end

        end
        
    end

    if gex_found == 0
        disp('No GEX file found in local directory')
    else
        disp(['Found GEX file in local directory: ',name])
    end
end

%Read File
file = readlines(name);

Number_of_Lines = numel(file);

for i = 1:Number_of_Lines
    curLine = char(file(i));

    %Skip if Line is empty
    if numel(curLine) > 0

        %If a square bracket is encountered make new struct field
        newsection = curLine(1) == '[';
        if newsection
            CurStruct = curLine(2:end-1);
        end

        %Check if the Line contains an '=' sign
        Check = logical(sum(curLine == '='));

        % if it does, split the line and enter information into current
        % field
        if Check
            %Split line at equal sign
            SplitLine = strsplit(curLine,'=');
            Right_Handside = strtrim(SplitLine{2});

            %Split right handside part by spaces to separate multiple
            %values
            SplitLine2 = strsplit(Right_Handside,' ');
            numericvalue = str2double(SplitLine2);

            %Determine if right hand side is numeric or string by looking
            %for NaN values after conversion to double
            right_handside_string = logical(sum(isnan(numericvalue)));

            switch right_handside_string
                case 1
                    fieldvalue = SplitLine2;
                case 0
                    fieldvalue = numericvalue;
            end

            %Lefthand side will be the field name
            fieldname = SplitLine{1};

            %Before assigning values to a field, check if the current field
            %is part of the waveform, loop geometry array or the gatetime array
            CheckWaveForm = logical(prod(fieldname(1:4) == 'Wave'));
            CheckGateArray = logical(prod(fieldname(1:4) == 'Gate'));

            CheckTxArray1 = logical(prod(fieldname(1:4) == 'TxLo'));
            CheckTxArray2 = logical(prod(fieldname(end-3:end-1)=='int'));
            CheckTxArray = and(CheckTxArray1,CheckTxArray2);

            %If any is true then store the information
            AnyArrayTrue = logical(sum([CheckWaveForm,CheckGateArray,CheckTxArray]));
            if AnyArrayTrue
                FirstRow_check = exist('temporary_array');
                
                if FirstRow_check
                    temporary_array = [temporary_array;fieldvalue];
                else
                    temporary_array = [];
                    temporary_array = [temporary_array;fieldvalue];
                end
            end
                S.(CurStruct).(fieldname) = fieldvalue;
            
        end
    else
        %Empty line! clear temporary array if it exists after placing it
        %inside field!

        isTempArray = exist('temporary_array');

        if isTempArray
            if CheckWaveForm
                S.(CurStruct).(fieldname(1:10)) = temporary_array;
            elseif CheckGateArray
                S.(CurStruct).GateArray = temporary_array;
            elseif CheckTxArray
                S.(CurStruct).TxArray = temporary_array;
            end
            clear temporary_array
        end
    end

end
end