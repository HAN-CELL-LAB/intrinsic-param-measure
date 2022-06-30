function itx2mat(file_name, file_path, dest_path, overwrite_opt) 
% itx2mat: converting ITX files to MAT files 
% USAGE
%     itx2mat -> prompt to open and select 
%     itx2mat(file_name)
%     itx2mat({file_1, file_2})
%     itx2mat(file_name, file_path)
%     itx2mat({file_1, file_2}, file_path)
%     itx2mat(file_name, file_path, dest_path)
%     itx2mat(file_name, file_path, dest_path, 'Y'/'N'/nan)
% OUPUT 
%     Converted MAT files are in the same folder as the ITX files. The MAT
%     files will be a struct
%     + itx_file:   location of the original ITX file
%     + recordings: struct array containing each sweep, each struct
%                   contains these fields: 
%                   + id:   unique segment id
%                   + dt:   time sampling interval (second) 
%                   + Fs:   sampling rate (Hz)
%                   + data: vector containing recorded data 

if ~exist('file_name', 'var')
    [file_name, file_path] = uigetfile('.itx', 'MultiSelect', 'on'); 
    if ~ischar(file_name) && ~iscell(file_name)
        fprintf('No file was chosen. Nothing was read. \n');    
        return;
    end     
    
end 

if ~exist('dest_path', 'var')
    dest_path = file_path;
end

if isempty(dest_path)
   dest_path = file_path;
end

if ~exist('overwrite_opt', 'var') 
    overwrite_opt = nan;
end

if iscell(file_name)
    for fn = 1:length(file_name)
        itx2mat(file_name{fn}, file_path, dest_path, overwrite_opt); 
    end
    return; 
end

if ~contains(file_name, '.itx')
    error('Need and ITX file');
end

itx_file = fullfile(file_path, file_name); 

file_pref = regexprep(file_name, '.itx$', ''); 
file_save = fullfile(dest_path, [file_pref '.mat']); 

if exist(file_save, 'file')
    fprintf('The file "%s" has already been converted to MAT file at "%s"\n', itx_file, file_save); 
    if isnan(overwrite_opt)
        continue_or_not = input('Do you want to overwrite (Y/N)? ', 's');
    else
        continue_or_not = overwrite_opt;
    end
    if strcmpi(continue_or_not, 'N')
        fprintf('... not overwriting %s.\n', file_save);
        return; 
    end
end

fid = fopen(itx_file,'r'); 
text_from_file = textscan(fid, '%s', 'delimiter', '\n');
text_from_file = text_from_file{1};

begin_indices = find(strcmp(text_from_file, 'BEGIN')); 
end_indices = find(strcmp(text_from_file, 'END')); 
if length(begin_indices) ~= length(end_indices)
    error('BEGIN and END number mismatch'); 
end
consec_begins = diff(begin_indices);
segment_durs = end_indices - begin_indices; 
if any(segment_durs(1:end-1) >= consec_begins)
    error('No nested BEGIN-END pairs allowed'); 
end 
num_segments = length(begin_indices); 
recordings = cell(num_segments, 1); 

% add on to find properties/metadata that Igor stupidly put at the end 
begin_indices = [begin_indices; length(text_from_file)]; 

fprintf('Beginning extraction of %s ... \n', itx_file); 
cur_tic = tic; 

for segment_ind = 1:num_segments
    begin_ind = begin_indices(segment_ind)+1; 
    end_ind = end_indices(segment_ind)-1; 
    
    % read the name after waves
    % this only accounts for anticipation of WAVES 
    % WAVES sb20190830_4_CFs_1_11_20_1p1
    segment_id = regexprep(text_from_file{begin_ind-2}, 'WAVES ', ''); 
    
    % read the data 
    read_data = cellfun(@str2double, text_from_file(begin_ind:end_ind), 'uni', 1);

    % only need the time axis, based on these 2. need to double check 
    % X SetScale/P x, 0.000000000E+00,  5.000000000E-05,"s", sb20190830_4_CFs_1_11_2_1p1
    % X SetAxis bottom,  0.000000000E+00,  1.000000000E+00; SetAxis left, -5.119843750E-09,  5.119843750E-09
    next_begin = begin_indices(segment_ind+1); 
    props = text_from_file(end_ind:next_begin); 
    find_SetScale = split(props(contains(props, 'SetScale/P x')), {',',';'});
    find_SetAxis = split(props(contains(props, 'SetAxis')), {',',';'});
    dt = str2double(find_SetScale{3}); 
    
    if length(find_SetAxis) >= 3 
        T = str2double(find_SetAxis{3});
    end

    note_msg = ''; 
    if ~isnan(T)
        if abs(T/length(read_data) - dt) > eps
            note_msg = sprintf('ERROR::mismatch_time'); 
        end
    end
    
    if ~isempty(note_msg)        
        fprintf('\t !!! at segment %d -> %s \n', segment_ind, note_msg);
    end
    
    recordings{segment_ind} = struct(...
        'id', segment_id, ...
        'dt', dt, ...
        'Fs', 1/dt, ...
        'T', T, ...
        'data', read_data, ...
        'note_msg', note_msg);
    
    
end
fclose(fid); 
recordings = vertcat(recordings{:});

save(file_save, 'recordings', 'itx_file'); 

fprintf('\t ... done. Saved as %s. \n\t ... Elapsed time = %.2f seconds.\n', file_save, toc(cur_tic)); 


end