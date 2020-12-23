clc; 

%% Check for necessary toolboxes
required_toolboxes = [ ...
    "Signal Processing Toolbox", ...
    "Curve Fitting Toolbox"];

ver_str = struct2array(ver); 

for i = 1:length(required_toolboxes)
    toolbox_ith = required_toolboxes(i);
    if ~contains(ver_str, toolbox_ith)
        warning('The "%s" needs to be installed. Please install it before moving on!', toolbox_ith);
    else 
        fprintf('%s satisfied.\n', toolbox_ith);
    end
end

%% Decide whether to save path 
answer = questdlg('Save the source matlab path?', ...
	'Save path option', ...
	'Yes', 'No', 'No');

if strcmpi(answer, 'Yes')
    savepath('matlab-src');
end