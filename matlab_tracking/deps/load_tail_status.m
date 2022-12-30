%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: D:\data\180503 23pentanedione_butanone_NaCl\status.xlsx
%    Worksheet: tails
%
% To extend the code for use with different selected data or a different
% spreadsheet, generate a function instead of a script.

% Auto-generated by MATLAB on 2018/09/14 17:13:29

%% Import the data
[~, ~, raw] = xlsread('D:\data\180503 23pentanedione_butanone_NaCl\status.xlsx','tails','A2:L103');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
stringVectors = string(raw(:,2));
stringVectors(ismissing(stringVectors)) = '';
raw = raw(:,[1,3,4,5,6,7,8,9,10,11,12]);

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Create table
status = table;

%% Allocate imported array to column variable names
status.animal = data(:,1);
status.animal_part = categorical(stringVectors(:,1));
status.id_run = data(:,2);
status.run = data(:,3);
status.offset = data(:,4);
status.offset_frames = data(:,5);
status.transferred = data(:,6);
status.id_flip = data(:,7);
status.use = data(:,8);
status.volumes_made = data(:,9);
status.annotations_done = data(:,10);
status.traces_done = data(:,11);

%% Clear temporary variables
clearvars data raw stringVectors R;