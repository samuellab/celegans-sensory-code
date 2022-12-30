%% Master analysis script for CTX Data

%% Step 0: Load metadata from Excel


CTX_load_excel;

animals = {''};
line = '';

for k = 1:length(animals)

% identify this animal's datasets
if strcmp(line,'ZM10104')
    row_nums = find(strcmp(sensory_status.animal,animals{k}));
    raw_data_root = cell2mat(sensory_status.raw_root(row_nums(1)));
    output_root = cell2mat(sensory_status.analyzed_root(row_nums(1)));
    idstack = sensory_status.idstack(row_nums(1));
    datasets = sensory_status.run(row_nums);
end

if strcmp(char(java.net.InetAddress.getLocalHost.getHostName),'SAMUELLAB21')
annotator_root = 'D:\Dropbox\Analysis Package CTX\annotator';
elseif strcmp(char(java.net.InetAddress.getLocalHost.getHostName),'CORE-SAMLAB')
annotator_root = 'E:\Dropbox\Analysis Package CTX\annotator';
end

if ~exist(fullfile(output_root,animals{k}))
    mkdir(fullfile(output_root,animals{k}))
end

cd(raw_data_root)

% Step 1: generate 4D GCaMP volumes, time file and stimulus

    % convert id ND2 to mat
    CTX_id_convert_ND2_to_TIF_mat(char(strcat(idstack(1),'.nd2')));
    % convert raw mat to uint8 mat
    CTX_id_uint8_volume(char(idstack(1)),raw_data_root,fullfile(output_root,animals{k}));
    % generate jpegs for each dataset from id stack
    CTX_id_jpg_generator(char(idstack(1)),animals{k},raw_data_root,fullfile(annotator_root,'images'));
    
for i = 1:length(datasets)
    % converts the ND2 into raw TIF, raw (double/uint16) mat, times.json
    cd(raw_data_root)
    CTX_convert_ND2_to_TIF_mat(char(strcat(datasets(i),'.nd2')));
    % convert raw mat to uint8 mat
    CTX_uint8_volume(char(datasets(i)),raw_data_root,fullfile(output_root,animals{k}));
    % generate jpegs for each run, populates both channels green
    CTX_dataset_jpg_generator(char(datasets(i)),i+1,animals{k},raw_data_root,fullfile(annotator_root,'images'));
    % scaffold with time and stimulus integrated
    CTX_times_stimuli(char(datasets(i)),raw_data_root,fullfile(output_root,animals{k}));
end

f = msgbox(['Step 1 completed for Animal' ' ' animals{k}]);

end
%% Step 2: annotate in annotator

% INSTRUCTIONS

% add animal to Datasets.json
% Launch Web Server for Chrome, folder: D:\Dropbox\Analysis Package CTX\annotator

%% Step 3: Track (tracks all datasets for a given animal)
CTX_load_excel;

animals = {''};
line = '';

for k = 1:length(animals)

% identify this animal's datasets
if strcmp(line,'ZM10104')
    row_nums = find(strcmp(sensory_status.animal,animals{k}));
    raw_data_root = cell2mat(sensory_status.raw_root(row_nums(1)));
    output_root = cell2mat(sensory_status.analyzed_root(row_nums(1)));
    idstack = sensory_status.idstack(row_nums(1));
    datasets = sensory_status.run(row_nums);
end

if strcmp(char(java.net.InetAddress.getLocalHost.getHostName),'SAMUELLAB21')
annotator_root = 'D:\Dropbox\Analysis Package CTX\annotator';
elseif strcmp(char(java.net.InetAddress.getLocalHost.getHostName),'samuel-lab-1')
annotator_root = 'I:\Dropbox\RV Data CTX\RV annotator';   
elseif strcmp(char(java.net.InetAddress.getLocalHost.getHostName),'CORE-SAMLAB')
annotator_root = 'E:\Dropbox\Analysis Package CTX\annotator';
end

if ~exist(fullfile(output_root,animals{k}))
    mkdir(fullfile(output_root,animals{k}))
end

cd(raw_data_root)

CTX_track_neurons(animals{k},line,annotator_root)

% Generate and plotting traces mat file (copy CTX versions at some point)
for i = 1:length(datasets)
   close all
   NG_generate_traces(char(datasets(i)),fullfile(output_root,animals{k}))
   CTX_plot_traces(char(datasets(i)),fullfile(output_root,animals{k}))
end
f = msgbox(['Step 3 completed for Animal' ' ' animals{k}]);

end

%% NOT functions, open these programs to run 
% need to check these programs

CTX_dataset_proofreader % once per dataset

CTX_avg_traces % once per animal