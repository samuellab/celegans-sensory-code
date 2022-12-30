function CTX_times_stimuli(dataset,raw_data_root,output_root)
% extract time, stimuli from text files

% raw_data_root = 'D:\ZM10104 2019\20190125 oct 10s';
% dataset = 'run201';
% output_root = 'D:\Dropbox\AL Data NG\ZM10104 (Sensory)\S_003';

% stimulus

cd(raw_data_root)
fileID = fopen(sprintf(strcat(dataset,'.txt')),'r');
formatSpec = '%s %f';
text_data = textscan(fileID,formatSpec,'HeaderLines',3,'Delimiter','\t');
fclose(fileID);
%cd(data_folder)
text_length = length(text_data{1});
odors = text_data{1}(2:2:text_length); % odor names
odor_times = cumsum(text_data{2}); % odor ON and odor OFF times
odor_ON = odor_times(1:2:text_length); % includes end time (last element)
odor_OFF = odor_times(2:2:text_length);


stimulus = cell(length(odors),3);
for i = 1:length(odors)
stimulus(i,:) = {odors(i), odor_ON(i), odor_OFF(i)};
end

% times
cd(fullfile(raw_data_root,dataset,'unprocessed_TIFs'))
times = load_json('times.json');

save(fullfile(output_root,strcat(dataset,'metadata.mat')),'times','stimulus')