function stimulus = get_stimulus(row, data_path)

animal_number = row.animal;
animal_part = row.animal_part;
run_number = row.run;

animal_folder = sprintf('animal_%03d', animal_number);

stimulus_filename_bare = sprintf('run%03d.txt', run_number);
stimulus_filename = fullfile(data_path, animal_folder, stimulus_filename_bare);

fileID = fopen(stimulus_filename,'r');
formatSpec = '%s %f';
text_data = textscan(fileID,formatSpec,'HeaderLines',3,'Delimiter','\t');
fclose(fileID);


text_length = length(text_data{1});
odors = text_data{1}(2:2:text_length); % odor names
odor_times = cumsum(text_data{2}); % odor ON and odor OFF times
odor_ON = odor_times(1:2:text_length); % includes end time (last element)
odor_OFF = odor_times(2:2:text_length);

fps = 4.3187;
odor_frames_ON = round(fps.*odor_ON);
odor_frames_OFF = round(fps.*odor_OFF);

stimulus = cell(length(odors),3);
for i = 1:length(odors)
stimulus(i,:) = {odors{i}, odor_frames_ON(i), odor_frames_OFF(i)};
end
