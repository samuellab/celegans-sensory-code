load_status;

data_location_l = '\\10.245.89.254\Albert\NeuroPAL Ultra Data';
data_location_w = 'D:\data\180503 23pentanedione_butanone_NaCl';

uncopied = status(~status.transferred & status.use,:);

for i = 1:size(uncopied, 1)
    tic
    row = uncopied(i,:)
    
    animal_string =  sprintf('animal_%03d', row.animal);
    animal_remote = fullfile(data_location_l,animal_string);
    animal_local = fullfile(data_location_w, animal_string);
    
    odor_string = sprintf('run%03d.txt', row.run);
    local_odors = fullfile(animal_local, odor_string);
    remote_odors = fullfile(animal_remote, odor_string);
    
    if row.animal_part == "tail"
        dataset_string = sprintf('run_tail_%03d.nd2', row.run);
    elseif row.animal_part == "head"
        dataset_string = sprintf('run%03d.nd2', row.run);
    end
    
    remote_run = fullfile(animal_remote, dataset_string);
    local_run = fullfile(animal_local, dataset_string);
    
    copyfile(remote_odors, local_odors);
        
    if ~exist(local_run, 'file')
        copyfile(remote_run, local_run);
    end
    toc
end