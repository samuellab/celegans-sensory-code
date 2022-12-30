local_root = 'D:\data\180503 23pentanedione_butanone_NaCl';
remote_root = 'G:\My Drive\data\171103 neuropal microfluidics\180503 23pentanedione_butanone_NaCl';
animals = 18;

% Pull traces
for a = animals
    name_string = sprintf('animal_%03d', a);
    local_animal = fullfile(local_root, name_string);
    
    tail_dir = fullfile(local_animal, [name_string '_tail']);
    if exist(tail_dir,'dir')
        pull_traces(tail_dir);
    end
    
    head_dir = fullfile(local_animal, [name_string '_head']);
    if exist(head_dir,'dir')
        pull_traces(head_dir)
    end
end

% Transfer
for a = animals
    name_string = sprintf('animal_%03d', a);
    local_animal = fullfile(local_root, name_string);
    
    tail_dir = fullfile(local_animal, [name_string '_tail']);
    if exist(tail_dir,'dir')
        listing = dir(fullfile(tail_dir, '*.mat'));
        
        for i = 1:length(listing)
            if listing(i).bytes < 10e6 % Threshold for 
                source_file = fullfile(tail_dir, listing(i).name);
                target_path = fullfile(remote_root, 'traces', name_string, ...
                    [name_string '_tail']);
                make_directory(target_path);
                target_file = fullfile(target_path, listing(i).name);
                copyfile(source_file, target_file);
            else % MIPs only for volume
                source_file = fullfile(tail_dir, listing(i).name);
                target_path = fullfile(remote_root, 'traces', name_string, ...
                    [name_string '_tail']);
                make_directory(target_path);
                target_file = fullfile(target_path, listing(i).name);
                A = MatfileArray(source_file, 'data');
                B = LazyArray(A, @max_intensity_z);
                mfile = matfile(target_file);
                mfile.data = zeros(size(B), 'uint8');
                for t = 1:size(B, 3)
                    mfile.data(:,:,t) = get_slice(B, t);
                end
            end
        end
    end
    
    head_dir = fullfile(local_animal, [name_string '_head']);
    if exist(head_dir,'dir')
        listing = dir(fullfile(head_dir, '*.mat'));
        
        for i = 1:length(listing)
            if listing(i).bytes < 10e6
                source_file = fullfile(head_dir, listing(i).name);
                target_path = fullfile(remote_root, 'traces', name_string, ...
                    [name_string '_head']);
                make_directory(target_path);
                target_file = fullfile(target_path, listing(i).name);
                copyfile(source_file, target_file);
            end
        end
    end
end

% Transfer MIP movies.

