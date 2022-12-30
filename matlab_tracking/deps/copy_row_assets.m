function copy_row_assets(row, data_path, remote_path, use_mips)

if nargin < 4
    use_mips = true;
end

animal_number = row.animal;
animal_part = row.animal_part;

animal_folder = sprintf('animal_%03d', animal_number);
dataset_id = sprintf('animal_%03d_%s', animal_number, animal_part);

root_folder = fullfile(data_path, animal_folder, dataset_id);
listing = dir(fullfile(root_folder, '*.mat'));

for i = 1:length(listing)
    source_file = fullfile(root_folder, listing(i).name);
    target_path = fullfile(remote_path, animal_folder, ...
        dataset_id);
    make_directory(target_path);
    target_file = fullfile(target_path, listing(i).name);
    if exist(target_file, 'file')
        continue
    elseif ~use_mips || listing(i).bytes < 20e6 % Threshold for non-image files
        copyfile(source_file, target_file);
    else % MIPs only for volume
        A = MatfileArray(source_file, 'data');
        B = LazyArray(A, @max_intensity_z);
        mfile = matfile(target_file);
        mfile.data = zeros(size(B), 'uint8');
        for t = 1:size(B, 3)
            mfile.data(:,:,t) = get_slice(B, t);
        end
    end
end

