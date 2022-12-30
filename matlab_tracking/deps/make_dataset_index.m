function make_dataset_index(asset_directory)
% MAKE_DATASET_INDEX(asset_directory)
%
%  This creates an 'index.json' file in the specified directory with
%  metadata for all of the datasets found in that directory.

listing = dir(asset_directory);
listing = listing(3:end);

index = struct();

make_filename = @(path, view, view_idx, channel_idx, t_idx) fullfile(...
    path, sprintf('%s_%d_%d_%d.jpg', view, view_idx, channel_idx, t_idx));

for i = 1:length(listing)
    
    if ~listing(i).isdir
        continue
    end
    id = listing(i).name;
    
    path = fullfile(asset_directory, id);
    dataset_index_location = fullfile(path, 'index.json');
    
    if exist(dataset_index_location, 'file')
        
        D = load_json(dataset_index_location);
        
    else
        
        D = struct();
        
        MIPZ = imread(make_filename(path, 'MIPZ', 0, 0, 0));
        MIPX = imread(make_filename(path, 'MIPX', 0, 0, 0));

        D.shape_x = size(MIPZ, 2);
        D.shape_y = size(MIPZ, 1);
        D.shape_z = size(MIPX, 2);
        D.shape_c = 0;
        D.shape_t = 0;
        
        % Pixel sizes for display.
        D.pixel_size_x = 2;
        D.pixel_size_y = 2;
        D.pixel_size_z = 8;
        
        save_json(D, dataset_index_location);
        
    end
    
    index.(id).id = id;
    index.(id) = merge_struct(index.(id), D);
    
end

save_json(index, fullfile(asset_directory, 'datasets.json'));