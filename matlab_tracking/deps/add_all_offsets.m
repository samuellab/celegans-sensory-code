load_status;

needs_offset = status(isnan(status.offset) & status.use,:);
filename = 'status.xlsx';
sheet = 1;

if ~exist('offset', 'var')
    offset = 0;
end

for i = 1:size(status, 1)
    row = status(i,:);
    if (~isnan(row.offset) || ~row.use)
        continue
    end
    i
    
    dead_frames = 2;
    
    animal_number = row.animal;
    animal_part = row.animal_part;
    id_run_number = row.id_run;
    id_flip = row.id_flip;
    run_number = row.run;
    frame_offsets = row.offset;
    frame_counts = row.frame_count;
    
    dataset_id = sprintf('animal_%03d_%s', animal_number, animal_part);

    annotation_file_w = 'D:\workspace\NeuroPAL\annotator\annotations.json';

    original_jpegs_w = sprintf(...
        'D:\\workspace\\NeuroPAL\\annotator\\images\\%s', dataset_id);

    data_location_w = 'D:\data\180503 23pentanedione_butanone_NaCl';
    data_location_l = '\\10.245.89.254\Albert\NeuroPAL Ultra Data';

    annotation_file = annotation_file_w;
    original_jpegs = original_jpegs_w;
    data_location = data_location_w;


    %
    if animal_part == "tail"
        id_run = sprintf('id_tail_%03d.nd2', id_run_number);
    elseif animal_part == "head"
        id_run = sprintf('id%03d.nd2', id_run_number);
    else
        error('Unknown part: %s', animal_part);
    end   

    if animal_part == "tail"
        run = sprintf('run_tail_%03d.nd2', run_number);
    elseif animal_part == "head"
        run = sprintf('run%03d.nd2', run_number);
    else
        error('Unknown part: %s', animal_part);
    end     
    


    data_root = fullfile(data_location, sprintf('animal_%03d', animal_number));

    output_folder = fullfile(data_root, dataset_id);
    make_directory(output_folder);

    annotator_folder = fullfile(data_location, 'annotator');
    make_directory(annotator_folder);

    annotator_image_folder = fullfile(annotator_folder, 'images');
    make_directory(annotator_image_folder);

    id_shape = [256 512 21];

    shape = [128 256 23 floor(frame_counts/23)];


    % Get frame 1 from the id run

    A = get_array_data(BioFormatsArray(fullfile(data_root, id_run)));

    gcamp = A(:,:,:,4,1);
    colors = A(:,:,:,2,1);

    gcamp_small = imresize(gcamp, [128 256]);
    colors_small = imresize(colors, [128 256]);


    f = fullfile(data_root, run);
    crop_dead = @(x) x(:,:,(1+dead_frames):end);
    crop = @(x) x(:, 1:225, :);
    get_clean_slice = @(x,t) crop_dead(get_slice(x,t));
    a = get_clean_slice(N,1);
    figure(6); 
    imshowpair(max_intensity_y(gcamp), max_intensity_y(crop(imresize(a, 2)*4)));
    figure(7);
    imshow(max_intensity_y(a)*10);
    figure(8);
    imshow(max_intensity_y(gcamp)*10);
    figure(9);
    imshow(max_intensity_z(a)*10);
    figure(10);
    imshow(max_intensity_z(gcamp)*10);
    
    
    scale_up = @(x) imresize(crop(x), 8);
    prep_z = @(x) scale_up(max_intensity_z(x));
    prep_y = @(x) scale_up(max_intensity_y(x));
    
    done = false;
    bad = false;
    while ~done
        N = ND2Array(f, shape, offset);
        a = get_clean_slice(N,1);
        
        figure(1); 
        imshowpair(prep_y(gcamp_small), prep_y(a));
        figure(2); 
        imshowpair(prep_y(colors_small), prep_y(a));
        figure(3); 
        imshowpair(prep_z(colors_small), prep_z(a));
        
        reply = input('j: lower, k: accept, l:higher: ', 's');
        
        if reply == 'j'
            offset = offset - 1;
        elseif reply == 'l'
            offset = offset + 1;
        elseif reply == 'k'
            done = true;
        elseif reply == 'x'
            bad = true;
            done = true;
        end
        
        offset = mod(offset, 23);
    end
    
    if bad
        xlRange = sprintf('N%d', i+1);
        xlswrite(filename,{'''needs additional offset'},sheet,xlRange);
        
        xlRange = sprintf('I%d', i+1);
        xlswrite(filename,false,sheet,xlRange);
        continue
    end
    
    xlRange = sprintf('E%d', i+1);
    xlswrite(filename,offset,sheet,xlRange)
    
end