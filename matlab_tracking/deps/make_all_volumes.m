load_status;

filename = 'status.xlsx';
sheet = 1;
dead_frames = 2;


for idx = 1:size(status, 1)
    row = status(idx,:);
    if (row.volumes_made || ~row.use)
        continue
    end
    idx
    row
    
    
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

    A = get_array_data(BioFormatsArray(fullfile(data_root, id_run)));

    gcamp = A(:,:,:,4,1);

    gcamp_small = imresize(gcamp, [128 256]);

    % Find a good range to make uint8

    convert_to_uint8 = {};

    quantile_target_lo = 0.9;
    quantile_target_hi = 0.9999;


    f = fullfile(data_root, run);
    N = ND2Array(f, shape, frame_offsets);
    crop_dead = @(x) x(:,:,(1+dead_frames):end);
    get_clean_slice = @(x,t) crop_dead(get_slice(x,t));

    size_T = size(N, 4);

    hi = 0;
    qs_lo = [];
    qs_hi = [];
    for i = 1:10

        t = randi(size_T);

        a = get_clean_slice(N, t);
        qs_lo(end+1) = quantile(a(:), quantile_target_lo);
        qs_hi(end+1) = quantile(a(:), quantile_target_hi);

    end

    lo = min(qs_lo);
    hi = max(qs_hi);
    convert_to_uint8 = @(x) uint8((double(x)-lo)/(hi-lo)*255);


    t = randi(size_T);

    a = get_clean_slice(N, t);

    az = max_intensity_z(a);
    sample_slice_for_jpeg = convert_to_uint8(az);

    f = fullfile(data_root, run);
    N = ND2Array(f, shape, frame_offsets);

    f_new = fullfile(output_folder, [run(1:end-3) 'mat']);

    file = matfile(f_new);

    file.Properties.Writable = true;
    S = size(N);
    S(3) = S(3) - dead_frames;
    file.data = zeros(S, 'uint8');
    for t = 1:size(N, 4)
        slice = get_clean_slice(N, t);
        slice = convert_to_uint8(slice);
        if id_flip
            slice = flip(slice, 1);
            slice = flip(slice, 3);
        end
        file.data(:,:,:,t) = slice;
    end
    
    xlRange = sprintf('J%d', idx+1);
    xlswrite(filename,true,sheet,xlRange)

    
end