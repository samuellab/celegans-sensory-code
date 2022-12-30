load_status;

filename = 'status.xlsx';
sheet = 1;
dead_frames = 2;

animal_numbers = unique(status.animal);


for animal_number = animal_numbers'
%     if animal_number~=15
%         continue
%     end

    disp(sprintf('Starting animal %d\n.', animal_number));

    animal_status = status(status.animal==animal_number & ...
        status.animal_part==part_to_use & status.use,:);
    
    if size(animal_status, 1) == 0
        continue
    end
    
    N_runs = size(animal_status, 1);
    
    dead_frames = 2;

    id_run_number = animal_status.id_run(1);
    id_flip = animal_status.id_flip(1);
    run_numbers = animal_status.run';
    frame_offsets = animal_status.offset';
    frame_counts = animal_status.frame_count';
    animal_part = part_to_use;

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

    runs = {};
    for i = 1:length(run_numbers)
        if animal_part == "tail"
            runs{i} = sprintf('run_tail_%03d.nd2', run_numbers(i));
        elseif animal_part == "head"
            runs{i} = sprintf('run%03d.nd2', run_numbers(i));
        else
            error('Unknown part: %s', animal_part);
        end     
    end

    data_root = fullfile(data_location, sprintf('animal_%03d', animal_number));

    output_folder = fullfile(data_root, dataset_id);
    make_directory(output_folder);
    
    image_folder = fullfile(output_folder, 'jpegs');
    if exist(image_folder, 'dir')
        disp(['animal ' num2str(animal_number) ' appears complete.']);
        continue
    end
    make_directory(image_folder);

    annotator_folder = fullfile(data_location, 'annotator');
    make_directory(annotator_folder);

    annotator_image_folder = fullfile(annotator_folder, 'images');
    make_directory(annotator_image_folder);

    id_shape = [256 512 21];
    
    A = get_array_data(BioFormatsArray(fullfile(data_root, id_run)));

    gcamp = A(:,:,:,4,1);

    gcamp_small = imresize(gcamp, [128 256]);
    
    
    make_original_filename = @(view, view_idx, channel_idx, t_idx) fullfile(...
        original_jpegs, ...
        sprintf('%s_%d_%d_%d.jpg', view, view_idx, channel_idx, t_idx));

    original_index_json = fullfile(original_jpegs, 'index.json');
    
    if ~exist(original_index_json, 'file')
        disp(sprintf('Animal %d has not yet been annotated by Ev.', ...
            animal_number));
        continue
    end
    
    data_idx = load_json(fullfile(original_jpegs, 'index.json'));

    % Update shape to include new volumes
    data_idx.shape_t = 1 + N_runs;
    data_idx.shape_c = 2;

    % Gather color volumes

    
    make_directory(image_folder);

    make_filename = @(view, view_idx, channel_idx, t_idx) fullfile(...
        image_folder, ...
        sprintf('%s_%d_%d_%d.jpg', view, view_idx, channel_idx, t_idx));

    process_image = @(x) uint8((x - quantile(x(:), 0.9))/2^5);
    NP_GCaMP = process_image(gcamp);
    if id_flip
        NP_GCaMP = flip(NP_GCaMP, 1);
        NP_GCaMP = flip(NP_GCaMP, 3);
    end

    write = @(x, view, view_idx, channel_idx, t_idx) imwrite(...
        x, ...
        make_filename(view, view_idx, channel_idx, t_idx), ...
        'Quality', 50);

    write_x = @(im, idx, c, t) write(im, 'X', idx, c, t);
    write_y = @(im, idx, c, t) write(im, 'Y', idx, c, t);
    write_z = @(im, idx, c, t) write(im, 'Z', idx, c, t);

    pluck_x = @(im, idx) squeeze(im(:,idx,:));
    pluck_y = @(im, idx) squeeze(im(idx,:,:))';
    pluck_z = @(im, idx) squeeze(im(:,:,idx));

    %

    t = 1; % NeuroPAL ID volume
    c = 1; % NeuroPAL / MIPs

    % Z
    for i = 1:data_idx.shape_z

        src = make_original_filename('Z',i,0,0);
        dst = make_filename('Z',i,c,t);
        copyfile(src, dst);

    end

    % Y
    for i = 1:data_idx.shape_y

        src = make_original_filename('Y',i,0,0);
        dst = make_filename('Y',i,c,t);
        copyfile(src, dst);

    end

    % X
    for i = 1:data_idx.shape_x

        src = make_original_filename('X',i,0,0);
        dst = make_filename('X',i,c,t);
        copyfile(src, dst);

    end

    %
    t = 1; % NeuroPAL ID volume
    c = 2; % GCaMP

    % Z
    for i = 1:data_idx.shape_z

        im = pluck_z(NP_GCaMP, i);
        write_z(im, i, c, t);


    end

    % Y
    for i = 1:data_idx.shape_y

        im = pluck_y(NP_GCaMP, i);
        write_y(im, i, c, t);

    end

    % X
    for i = 1:data_idx.shape_x

        im = pluck_x(NP_GCaMP, i);
        write_x(im, i, c, t);

    end

    % run slices

    % Times to take out of runs
    run_t_idxs = [1];
    t = 1;

    for run_idx = 1:length(runs)

        filename = fullfile(output_folder, [runs{run_idx}(1:end-3) 'mat']);
        S = load(filename);
        run_data = 2 * S.data(:,:,:,run_t_idxs);
        clear S;

        for j = 1:length(run_t_idxs)

            t = t + 1;

            A = run_data(:,:,:,j);

            % Z
            mipz = max_intensity_z(A);
            for i = 1:data_idx.shape_z

                im = pluck_z(A, i);
                write_z(im, i, 2, t);
                write_z(mipz, i, 1, t);

            end

            % Y
            mipy = max_intensity_y(A);
            for i = 1:data_idx.shape_y

                im = pluck_y(A, ceil(i/2));
                write_y(im, i, 2, t);
                write_y(mipy, i, 1, t);

            end

            % X
            mipx = max_intensity_x(A);
            for i = 1:data_idx.shape_x

                im = pluck_x(A, ceil(i/2));
                write_x(im, i, 2, t);
                write_x(mipx, i, 1, t);

            end

        end

    end

    save_json(data_idx, fullfile(image_folder, 'index.json'));

    % Update index


    copyfile(image_folder, fullfile(annotator_image_folder, dataset_id));

    make_dataset_index(annotator_image_folder);
    copyfile(fullfile(annotator_image_folder, 'datasets.json'), ...
        fullfile(annotator_folder, 'datasets.json'));

    % Update annotations
    all_annotations_struct = load_json(annotation_file);
    fns = fieldnames(all_annotations_struct);

    output_annotations = struct();
    for i = 1:length(fns)
        if strcmp(all_annotations_struct.(fns{i}).dataset_id, dataset_id)
            a = all_annotations_struct.(fns{i});
            a.c = 0;
            a.t = 1;
            output_annotations.(fns{i}) = a;
        end
    end

    output_filename = fullfile(annotator_folder, 'annotations.json');

    if exist(output_filename, 'file')
        current_annotations = load_json(output_filename);
    else
        current_annotations = struct();
    end

    current_annotations = merge_struct(current_annotations, output_annotations);
    save_json(current_annotations, output_filename);
    
end