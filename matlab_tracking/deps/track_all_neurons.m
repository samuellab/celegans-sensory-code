load_status;

to_track = status(status.use & status.annotations_done ...
    & ~status.neurons_tracked, :);

data_location = 'D:\data\180503 23pentanedione_butanone_NaCl';

xl_filename = 'status.xlsx';
xl_sheet = 1;

annotator_folder = fullfile(data_location, 'annotator');
annotation_file = fullfile(annotator_folder, 'annotations.json');
all_annotations_struct = load_json(annotation_file);

fns = fieldnames(all_annotations_struct);

all_annotations_struct_array = all_annotations_struct.(fns{1});
for i = 2:length(fns)
    all_annotations_struct_array(end+1) = all_annotations_struct.(fns{i});
end

all_annotations = struct2table(all_annotations_struct_array);

all_annotations.neuron_id = categorical(all_annotations.neuron_id);
all_annotations.dataset_id = categorical(all_annotations.dataset_id);

translate_coord = @(x, max, scale) x/(scale*max)*(max-1)+1;

all_annotations.x = translate_coord(all_annotations.x, 256, 4);
all_annotations.y = translate_coord(all_annotations.y, 128, 4);
all_annotations.z = translate_coord(all_annotations.z, 21, 8);

all_annotations.Properties.RowNames = all_annotations.id;

all_id_annotations = all_annotations(all_annotations.t==1, :);

perfect_confidence = 0.98;
img_reg_frame_size = [31 31 5];

for idx = 1:size(status, 1)
    row = status(idx, :);
    if ~row.annotations_done || row.neurons_tracked || ~row.use
        continue
    end

    tic
    animal_number = row.animal;
    animal_part = row.animal_part;
    run_number = row.run;
    id_run_number = row.id_run;
    id_flip = row.id_flip;

    if animal_part == "tail"
        run_file = sprintf('run_tail_%03d.mat', run_number);
    elseif animal_part == "head"
        run_file = sprintf('run%03d.mat', run_number);
    end

    if animal_part == "tail"
        id_run_file = sprintf('id_tail_%03d.nd2', id_run_number);
    elseif animal_part == "head"
        id_run_file = sprintf('id%03d.nd2', id_run_number);
    end

    dataset_id = sprintf('animal_%03d_%s', animal_number, animal_part);
    data_root = fullfile(data_location, sprintf('animal_%03d', animal_number));

    output_folder = fullfile(data_root, dataset_id);

    selected_rows = all_id_annotations.dataset_id == dataset_id;
    id_annotations = all_id_annotations(selected_rows, :);

    selected_rows = all_annotations.dataset_id == dataset_id;
    annotations = all_annotations(selected_rows, :);

    animal_status = status(status.animal==animal_number ...
        & status.animal_part==animal_part & status.use,:);
    annotator_t_idx = find(animal_status.run == run_number) + 1;

    % Get frame 1 from the id run

    A = get_array_data(BioFormatsArray(fullfile(data_root, id_run_file)));

    gcamp = A(:,:,:,4,1);

    if id_flip
       gcamp = flip(gcamp,1); 
       gcamp = flip(gcamp,3);
    end

    gcamp_small = imresize(gcamp, [128 256]);

    % View ID annotations

    run_annotations = id_annotations(:, {'x', 'y'});
    h = tight_figure(2); imshow(max_intensity_z(gcamp_small)*3); hold on; 
    plot(run_annotations.x, run_annotations.y, '.m', 'MarkerSize', 10);
    export_fig(fullfile(output_folder, 'id_annotations.png'), h);

    query = annotations.t == annotator_t_idx;

    run_annotations = annotations(query, :);
    run_annotations.t = kron(1, ones([size(run_annotations, 1), 1]));
    run_annotations.c = kron(1, ones([size(run_annotations, 1), 1]));

    run_annotations = update_annotation_ids(run_annotations);

    fiducial_annotations = run_annotations;

    A = run_annotations;
    src_annotations = A(A.t==1 & A.confidence > perfect_confidence, :);
    size_A = size(A, 1);

    image_file = fullfile(output_folder, run_file);
    X = MatfileArray(image_file, 'data');
    Z = LazyArray(X, @smooth3);
    size_T = size(X, 4);

    source_volume = get_slice(X, 1);

    parfor a_idx = 1:size_A
        source_annotation = src_annotations(a_idx, :);
        prev_annotation = source_annotation;

        for t = 2:size_T
            prev_volume = get_slice(X, t-1);
            target_volume = get_slice(X, t);

            new_row_image_p = register_annotation_image_based(...
                prev_annotation, prev_volume, target_volume, ...
                img_reg_frame_size, t);

            new_row_image = register_annotation_image_based(...
                source_annotation, source_volume, target_volume, ...
                img_reg_frame_size, t, 'guess_annotation', ...
                new_row_image_p);

            A = [A; new_row_image];

            prev_annotation = new_row_image;
        end

    end

    fiducial_annotations = A;

    A = fiducial_annotations;
    src_annotations = A(A.t==1 & A.confidence > perfect_confidence, :);

    t_ref = 1;

    query = (A.t==t_ref) & (A.confidence > perfect_confidence);
    complete_annotations = A(query, :);

    neuron_ids_todo = setdiff(id_annotations.neuron_id, ...
        complete_annotations.neuron_id);
    queryfn = @(x) any(x==neuron_ids_todo);
    query = arrayfun(queryfn, id_annotations.neuron_id);
    todo = id_annotations(query, :);

    size_A = size(todo, 1);

    A_ = A; % To help parfor variable classification.

    for a_idx = 1:size_A
        source_id_annotation = todo(a_idx, :);

        % Geom for the first time point.
        t = 1;
        target_reference_annotations = A_(A_.t==t, :);

        try
            new_row_geom = register_annotation_geometrically(...
                source_id_annotation, ...
                id_annotations, ...
                target_reference_annotations);
        catch e
            xlRange = sprintf('K%d', idx+1);
            xlswrite(xl_filename,false,xl_sheet,xlRange);
        end
        A = [A; new_row_geom];

        source_annotation = new_row_geom;
        source_reference_annotations = target_reference_annotations;
        source_volume = get_slice(Z, t);

        parfor t = 2:size_T

            target_reference_annotations = A_(A_.t==t, :);

            target_volume = get_slice(Z, t);

            new_row_geom = register_annotation_geometrically(...
                source_annotation, ...
                source_reference_annotations, ...
                target_reference_annotations, ...
                'force_translation', true);

            new_row_image = register_annotation_image_based(...
                source_annotation, source_volume, target_volume, ...
                img_reg_frame_size, t, 'guess_annotation', ...
                new_row_geom);

            A = [A; new_row_image];

        end

    end

    new_run_annotations = A;

    % Visualize
    slice = get_slice(Z, 1);

    dots = new_run_annotations; 

    im_z = autoscale(max_intensity_z(slice));

    h = tight_figure(1);
    imshow(im_z); hold on;
    plot(dots.x, dots.y, '.m', 'MarkerSize', 1);
    png_filename = sprintf('%s_annotations.png', run_file(1:end-4));
    export_fig(fullfile(output_folder, png_filename), h);

    annotations = new_run_annotations;
    annotations_filename = sprintf('%s_annotations.mat', ...
        run_file(1:end-4));
    annotation_output_file = fullfile(output_folder, annotations_filename);
    save(annotation_output_file, 'annotations');
    toc
    
    xlRange = sprintf('L%d', idx+1);
    xlswrite(xl_filename,true,xl_sheet,xlRange);
end