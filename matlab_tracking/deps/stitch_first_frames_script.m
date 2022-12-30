animal_number = 18;
animal_part = 'tail';


%
load_tail_status;

animal_status = status(status.animal == animal_number & status.use, :);


id_run_number = animal_status.id_run(1);
id_flip = animal_status.id_flip(1);
run_numbers = animal_status.run';
frame_offsets = animal_status.offset';
frame_counts = animal_status.offset_frames';

%
if strcmp(animal_part, 'tail')
    id_run = sprintf('id_tail_%03d.nd2', id_run_number);
elseif strcmp(animal_part, 'head')
    id_run = sprintf('id%03d.nd2', id_run_number);
else
    error('Unknown part: %s', animal_part);
end   

runs = {};
for i = 1:length(run_numbers)
    if strcmp(animal_part, 'tail')
        runs{i} = sprintf('run_tail_%03d.mat', run_numbers(i));
    elseif strcmp(animal_part, 'head')
        runs{i} = sprintf('run%03d.mat', run_numbers(i));
    else
        error('Unknown part: %s', animal_part);
    end     
end

dataset_id = sprintf('animal_%03d_%s', animal_number, animal_part);

data_location_w = 'D:\data\180503 23pentanedione_butanone_NaCl';
data_location_g = 'G:\My Drive\data\171103 neuropal microfluidics\180503 23pentanedione_butanone_NaCl';

data_location = data_location_w;

data_root = fullfile(data_location, sprintf('animal_%03d', animal_number));

output_folder = fullfile(data_root, dataset_id);
make_directory(output_folder);

annotator_folder = fullfile(data_location, 'annotator');
annotation_file = fullfile(annotator_folder, 'annotations.json');

annotator_image_folder = fullfile(annotator_folder, 'images');
make_directory(annotator_image_folder);



% Place annotations in a table

all_annotations_struct = load_json(annotation_file);

fns = fieldnames(all_annotations_struct);

all_annotations = all_annotations_struct.(fns{1});
for i = 2:length(fns)
    all_annotations(end+1) = all_annotations_struct.(fns{i});
end

annotations = struct2table(all_annotations(...
    strcmp({all_annotations.dataset_id}, dataset_id)));

translate_coord = @(x, max, scale) x/(scale*max)*(max-1)+1;

annotations.x = translate_coord(annotations.x, 256, 4);
annotations.y = translate_coord(annotations.y, 128, 4);
annotations.z = translate_coord(annotations.z, 21, 8);

annotations.Properties.RowNames = annotations.id;

id_annotations = annotations(annotations.t==1, :);

% Get frame 1 from the id run

A = get_array_data(BioFormatsArray(fullfile(data_root, id_run)));

gcamp = A(:,:,:,4,1);

if id_flip
   gcamp = flip(gcamp,1); 
   gcamp = flip(gcamp,3);
end

gcamp_small = imresize(gcamp, [128 256]);



% View ID annotations

frame_1 = id_annotations(:, {'x', 'y'});
tight_figure(2); imshow(max_intensity_z(gcamp_small)*3); hold on; 
plot(frame_1.x, frame_1.y, '.m', 'MarkerSize', 10);

% Set up run annotations for each individual run

run_annotations = {};

for i = 1:length(runs)

    run_annotations{i} = annotations({},:);
    query = annotations.t == 1 + i;

    frame_1 = annotations(query, :);
    frame_1.t = kron(1, ones([size(frame_1, 1), 1]));
    frame_1.c = kron(i, ones([size(frame_1, 1), 1]));
    run_annotations{i} = [run_annotations{i}; frame_1];

    run_annotations{i} = update_annotation_ids(run_annotations{i});
end

%% Track reference neurons through time series using img

visualize = false;

perfect_confidence = 0.98;

fiducial_annotations = run_annotations;

img_reg_frame_size = [31 31 5];

tic

for run_idx = 1:length(runs)

    A = run_annotations{run_idx};
    src_annotations = A(A.t==1 & A.confidence > perfect_confidence, :);
    size_A = size(A, 1);

    image_file = fullfile(output_folder, runs{run_idx});
    X = MatfileArray(image_file, 'data');
    source_volume = get_slice(X, 1);
    size_T = size(X, 4);
    
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

            if visualize

                show_annotation(source_annotation, source_volume, 'figure', 1);
                show_annotation(source_annotation, source_volume, ...
                    'size', [21 21 9], 'figure', 2);
                show_annotations(source_annotation, source_volume, ...
                    'figure', 3);
                
                show_annotation(new_row_image_p, target_volume, 'figure', 4);
                show_annotation(new_row_image_p, target_volume, ...
                    'size', [21 21 9], 'figure', 5);
                show_annotations(new_row_image_p, target_volume, ...
                    'figure', 6);

                show_annotation(new_row_image, target_volume, 'figure', 7);
                show_annotation(new_row_image, target_volume, ...
                    'size', [21 21 9], 'figure', 8);
                show_annotations(new_row_image, target_volume, ...
                    'figure', 9);
                drawnow;
%                pause;
%                 show_annotation(new_row_bright, target_volume, 'figure', 1);
%                 show_annotation(new_row_bright, target_volume, 'size', [21 21 9], ...
%                     'figure', 2);
%                 pause;

            end
            
            prev_annotation = new_row_image;
        end

    end

    fiducial_annotations{run_idx} = A;

end

toc

%% Visualize

run_idx = 1;

image_file = fullfile(output_folder, runs{run_idx});
A = MatfileArray(image_file, 'data');
slice = get_slice(A, 1);

dots = fiducial_annotations{run_idx}; 

im_z = max_intensity_z(slice);
im_y = max_intensity_y(slice);

figure(1); clf;
subplot(211);
imshow(im_z*5); hold on;
plot(dots.x, dots.y, '.m', 'MarkerSize', 10);

subplot(212);
imshow(im_y*5); hold on;
plot(dots.x, dots.z, '.m', 'MarkerSize', 10);

%% Visualize in time

run_idx = 1;

image_file = fullfile(output_folder, runs{run_idx});
fid_dots = fiducial_annotations{run_idx}; 
A = MatfileArray(image_file, 'data');
slice = get_slice(A, 1);

for t = 1:5:size(A, 4)
    slice = get_slice(A, t);

    dots = fid_dots(fid_dots.t==t,:);

    im_z = max_intensity_z(slice);
    im_y = max_intensity_y(slice);

    figure(2); clf;
    subplot(211);
    imshow(im_z*5); hold on;
    plot(dots.x, dots.y, '.r', 'MarkerSize', 10);

    subplot(212);
    imshow(im_y*5); hold on;
    plot(dots.x, dots.z, '.r', 'MarkerSize', 10);
    pause(0.001);
    
end

%% Track the rest of the neurons through time series using geom

perfect_confidence = 0.98;

new_run_annotations = fiducial_annotations;

tic


for run_idx = 1:length(runs)

    A = fiducial_annotations{run_idx};
    src_annotations = A(A.t==1 & A.confidence > perfect_confidence, :);
    size_A = size(A, 1);

    image_file = fullfile(output_folder, runs{run_idx});
    X = MatfileArray(image_file, 'data');
    X_S = CachedArray(LazyArray(X, @smooth3));
    Z = LazyArray(X, @smooth3);
    %Z = get_array_data(X_S);
    size_T = size(X, 4);

    t_ref = 1;

    query = (A.t==t_ref) & (A.confidence > perfect_confidence);
    complete_annotations = A(query, :);

    neuron_ids_todo = setdiff(id_annotations.neuron_id, ...
        complete_annotations.neuron_id);
    queryfn = @(x) any(strcmp(x, neuron_ids_todo));
    query = cellfun(queryfn, id_annotations.neuron_id);
    todo = id_annotations(query, :);
    
    size_A = size(todo, 1);
    
    A_ = A; % To help parfor variable classification.
    
    for a_idx = 1:size_A
        source_id_annotation = todo(a_idx, :);
        
        % Geom for the first time point.
        t = 1;
        target_reference_annotations = A_(A_.t==t, :);
        
        new_row_geom = register_annotation_geometrically(...
            source_id_annotation, ...
            id_annotations, ...
            target_reference_annotations);
        
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

    new_run_annotations{run_idx} = A;

end

toc


%% Visualize

run_idx = 1;

image_file = fullfile(output_folder, runs{run_idx});
A = MatfileArray(image_file, 'data');
slice = get_slice(A, 1);

dots = new_run_annotations{run_idx}; 

im_z = max_intensity_z(slice);
im_y = max_intensity_y(slice);

figure(1); clf;
subplot(211);
imshow(im_z*5); hold on;
plot(dots.x, dots.y, '.m', 'MarkerSize', 10);

subplot(212);
imshow(im_y*5); hold on;
plot(dots.x, dots.z, '.m', 'MarkerSize', 10);

%% Visualize in time

run_idx = 1;

image_file = fullfile(output_folder, runs{run_idx});
all_dots = new_run_annotations{run_idx}; 
A = MatfileArray(image_file, 'data');

for t = 1:10:size(A, 4)
    slice = get_slice(A, t);

    dots = all_dots(all_dots.t==t,:);

    im_z = max_intensity_z(slice);
    im_y = max_intensity_y(slice);

    figure(1); clf;
    subplot(211);
    imshow(im_z*5); hold on;
    plot(dots.x, dots.y, '.m', 'MarkerSize', 10);

    subplot(212);
    imshow(im_y*5); hold on;
    plot(dots.x, dots.z, '.m', 'MarkerSize', 10);
    pause(0.01);
    
end

%% Update spots using image registration

annotations = new_run_annotations;
annotation_output_file = fullfile(output_folder, 'annotations.mat');
save(annotation_output_file, 'annotations');
