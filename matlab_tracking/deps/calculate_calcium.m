function calculate_calcium(row, data_path)

feature_size = [5, 5, 3];
points_to_keep = 25;

animal_number = row.animal;
animal_part = row.animal_part;
run_number = row.run;

animal_folder = sprintf('animal_%03d', animal_number);
dataset_id = sprintf('animal_%03d_%s', animal_number, animal_part);

if animal_part == "tail"
    run_file = sprintf('run_tail_%03d.mat', run_number);
elseif animal_part == "head"
    run_file = sprintf('run_%03d.mat', run_number);
end

root_folder = fullfile(data_path, animal_folder, dataset_id);

annotations_filename_bare = sprintf('%s_annotations.mat', ...
    run_file(1:end-4));
annotations_filename = fullfile(root_folder, annotations_filename_bare);

S = load(annotations_filename);
annotations = S.annotations;
clear S;

annotations_gcamp = annotations;

tic

A = annotations_gcamp;
A{1, 'gcamp'} = 0;

filename = fullfile(root_folder, run_file);
S = load(filename);
X = S.data;
X_S = LazyArray(X, @smooth3);
clear S;

img_size_all = size(X);
img_size = img_size_all(1:3);

parfor i = 1:size(A, 1)

    a = A(i, :);

    img = get_slice(X_S, a.t);

    coord = [a.y, a.x, a.z];

    if any(coord < 1) || any(coord > img_size)
        a.gcamp = nan;
    else

        vol = get_centered_section(img, coord, feature_size);

        values = sort(vol(:), 'descend');
        gcamp = mean(values(1:points_to_keep));

        a.gcamp = gcamp;

    end

    A(i, :) = a;

end

annotations_gcamp = A;

toc

output_filename_bare = sprintf('%s_annotations_gcamp.mat', ...
    run_file(1:end-4));
output_filename = fullfile(root_folder, output_filename_bare);

save(output_filename, 'annotations_gcamp');

%fps = 99.33/23;
%ft = 1/fps;
%frame = @(t) round(t/ft);

% stimulus = {
%     '10^{-4} 2,3-pentanedione', frame(60), frame(60+10);
%     '10^{-4} 2-butanone', frame(60+10+50), frame(60+10+50+10);
%     '200 uM NaCl', frame(60+10+50+10+50), frame(60+10+50+10+50+20)};
stimulus = get_stimulus(row, data_path);

A = annotations_gcamp;

neuron_names = unique(A.neuron_id);

size_N = length(neuron_names);
size_T = max(A.t);

gcamp = NaN(size_N, size_T);

positions = NaN(size_N, 3, size_T);

for n = 1:length(neuron_names)
    queryfn = @(x) x == neuron_names(n);
    query = arrayfun(queryfn, A.neuron_id);
    dots = sortrows(A(query, :), 't');

    gcamp(n, :) = dots.gcamp';

    for t = 1:size_T
        frame = dots(dots.t==t, :);
        positions(n, :, t) = [frame.x, frame.y, frame.z];
    end
end

filename = fullfile(root_folder, [run_file(1:end-4) '_traces.mat']);
save(filename, 'neuron_names', 'gcamp', ...
    'positions', 'fps', 'stimulus');

