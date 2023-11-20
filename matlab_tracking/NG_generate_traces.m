function NG_generate_traces(dataset,output_root)
% NG_generate_traces

%dataset = 'run201';
%output_root = 'D:\Dropbox\AL Data NG\ZM10104 (Sensory)\S_003';


feature_size = [5, 5, 3];
points_to_keep = 10;

run_file = strcat(dataset,'gcamp_vol_8.mat');

annotations_filename = fullfile(output_root, strcat(dataset,'_annotations.mat'));

S = load(annotations_filename);
annotations = S.annotations;
clear S;

annotations_gcamp = annotations;

tic


A = annotations_gcamp;
A{1, 'gcamp'} = 0;

filename = fullfile(output_root, run_file);
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

save(fullfile(output_root,strcat(dataset,'_annotations_gcamp.mat')), 'annotations_gcamp');

load(fullfile(output_root,strcat(dataset,'metadata.mat')));

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

filename = fullfile(output_root,strcat(dataset,'_traces.mat'));
save(filename, 'neuron_names', 'gcamp', ...
    'positions', 'times', 'stimulus');
