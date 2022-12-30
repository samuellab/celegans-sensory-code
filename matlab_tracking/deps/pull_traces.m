function pull_traces(folder)

feature_size = [5, 5, 3];
points_to_keep = 25;

%folder = 'D:\data\180503 23pentanedione_butanone_NaCl\animal_001\animal_001_tail';

runs = dir(fullfile(folder, 'run*.mat'));

S = load(fullfile(folder, 'annotations.mat'));
annotations = S.annotations;
clear S;

annotations_gcamp = annotations;

for idx = 1:length(annotations)
    tic

    A = annotations_gcamp{idx};
    A{1, 'gcamp'} = 0;

    filename = fullfile(folder, runs(idx).name);
    S = load(filename);
    X = S.data;
    X_S = LazyArray(X, @smooth3);
    clear S;


    parfor i = 1:size(A, 1)

        a = A(i, :);

        img = get_slice(X_S, a.t);

        vol = get_centered_section(img, [a.y, a.x, a.z], feature_size);

        values = sort(vol(:), 'descend');
        gcamp = mean(values(1:points_to_keep));

        a.gcamp = gcamp;

        A(i, :) = a;

    end

    annotations_gcamp{idx} = A;

    toc
end


output_file = fullfile(folder, 'annotations_gcamp.mat');
save(output_file, 'annotations_gcamp');

fps = 99.33/23;
ft = 1/fps;
frame = @(t) round(t/ft);

stimulus = {
    '10^{-4} 2,3-pentanedione', frame(60), frame(60+10);
    '10^{-4} 2-butanone', frame(60+10+50), frame(60+10+50+10);
    '200 uM NaCl', frame(60+10+50+10+50), frame(60+10+50+10+50+20)};

for idx = 1:length(annotations)

    A = annotations_gcamp{idx};

    neuron_names = unique(A.neuron_id);

    size_N = length(neuron_names);
    size_T = max(A.t);

    gcamp = NaN(size_N, size_T);

    positions = NaN(size_N, 3);

    for n = 1:length(neuron_names)
        queryfn = @(x) strcmp(neuron_names{n}, x);
        query = cellfun(queryfn, A.neuron_id);
        dots = sortrows(A(query, :), 't');

        gcamp(n, :) = dots.gcamp';

        frame_1 = dots(dots.t==1, :);

        positions(n, :) = [frame_1.x, frame_1.y, frame_1.z];

    end
    
    filename = fullfile(folder, ['traces_' runs(idx).name]);
    save(filename, 'neuron_names', 'gcamp', ...
        'positions', 'fps', 'stimulus');

end