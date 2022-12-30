function show_annotations(annotations, A, varargin)

default_options = struct(...
    'size', [9 9 5], ...
    'index_t', [], ...
    'figure', 1 ...
);

input_options = varargin2struct(varargin{:}); 
options = merge_struct(default_options, input_options);

if isempty(options.index_t)
    options.index_t = ~(size(A, 4) == 1);
end

coords = annotations{:, {'y', 'x', 'z'}};

if options.index_t
    t = annotations.t(1);
    vol = get_slice(A, t);
else
    vol = A;
end


z = max_intensity_z(vol);
y = max_intensity_y(vol);

figure(options.figure); clf;
subplot(211);
imshow(autoscale(z)); hold on;
plot(coords(:, 2), coords(:, 1), 'mo');

subplot(212);
imshow(autoscale(y)); hold on;
plot(coords(:, 2), coords(:, 3), 'mo');

