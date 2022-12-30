function show_annotation(annotation, A, varargin)

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

coords = annotation{:, {'y', 'x', 'z'}};
coords = round(coords);

if options.index_t
    t = annotation.t;
    slice = get_slice(A, t);
else
    slice = A;
end

offset = floor(options.size/2);

vol = get_image_section(coords-offset, options.size-1, slice);

z = max_intensity_z(vol);
y = max_intensity_y(vol);

figure(options.figure); clf;
subplot(211);
imshow(autoscale(z)); hold on;
plot(offset(2)+1, offset(1)+1, 'mo');

subplot(212);
imshow(autoscale(y)); hold on;
plot(offset(2)+1, offset(3)+1, 'mo');

