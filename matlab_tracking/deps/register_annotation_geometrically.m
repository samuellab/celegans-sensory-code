function annotation = register_annotation_geometrically(...
    source_annotation, ...
    source_reference_annotations, ...
    target_reference_annotations, ...
    varargin)

TRA = target_reference_annotations;
assert(all(TRA.t == TRA.t(1)));
assert(all(TRA.c == TRA.c(1)));

default_options = struct(...
    'confidence_threshold', 0.95, ...
    'metric', [1, 1, 2], ...
    'rigid_confidence_out', 0.5, ...
    'translation_confidence_out', 0.5, ...
    'copy_confidence_out', 0.1, ...
    'force_translation', false ...
);

input_options = varargin2struct(varargin{:}); 
options = merge_struct(default_options, input_options);
g = diag(options.metric);

source_id = source_annotation.id{1};
ann_src_coord = source_annotation{:, {'x', 'y', 'z'}} * g;


good_neurons = TRA.neuron_id;
queryfn = @(x) any(x==good_neurons);
query = arrayfun(queryfn, source_reference_annotations.neuron_id);
good_src = source_reference_annotations(query, :);

src_coord = source_annotation(source_id, {'x', 'y', 'z'});
d = [];
for i = 1:size(good_src, 1)
    ref_src_coord = good_src(i, {'x', 'y', 'z'});
    d(i) = src_coord{:,:} * g * ref_src_coord{:,:}';
end
good_src.d = d';
good_src = sortrows(good_src, 'd');

if ~options.force_translation && size(good_src, 1) >= 3 % rigid
    
    good_src_3 = good_src(1:3, :);
    good_src_coords = good_src_3{:, {'x', 'y', 'z'}} * g;
    
    good_dst_coords = NaN(size(good_src_coords));
    for i = 1:size(good_src_3, 1)
        queryfn = @(x) x == good_src_3{i, 'neuron_id'};
        dst_row = TRA(arrayfun(queryfn, TRA.neuron_id), :);
        good_dst_coords(i, :) = dst_row{:, {'x', 'y', 'z'}} * g;
    end

    good_src_2d = good_src_coords(:, 1:2);
    good_dst_2d = good_dst_coords(:, 1:2);
    xy_tform = fitgeotrans(good_src_2d, good_dst_2d, ...
         'nonreflectivesimilarity');
    
    dz = mean(good_dst_coords(:, 3) - good_src_coords(:, 3));
    
    tform = eye(4);
    tform(1:2, 1:2) = xy_tform.T(1:2, 1:2);
    tform(4, 1:2) = xy_tform.T(3, 1:2);
    tform(4, 3) = dz;
    %tform = estimateRigidTransform(good_dst_coords', good_src_coords')';
    
    ann_dst_coord_h = [ann_src_coord 1] * tform;
    ann_dst_coord = ann_dst_coord_h(1:3) / g;
    
    confidence = options.rigid_confidence_out;
    
elseif size(good_src, 1) >= 1 % translation
    
    good_src_1 = good_src(1, :);
    queryfn = @(x) x == good_src_1.neuron_id;
    dst_row = TRA(arrayfun(queryfn, TRA.neuron_id), :);
    good_dst_coord = dst_row{:, {'x', 'y', 'z'}} * g;
    
    delta = good_dst_coord - good_src_1{:, {'x', 'y', 'z'}} * g;
    ann_dst_coord = (ann_src_coord + delta)/g;
    
    confidence = options.translation_confidence_out;

else % copy
    
    confidence = options.copy_confidence_out;
    
end

new_row = source_annotation;

new_row.x = ann_dst_coord(1);
new_row.y = ann_dst_coord(2);
new_row.z = ann_dst_coord(3);
new_row.t = TRA.t(1);
new_row.c = TRA.c(1);
new_row.confidence = confidence;
new_row.trace = {[]};

new_row = update_annotation_ids(new_row);

annotation = new_row;