function annotations = update_annotation_ids(annotations)

for i = 1:size(annotations, 1)
    id_tokens = {char(annotations.dataset_id(i)), ...
        char(annotations.neuron_id(i)), ...
        num2str(annotations.c(i)), ...
        num2str(annotations.t(i))};
    annotations.id{i} = strjoin(id_tokens, '_');
end

annotations.Properties.RowNames = annotations.id;