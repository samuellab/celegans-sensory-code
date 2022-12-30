function CTX_uint8_volume(dataset,raw_data_root,output_folder)
% Make mat file: convert raw mat files into uint8
% function which takes dataset name, animal location

load(fullfile(raw_data_root,dataset,'unprocessed_TIFs\gcamp_vol_16.mat'))

% Find a good range to make uint8

    convert_to_uint8 = {};

    quantile_target_lo = 0.9;
    quantile_target_hi = 0.9999;

    size_T = size(green_vol, 4);

    hi = 0;
    qs_lo = [];
    qs_hi = [];
    for i = 1:10

        t = randi(size_T);

        a = get_slice(green_vol,t,4);
        qs_lo(end+1) = quantile(a(:), quantile_target_lo);
        qs_hi(end+1) = quantile(a(:), quantile_target_hi);

    end

    lo = min(qs_lo);
    hi = max(qs_hi);
    convert_to_uint8 = @(x) uint8((double(x)-lo)/(hi-lo)*255);

    f_new = fullfile(output_folder, [dataset 'gcamp_vol_8.mat']);

    file = matfile(f_new);

    file.Properties.Writable = true;
    S = size(green_vol);
    file.data = zeros(S, 'uint8');
    for t = 1:size(green_vol, 4)
        slice = get_slice(green_vol, t);
        slice = convert_to_uint8(slice);
        file.data(:,:,:,t) = slice;
    end



