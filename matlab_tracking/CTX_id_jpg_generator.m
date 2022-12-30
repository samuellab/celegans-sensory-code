function CTX_id_jpg_generator(dataset,animal,raw_data_root,output_folder)
% CTX_jpg_generator

% dataset = 'run201';
% run_num = 1;
% raw_data_root = 'D:\ZM10104 2019\20190125 oct 10s';
% output_root = 'D:\Dropbox\AL Data NG\ZM10104 (Sensory)\S_003';

% HARDCODED FOR NOW
shape_x = 512;
shape_y = 256;
shape_z = 21;

image_folder = fullfile(output_folder,animal);%fullfile(raw_data_root,dataset, 'jpegs');
    make_directory(image_folder);

 % Load two mat volumes
 load(fullfile(raw_data_root,dataset,'unprocessed_TIFs\gcamp_vol_16.mat'))
 load(fullfile(raw_data_root,dataset,'unprocessed_TIFs\red_vol_16.mat'))
 
 % Gather color volumes
     make_filename = @(view, view_idx, channel_idx, t_idx) fullfile(...
        image_folder, ...
        sprintf('%s_%d_%d_%d.jpg', view, view_idx, channel_idx, t_idx));

    process_image = @(x) uint8((x - quantile(x(:), 0.9))/1);

    write = @(x, view, view_idx, channel_idx, t_idx) imwrite(...
        x, ...
        make_filename(view, view_idx, channel_idx, t_idx), ...
        'Quality', 50);

    write_x = @(im, idx, c, t) write(im, 'X', idx, c, t);
    write_y = @(im, idx, c, t) write(im, 'Y', idx, c, t);
    write_z = @(im, idx, c, t) write(im, 'Z', idx, c, t);

    pluck_x = @(im, idx) squeeze(im(:,idx,:));
    pluck_y = @(im, idx) squeeze(im(idx,:,:))';
    pluck_z = @(im, idx) squeeze(im(:,:,idx));

% run slices, we only need for 1 run

    % Times to take out of runs
    run_t_idxs = [1];
    t = 1;
% load mat file
        filename = fullfile(raw_data_root,dataset, 'unprocessed_TIFs\gcamp_vol_16.mat');%fullfile(output_folder, [runs{run_idx}(1:end-3) 'mat']);
        S = load(filename);
        green_data = process_image(S.green_vol(:,:,:,run_t_idxs));
        filename = fullfile(raw_data_root,dataset, 'unprocessed_TIFs\red_vol_16.mat');%fullfile(output_folder, [runs{run_idx}(1:end-3) 'mat']);
        R = load(filename);
        red_data = process_image(R.red_vol(:,:,:,run_t_idxs));
        clear S;
        clear R;

        % channels 1 green, 2 red, 3 green MIP
        for j = 1:length(run_t_idxs)

            %t = t + 1;

            A = green_data(:,:,:,j);
            B = red_data(:,:,:,j);

            % Z
            mipz = max_intensity_z(A);
            for i = 1:shape_z

                im = pluck_z(A, i);
                write_z(im, i, 1, t);
                write_z(mipz, i, 3, t);

                im = pluck_z(B, i);
                write_z(im, i, 2, t);
            end

            % Y
            mipy = max_intensity_y(A);
            for i = 1:shape_y

                im = pluck_y(A, ceil(i/2));
                write_y(im, i, 2, t);
                write_y(mipy, i, 1, t);

                im = pluck_y(B, i);
                write_y(im, i, 2, t);
            end

            % X
            mipx = max_intensity_x(A);
            for i = 1:shape_x

                im = pluck_x(A, ceil(i/2));
                write_x(im, i, 2, t);
                write_x(mipx, i, 1, t);

                im = pluck_x(B, i);
                write_x(im, i, 2, t);
            end

        end
        