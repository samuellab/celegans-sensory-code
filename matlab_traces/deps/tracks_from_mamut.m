function tracks = tracks_from_mamut(mamut_file, pixel_size)
% tracks = TRACKS_FROM_MAMUT(mamut_file)
%
%   Convert tracks from a mamut file into a matlab table.
%
% Author: Vivek Venkatachalam (vivekv2@gmail.com)

if nargin < 2
    pixel_size = [1 1 2];
end

if exist(mamut_file, 'dir')
    [p,n,e] = fileparts(mamut_file);
    mamut_file = fullfile(p, [n e], [n '-mamut.xml']);
end

hash = Simulink.getFileChecksum(mamut_file);
[p, n, ~] = fileparts(mamut_file);
cached_file = fullfile(p, ['_' n '-' hash '.mat']);

if exist(cached_file, 'file')
    S = load(cached_file);
    tracks = S.tracks;
else

    tracking_data = readMamutXML(mamut_file);

    tracks = array2table(tracking_data,...
        'VariableNames', ...
        {'id', 'type', 'y', 'x', 'z', 'radius', 'parent_id', 'time', ...
         'confidence', 'skeletonId'});

    tracks.time = tracks.time + 1;
    tracks.x = tracks.x/pixel_size(2) + 1;
    tracks.y = tracks.y/pixel_size(1) + 1;
    tracks.z = tracks.z/pixel_size(3) + 1;
    
    save(cached_file, 'tracks');
end