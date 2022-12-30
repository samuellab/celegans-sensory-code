function [traces, scales, offsets] = cleanTraces(traces, fps, varargin)
%CLEANTRACES Clean up the neural traces by removing outliers & bleaching,
%            then scaling the traces to [0,1].
%
%   [TRACES, SCALES, OFFSETS] = CLEANTRACES(TRACES, FPS)
%   [TRACES, SCALES, OFFSETS] = CLEANTRACES(TRACES, FPS, SIGMA)
%   [TRACES, SCALES, OFFSETS] = CLEANTRACES(TRACES, FPS, SIGMA, DEBLEACH)
%   [TRACES, SCALES, OFFSETS] = CLEANTRACES(TRACES, FPS, SIGMA, DEBLEACH,
%                                           INTERP)
%   [TRACES, SCALES, OFFSETS] = CLEANTRACES(TRACES, SIGMA, DEBLEACH, FPS,
%                                           INTERP, SMOOTH, WINDOW)
%
%   Inputs:
%   traces   - the neural traces to clean
%   fps      - frames/second
%   sigma    - remove outliers wherein a single frame is sigma standard
%              deviations from the mean ([] = don't remove outliers).
%              default: 10
%   debleach - detrend the traces to remove bleaching
%              default: 2
%              0 = leave the data as is (don't detrend)
%              1 = detrend the data by computing a global bleach curve
%              2 = detrend the data by computing individual bleach curves
%              3 = detrend the data by computing individual bleach curves
%                  + compute dF/F0 using F0 = 5th percentile.
%   interp   - interpolate missing data with the specified method
%              (see interp1)
%              default: [] (no interpolation)
%   smooth   - smooth the traces with the specified method
%              (see smoothdata - you MUST specify the window size as well)
%              causal = causal bandpass
%              high = high pass
%              low = low pass
%              default: [] (no smoothing)
%   window   - the window size to use for filtering
%              default: [] (no smoothing)
%
%   Outputs:
%   traces   - the scaled neural traces, cleaned up
%   scales   - the scales used to normalize to unit traces
%   offsets  - the offsets used to normalize to unit traces

% Are we removing extreme outliers?
sigma_threshold = 10;
if ~isempty(varargin) && ~isempty(varargin{1})
    sigma_threshold = varargin{1};
end

% Are we detrending the traces to remove bleaching?
detrend_mode = 2;
if length(varargin) > 1
    detrend_mode = varargin{2};
end

% Are we interpolating missing data?
interp_method = [];
if length(varargin) > 2
    interp_method = varargin{3};
end

% Are we smoothing the data?
smooth_method = [];
smooth_window = [];
if length(varargin) > 4
    smooth_method = varargin{4};
    smooth_window = varargin{5};
end

% Clean up the data by removing the first & last frames (outliers).
traces(:,1) = nan;
traces(:,end) = nan;

% Dampen extreme outliers: a single frame wherein the signal change is
% larger than sigma_threshold standard deviations from the mean.
if ~isempty(sigma_threshold) && abs(sigma_threshold) > 0
    
    % Compute the thresholds for ouliers.
    extreme_threshold = sigma_threshold * nanstd(traces, 0, 2) + ...
        nanmean(traces, 2);
    diff_traces = diff(traces, 1, 2);
    
    % Find the outliers.
    extreme_max_traces = diff_traces > extreme_threshold;
    extreme_min_traces = diff_traces < -extreme_threshold;
    [iExtreme_neurons, iExtreme_frames] = find(...
        (extreme_max_traces(:,1:(end-1)) & extreme_min_traces(:,2:end)) ...
        | (extreme_min_traces(:,1:(end-1)) & extreme_max_traces(:,2:end)));
    iExtreme_frames = iExtreme_frames + 1;
    
    % Remove the outliers.
    for i = 1:length(iExtreme_neurons)
        traces(iExtreme_neurons(i),iExtreme_frames(i)) = nan;
    end
    
    % Median filter nearest neighbors..
    %traces = medfilt1(traces, 3, [], 2);
end

% Initialize the trace offsets & scaling.
offsets = zeros(size(traces,1),1);
detrend_offsets = zeros(size(traces,1),1);
scales = ones(size(traces,1),1);

% Detrend the traces to remove bleaching.
x = 1:size(traces,2);
if detrend_mode > 0

    % Determine F0 for the traces.
    % Note: this must go here, after de-bleaching you may invert traces
    % if F0 < 0;
    F0 = prctile(traces, 5, 2);
    
    % Determine the filter order for median filtering (10s).
    filt_order = round(10 * fps);
    
    % Do we have enough data to compute the bleach curve?
    detrend_threshold = 0.1 * size(traces,2); % 10% of the data

    % Assume a uniform bleach curve & detrend.
    if detrend_mode == 1
        
        % Scale the traces.
        offsets = nanmin(traces, [], 2);
        traces = traces - offsets;
        scales = nanmax(traces, [], 2);
        traces = traces ./ scales;
        
        % Compute the global bleach curve.
        y = nanmean(traces, 1);
        y_filt = medfilt1(y, filt_order, 'omitnan');
        y_filt_data = ~isnan(y_filt);
        y_data = ~isnan(y);
        if sum(y_data) > detrend_threshold
            %[p, ~, mu] = polyfit(x(y_data), y(y_data), 2);
            %f_y = polyval(p, x, [] ,mu);
            f = fit(x(y_filt_data)', y_filt(y_filt_data)', 'exp1');

            % Detrend the traces to remove bleaching.
            if f.b < 0 % bleach curves must decay
                f_y = f.a * exp(f.b * x);
                detrend_offsets(:) = f.a;
                traces = traces - f_y;
            end
        end
        
    % Compute individual bleach curves, per neuron, & detrend.
    else
        
        % Compute individual bleach curves, per neuron, & detrend.
        for i = 1:size(traces,1)
            yi = traces(i,:);
            yi_filt = medfilt1(yi, filt_order, 'omitnan');
            yi_filt_data = ~isnan(yi_filt);
            
            % Do we have enough data to use a local bleach curve?
            if sum(yi_filt_data) > detrend_threshold
                fi = fit(x(yi_filt_data)', yi_filt(yi_filt_data)', 'exp1');
                
                % Detrend the traces to remove bleaching.
                if fi.b < 0 % bleach curves must decay
                    fi_yi = fi.a * exp(fi.b * x);
                    detrend_offsets(i) = fi.a;
                    traces(i,:) = (yi - fi_yi);
                end
            end
        end
    end
    
    % Compute dF/F0 using F0 = 5th percentile.
    if detrend_mode == 3
        F0 = zeros(size(traces,1),1) + nanmedian(F0); % use the median F0
        scales = F0;
        %scales(scales < 1) = 1; % dNMF ~ [0,1]
        offsets = zeros(size(F0));
        traces = (traces - offsets) ./ scales;
    end
end

% Interpolate missing data.
if ~isempty(interp_method)
    for i = 1:size(traces,1)
        
        % Find the missing data.
        nan_data = isnan(traces(i,:));
        
        % Is there any data to work with?
        if sum(nan_data) < size(traces,2)
            
            % Interpolate the missing data.
            y_data = ~nan_data;
            traces(i,nan_data) = interp1(x(y_data), traces(i,y_data), ...
                x(nan_data), interp_method);
        end
    end
end

% Smooth the data.
if ~isempty(smooth_method) && ~isempty(smooth_window)
    
    % Causal bandpass filter.
    if strcmpi(smooth_method, 'causal')
        traces = causalBandpassFilter(traces, ...
            smooth_window(1), smooth_window(2), smooth_window(3));

    % Highpass filter.
    elseif strcmpi(smooth_method, 'high')
        traces = highpassFilter(traces, smooth_window(1), smooth_window(2));
        
    % Lowpass filter.
    elseif strcmpi(smooth_method, 'low')
        traces = lowpassFilter(traces, smooth_window(1), smooth_window(2));
        
    % Standard filter.
    else
        traces = smoothdata(traces, 2, smooth_method, smooth_window, ...
            'includenan');
    end
end

% Re-scale the traces.
if detrend_mode < 3
    new_offsets = nanmin(traces, [], 2);
    traces = traces - new_offsets;
    new_scales = nanmax(traces, [], 2);
    traces = traces ./ new_scales;
    
    % Compute the compounded scales & offsets.
    offsets = offsets + (detrend_offsets + new_offsets) .* scales;
    scales = scales .* new_scales;
    
    % Rescale the traces to [0.05, 0.95].
    %traces = traces * 0.9 + 0.05;
end
end
