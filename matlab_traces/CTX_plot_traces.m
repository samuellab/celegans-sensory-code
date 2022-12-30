function CTX_plot_traces(dataset,analyzed_root)

% NG_plot_traces('run102','D:\Dropbox\AL Data NG\ZM10104 (Sensory)\S_118')

%dataset = 'run103';
%analyzed_root = 'D:\Dropbox\AL Data NG\322 (Inter1)\I3_004';

output_folder = fullfile(analyzed_root,dataset);
make_directory(output_folder)

load(fullfile(analyzed_root,strcat(dataset,'_traces.mat')))
time = times.times;
% run cleaning
[gcamp_cleaned,scales,offsets] = cleanTraces(gcamp,2.3);
% using 5% currently

%
    exp_neurons = neuron_names;
    if iscategorical(exp_neurons)
        exp_neurons = cellstr(exp_neurons);
    end

% Plot traces
        % Demarcate the traces.
        for j=1:length(exp_neurons)
            if mod(j,2) == 0
                fill([0 0 max(time) max(time)], [j-1, j, j, j-1], ...
                    [0.9 0.9 0.9]);
                hold on
            end
        end
  
for i =1:length(exp_neurons)
    
plot(time,gcamp_cleaned(i,:)+i)
hold on
 mean_fbaseline = zeros(size(gcamp_cleaned(i,:))) + i;
        plot(time, mean_fbaseline, 'Color', 'k')
end

xlim([0 max(time)])
ylim([1 length(exp_neurons)+1])
yticks((1:length(exp_neurons)) + 0.5)
yticklabels(exp_neurons);
xlabel('Time (s)')
title(strcat(analyzed_root,dataset))

%stimuli y
ylimits = [0 length(exp_neurons)+1];
% Plot the stimuli.
num_stims = length(stimulus);
for i=1:num_stims
        
        % Plot the stimulus.
        stim_trace_x = [cell2mat(stimulus(i,2)), cell2mat(stimulus(i,2)), cell2mat(stimulus(i,3)), cell2mat(stimulus(i,3))];
        stim_trace_y = [ylimits(1), ylimits(2), ylimits(2), ylimits(1)];
        fill(stim_trace_x, stim_trace_y, [0.4 0.4 0.4], ...
            'FaceAlpha', 0.6, ...
            'LineStyle', ':');
        
        % Label the stimulus.
%         if i == num_stims
%             stim_mid = stim.mid(j) + stim.offset - min_trace_time;
%             text(stim_mid, ylimits(2), stim_list(iStim).name, ...
%                 'FontName', font, 'FontSize', font_size, ...
%                 'HorizontalAlignment', 'center', ...
%                 'VerticalAlignment', 'bottom');
%         end

end

savefig(fullfile(output_folder,strcat(dataset,'.fig')))
 saveas(gcf,fullfile(output_folder,strcat(dataset,'.png')))
