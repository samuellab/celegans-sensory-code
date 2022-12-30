% CTX dataset proofreader
% proofreads IDs for a single dataset

%% load existing traces file
%dataset = 'run601';
%analyzed_root = 'D:\Dropbox\AL Data NG\322 (Inter1)\I3_002';

clear all

dataset = 'run101';
%analyzed_root = 'E:\Dropbox\AL Data CTX\ZM10104 (Sensory)\A5_006';
analyzed_root = 'E:\Dropbox\HC Data CTX\ZM10104 (Sensory)\HC_D_001';
%analyzed_root = 'E:\Dropbox\RV Data CTX\ZM10104 (Sensory)\RV_F_011';

output_folder = fullfile(analyzed_root,dataset);
make_directory(output_folder)

% clean traces, generate deltaF/F_o
load(fullfile(analyzed_root,strcat(dataset,'_traces.mat')))
time = times.times;
% run cleaning
[gcamp_cleaned,scales,offsets] = cleanTraces(gcamp,2.3,10,2);
% currently cleaning scales from 0 to 1, should we be using deltaF/F
% instead? Can do this using the scales and offsets:
% to get to deltaF/F_0: F_0 = offset, deltaF = gcamp_cleaned.*scale 
corrected_F_new = zeros(size(gcamp));
F_0 = zeros(size(gcamp,1));
for j = 1: size(gcamp,1)
    % bleach corrected raw fluorescence
    % compute deltaF_F0 when averaging
    corrected_F_new(j,:) = gcamp_cleaned(j,:).*scales(j)+offsets(j);
    % First guess for F_0 is median
    F_0(j) = scales(j,:);
end
%
    tracked_IDs = neuron_names;
    if iscategorical(tracked_IDs)
        tracked_IDs = cellstr(tracked_IDs);
    end

% clone important variables, save new structure, generate flag
proofread_file = fullfile(analyzed_root,strcat(dataset,'_prfrd_data.mat'));

if exist(proofread_file)>=1
load(proofread_file)
else
proofread_neurons = tracked_IDs;
use_flag = zeros(size(tracked_IDs)); %(0 auto, 1 good, -1 don't use)
end
time = times.times;

corrected_F = corrected_F_new;


% Save raw traces, cleaned traces, times, proofread neuron names, raw tracked names, stimulus
save(proofread_file, 'tracked_IDs','proofread_neurons','use_flag', 'corrected_F','F_0', ...
    'positions', 'times', 'stimulus');
% plot all traces
% Plot traces
        % Demarcate the traces.
        figure(3)
%         for j=1:length(proofread_neurons)
%             if mod(j,2) == 0
%                 fill([0 0 max(time) max(time)], [j-1, j, j, j-1], ...
%                     [0.9 0.9 0.9]);
%                 hold on
%             end
%         end
  
for i =1:length(proofread_neurons)
    
plot(time,gcamp_cleaned(i,:)/(max(gcamp_cleaned(i,:)))+i,'LineWidth',2)
hold on
  mean_fbaseline = zeros(size(gcamp_cleaned(i,:))) + i;
         plot(time, mean_fbaseline, 'Color', 'k')
end

xlim([0 max(time)])
ylim([1 length(proofread_neurons)+1])
yticks((1:length(proofread_neurons)) + 0.5)
yticklabels(proofread_neurons);
xlabel('Time (s)')
title(dataset)

set(gca,'FontSize',16)
%stimuli y
ylimits = [0 length(proofread_neurons)+1];
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

% add an auto-save figure feature

%% Proofreader GUI
% use GUI in concert with annotator window and full traces
% display deltaF/F0 trace


r_max = length(tracked_IDs);

h = figure();


r = 1;
while r <= r_max
   plot(time,corrected_F(r,:)) 
   hold on
   if ~isnan(F_0(r))
   yline(F_0(r),'-.b')
   ylim([0 1.2*max(corrected_F(r,:))])
   end
   hold on
   xlim([0 max(time)])

%yticklabels(exp_neurons);
xlabel('Time (s)')
title(strcat('Tracked ID:',tracked_IDs{r},'Proofread ID:',proofread_neurons{r},', F_0 =',num2str(F_0(r))))
% change plot size
if use_flag(r) == 0
set(gcf,'color','y')
elseif use_flag(r) == 1
set(gcf,'color','g')
elseif use_flag(r) == -1
set(gcf,'color','r')
end

%stimuli y
ylimits = [0 1.2*max(corrected_F(r,:))];
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
hold off

was_a_key = waitforbuttonpress;
    if was_a_key && strcmp(get(h, 'CurrentKey'), 'leftarrow')
      r = max(1,r - 1);
    elseif was_a_key && strcmp(get(h, 'CurrentKey'), 'rightarrow')
      r = min(r + 1,r_max);
    elseif was_a_key && strcmp(get(h, 'CurrentKey'), 'e') % edit neurons
prompt = {'What is the correct neuron?'};
dlgtitle = 'Input';
dims = [1 35];
definput = {proofread_neurons{r}};
answer = char(inputdlg(prompt,dlgtitle,dims,definput));
% clear other places where this neuron is assigned
if sum(contains(proofread_neurons,answer))>=1
proofread_neurons{contains(proofread_neurons,answer)}='reassigned';
end
proofread_neurons{r} = answer;
    elseif was_a_key && strcmp(get(h, 'CurrentKey'), 'g') % flag positive
        use_flag(r) = 1;
    elseif was_a_key && strcmp(get(h, 'CurrentKey'), 'f') % flag neutral
        use_flag(r) = 0;
    elseif was_a_key && strcmp(get(h, 'CurrentKey'), 'd') % flag negative
        use_flag(r) = -1;
    elseif was_a_key && strcmp(get(h, 'CurrentKey'), 'b') % adjust baseline F_0
prompt = {'How much should the baseline F_0 be adjusted?'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'0'};
delta = str2double(inputdlg(prompt,dlgtitle,dims,definput));
F_0(r) = F_0(r) + delta;
    elseif was_a_key && strcmp(get(h, 'CurrentKey'), 's') % save data
        save(proofread_file, 'tracked_IDs','proofread_neurons','use_flag', 'corrected_F','F_0', ...
    'positions', 'times', 'stimulus');
    else % add save button
        return
    end

end
