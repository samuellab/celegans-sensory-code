% Average traces in a single animal with 10 s prestim
% FILTER BY USE CASE is ON FOR NOW
%% Step 0: Load metadata from Excel (cell includes all steps)
clear all

CTX_load_excel;

animals = {''};

for k2 = 1:length(animals)

animal = animals{k2};
line = 'ZM10104'; %'3.2.2' '4.2.2'
% identify this animal's datasets
if strcmp(line,'ZM10104')
    row_nums = find(strcmp(sensory_status.animal,animal));
    raw_data_root = cell2mat(sensory_status.raw_root(row_nums(1)));
    output_root = cell2mat(sensory_status.analyzed_root(row_nums(1)));
    datasets = sensory_status.run(row_nums);
elseif strcmp(line,'3.2.2')
    row_nums = find(strcmp(inter1_status.animal,animal));
    raw_data_root = cell2mat(inter1_status.raw_root(row_nums(1)));
    output_root = cell2mat(inter1_status.analyzed_root(row_nums(1)));
    datasets = inter1_status.run(row_nums);
elseif strcmp(line,'4.2.2')
    row_nums = find(strcmp(inter2_status.animal,animal));
    raw_data_root = cell2mat(inter2_status.raw_root(row_nums(1)));
    output_root = cell2mat(inter2_status.analyzed_root(row_nums(1)));
    datasets = inter2_status.run(row_nums);
end

analyzed_root = fullfile(output_root,animal);

% Load proofread sets, compute delta F/F_0
frame_window = 120;
% data is a dataset x neuron structure

for i = 1:length(datasets)
    if exist(fullfile(analyzed_root,strcat(char(datasets(i)),'_prfrd_data.mat')))
    load(fullfile(analyzed_root,strcat(char(datasets(i)),'_prfrd_data.mat')))
    time = times.times;
    deltaF_Fo = zeros(size(corrected_F));
for j = 1: size(corrected_F,1) % loop for trace by trace, compute deltaF/Fo
    deltaF_Fo(j,:) = (corrected_F(j,:)-F_0(j))./F_0(j);
    data(i,j).neuron = proofread_neurons{j};
    data(i,j).use_flag = use_flag(j);
    data(i,j).full_deltaF_Fo = deltaF_Fo(j,:);
    
for k = 1:size(stimulus,1)
    start_index = find(time > (stimulus{k,2}-10), 1 );
%data(i,j).traces_int(k,:) = deltaF_Fo(j,start_index:(start_index+frame_window));
%data(i,j).times_int(k,:) = time(start_index:(start_index+frame_window));

% offset
data(i,j).traces_int(k,:) = deltaF_Fo(j,start_index-1:(start_index+frame_window-1));
data(i,j).times_int(k,:) = time(start_index-1:(start_index+frame_window-1));
end 
    
end
    end
end
% maybe we should save this?

%
% filter by use case ON currently
%data_use = data;
data_use = data([data.use_flag] == 1);
%
% determine number of neurons with data
neuron_list = unique({data_use.neuron});
for i = 1:length(neuron_list)
    data_avg(i).neuron = neuron_list{i};
    data_avg(i).fps = 2.5037; %(average time delta)
    data_temp = data_use(strcmp({data_use.neuron},neuron_list{i}));
    % allocate temporary matrix
    stim_num = size(data_temp(1).traces_int,1);
    traces_temp = zeros(length(data_temp)*stim_num,size(data_temp(1).traces_int,2));
    for j = 1:length(data_temp)
        for k = 1:stim_num
        traces_temp(k + (j-1)*stim_num,:) = data_temp(j).traces_int(k,:);
        end
    end
    data_avg(i).all_traces_int = traces_temp;
    data_avg(i).mean_trace = nanmean(traces_temp,1);
    data_avg(i).std_trace = nanstd(traces_temp,0,1);
    data_avg(i).N = size(traces_temp,1);
end

% plot all
ap_time = (0:(length(data_avg(1).mean_trace)-1))./data_avg(1).fps;

output_folder = fullfile(analyzed_root,'Traces_prestim');
make_directory(output_folder)

for n = 1:length(data_avg)

fig1 = figure(1); % plot all
for i=1:size(data_avg(n).all_traces_int,1)
   plot(ap_time,data_avg(n).all_traces_int(i,:))
   hold on
end
xlim([0 max(ap_time)])
title(strcat(data_avg(n).neuron,',',char(stimulus{1,1})))
xlabel('time (s)')
ylabel('\Delta F/F_0')
saveas(fig1,fullfile(output_folder,strcat(animal,'_',data_avg(n).neuron,'_ind_prestim.png')))
fig2 = figure(2); % plot mean
errorbar(ap_time,data_avg(n).mean_trace,data_avg(n).std_trace./sqrt(data_avg(n).N))
hold on
        % Plot the stimulus.
        stim_trace_x = [10 10 20 20];
        stim_trace_y = [-10 10 10 -10];
        fill(stim_trace_x, stim_trace_y, [0.4 0.4 0.4], ...
            'FaceAlpha', 0.6, ...
            'LineStyle', ':');
ylim([min(0,min(data_avg(n).mean_trace)) max(data_avg(n).mean_trace)*1.2])
xlim([0 max(ap_time)])
title(strcat(data_avg(n).neuron,',',char(stimulus{1,1})))
xlabel('time (s)')
ylabel('\Delta F/F_0')
saveas(fig2,fullfile(output_folder,strcat(animal,'_',data_avg(n).neuron,'_avg_prestim.png')))
close all
end
% for 1 animal, mean/std, N, near/far side per neuron

% find neurons which match string
% filter by neuron
% data(strcmp({data.neuron},'AWBR'))



%
% save data avg, stimulus structure
 save(fullfile(analyzed_root,'avg_data_prestim.mat'),'data','data_avg','stimulus');

 end

%% Plot one
n = 11;

ap_time = (0:(length(data_avg(1).mean_trace)-1))./data_avg(1).fps;

fig1 = figure(3); % plot all
for i=1:size(data_avg(n).all_traces_int,1)
   plot(ap_time,data_avg(n).all_traces_int(i,:))
   hold on
end
xlim([0 max(ap_time)])
title(strcat(data_avg(n).neuron,',',char(stimulus{1,1})))
xlabel('time (s)')
ylabel('\Delta F/F_0')

fig2 = figure(4); % plot mean
errorbar(ap_time,data_avg(n).mean_trace,data_avg(n).std_trace./sqrt(data_avg(n).N))
hold on
        % Plot the stimulus.
        stim_trace_x = [0 0 10 10];
        stim_trace_y = [-10 10 10 -10];
        fill(stim_trace_x, stim_trace_y, [0.4 0.4 0.4], ...
            'FaceAlpha', 0.6, ...
            'LineStyle', ':');
ylim([min(0,min(data_avg(n).mean_trace)) max(data_avg(n).mean_trace)*1.2])
xlim([0 max(ap_time)])
title(strcat(data_avg(n).neuron,',',char(stimulus{1,1})))
xlabel('time (s)')
ylabel('\Delta F/F_0')
