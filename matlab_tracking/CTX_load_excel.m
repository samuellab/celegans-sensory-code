% load data status information from Excel Sheet

if strcmp(char(java.net.InetAddress.getLocalHost.getHostName),'SAMUELLAB21')
file_location = 'D:\Dropbox\AL Data CTX\CTX Data Status.xlsx';
%elseif strcmp(char(java.net.InetAddress.getLocalHost.getHostName),'SAMUELLAB21') % Nic Tan
%file_location = 'D:\Dropbox\NT Data NG\NT Data Status.xlsx';
elseif strcmp(char(java.net.InetAddress.getLocalHost.getHostName),'samuel-lab-1') % Ray
file_location = 'I:\Dropbox\RV Data CTX\RV Data Status.xlsx';
elseif strcmp(char(java.net.InetAddress.getLocalHost.getHostName),'CORE-SAMLAB') % CORE machine
    file_location = 'E:\Dropbox\AL Data CTX\CTX Data Status M2.xlsx';
% file_location = 'E:\Dropbox\RV Data CTX\RV Data Status.xlsx'; % manually comment out which one to use
end

[~, ~, sensory_raw] = xlsread(file_location,'ZM10104 (Sensory)','A2:N1700');

sensory_raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),sensory_raw)) = {''};

%% Create table
sensory_status = {};

%% Allocate imported array to column variable names
sensory_status.raw_root = sensory_raw(:,1);
sensory_status.analyzed_root = sensory_raw(:,2);
sensory_status.animal = sensory_raw(:,3);
sensory_status.stimulus = sensory_raw(:,4);
sensory_status.idstack = sensory_raw(:,5);
sensory_status.run = sensory_raw(:,6);
sensory_status.volumes = sensory_raw(:,7);
sensory_status.IDed = sensory_raw(:,8);
sensory_status.tracked = sensory_raw(:,9);
sensory_status.traces = sensory_raw(:,10);
sensory_status.use = sensory_raw(:,13);

% to access info in a cell:
% sensory_status.field{row}
% to check if odors were delivered:
% ~isempty(strfind(sensory_status.stimulus{1}, '1-octanol')) returns 1 if
% odor was used, 0 if not

%% Clear temporary variables
clearvars sensory_raw;