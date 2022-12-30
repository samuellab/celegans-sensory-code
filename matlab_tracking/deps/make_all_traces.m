load_status;

filename = 'status.xlsx';
sheet = 1;

to_pull = status(status.neurons_tracked & ~status.traces_done, :);

local_root = 'D:\data\180503 23pentanedione_butanone_NaCl';
%remote_root = 'G:\My Drive\data\171103 neuropal microfluidics\180503 23pentanedione_butanone_NaCl';
remote_root = 'C:\Users\vivek\Dropbox\data\180503 23pentanedione_butanone_NaCl';

for idx = 1:size(status, 1)
    row = status(idx, :);
    if ~row.neurons_tracked || row.traces_done
        continue
    end
    
    calculate_calcium(row, local_root);
    copy_row_assets(row, local_root, remote_root);
    
    xlRange = sprintf('M%d', idx+1);
    xlswrite(filename,true,sheet,xlRange)
end