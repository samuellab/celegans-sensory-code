classdef BigDataViewerArray < ArrayBase
    % BIGDATAVIEWERARRAY Immutable array that indexes into a BigDataViewer
    % file on disk
    %
    % Authors: Vivek Venkatachalam (vivekv2@gmail.com) 
    %          Ben Sorscher        (bsorsch@gmail.com)
    
    properties
        
        % Directory containing the BigDataViewer HDF5 file along with the
        % XML metadata file.
        Path
        
    end
    
    methods
        
        function obj = BigDataViewerArray(path)
            % x = BIGDATAVIEWERARRAY(path)
            %
            %   Creates a read-only reference to the big data viewer file
            %   located at the specified location. There must be both an
            %   HDF5 file and an associated XML file that conform to the
            %   BigDataViewer specification.
            
            obj.Path = path;
            
            S = xml2struct(obj.get_xml_file());
            
            size_C = length(...
                S.SpimData.SequenceDescription.ViewSetups.ViewSetup);
            
            size_raw = strsplit(...
                S.SpimData.SequenceDescription.ViewSetups.ViewSetup{1} ...
                .size.Text, ' ');
            
            size_Y = str2num(size_raw{1});
            size_X = str2num(size_raw{2});
            size_Z = str2num(size_raw{3});
            
            time_span = S.SpimData.SequenceDescription.Timepoints;
            
            start_time = str2num(time_span.first.Text);
            end_time = str2num(time_span.last.Text);
            
            size_T = end_time - start_time + 1;
            
            obj.Size = [size_Y size_X size_Z size_C size_T];
            
            obj.ElementClass = 'int16';
            
        end
        
        function y = get_h5_file(this)
            [p,d,~] = fileparts(this.Path);
            y = fullfile(p,d,[d '.h5']);
        end
        
        function y = get_xml_file(this)
            [p,d,~] = fileparts(this.Path);
            y = fullfile(p,d,[d '.xml']);
        end
        
        function y = get_mamut_file(this)
            [p,d,~] = fileparts(this.Path);
            y = fullfile(p,d,[d '-mamut.xml']);
        end
        
        function [varargout] = subsref(this, S)

            % Determine which slices we will need to transform
            requested = S.subs{ndims(this)};

            % Expand ':' for sliced dimension.
            if ischar(requested) && requested == ':'
                requested = 1:this.Size(end);
            end

            data = zeros([this.Size(1:end-1) length(requested)], ...
                this.ElementClass);

            idx = num2cell(repmat(':', 1, length(this.Size)));
            for i = 1:length(requested)

                idx{end} = i;
                data(idx{:}) = this.get_slice(requested(i));

            end

            new_S = S;
            new_S.subs{ndims(this)} = ':';
            varargout{1} = subsref(data, new_S);

        end
        
        function data = get_channel_slice(this, c, t)
            data = h5read(this.get_h5_file(), ...
                sprintf('/t%05d/s%02d/0/cells', t-1, c-1));
        end

        function data = get_slice(this, t)

            assert(numel(t)==1, ...
                'get_slice can only be called on single slices');

            S = this.Size;
            data = zeros(S(1), S(2), S(3), S(4), this.ElementClass);
            
            for c = 1:S(4)
                data(:,:,:,c) = this.get_channel_slice(c,t);
            end

        end
        
    end
    
end

