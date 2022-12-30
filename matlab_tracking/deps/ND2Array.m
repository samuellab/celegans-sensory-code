classdef ND2Array < ArrayBase
    % ND2ARRAY Immutable array that indexes into an ND2 file.
    %
    % Author: Vivek Venkatachalam (vivekv2@gmail.com)

    properties (SetAccess = immutable)

        % Memory Mapped object on disk
        View

    end

    methods

        function obj = ND2Array(filename, shape, phase)
            % x = ND2ARRAY(filename, type, shape)
            %
            %   Creates read-only reference to data in an ND2 file.
            
            if nargin == 2
                phase = 0;
            end
            
            fid = fopen(filename, 'r');
            fseek(fid, 43500*2, 'bof');
            while ~fread(fid, 1, 'uint16')
            end
            fread(fid, 2, 'uint16');
            fread(fid, (shape(2)*shape(1)+2048)*phase, 'uint16');
            offset = ftell(fid);
            fclose(fid);
            
            obj.View = memmapfile(filename, ...
                'Offset', offset, ...
                'Format', ...
                    {'uint16', [shape(2) shape(1)], 'data'; ...
                     'uint16', 2048, 'junk'});
                
            obj.Size = shape;

            obj.ElementClass = 'uint16';

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

        function data = get_slice(this, t)

            assert(numel(t)==1, ...
                'get_slice can only be called on single slices');
            
            start_idx = 1 + (t-1)*this.Size(3);
            end_idx = start_idx + this.Size(3)-1;
            idx = start_idx:end_idx;
            data = zeros(this.Size(1:3), this.ElementClass);
            for i = 1:length(idx)
                data(:,:,i) = this.View.Data(idx(i)).data';
            end

        end

    end
    
    methods(Static)
        
        function test(x)
            disp(x);
        end
    
        function obj = using_matfile(location, datatype)
            
            filename = [location '.bin'];
            matfile_name = [location '.mat'];
            S = load(matfile_name);
            
            if isfield(S, 'C')
                shape = [S.W S.H S.C S.D];
            else
                shape = [S.W S.H S.D];
            end
            
            if nargin == 2
                dtype = datatype;
            elseif isfield(S, 'dtype')
                dtype = S.dtype;
            else
                dtype = 'uint16';
            end
            
            obj = MemoryMappedArray(filename, dtype, shape);
            
        end
        
    end

end
