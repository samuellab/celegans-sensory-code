classdef MemoryMappedArray < ArrayBase
    % MATFILEARRAY Immutable array that indexes into a matfile on disk.
    %
    % Author: Vivek Venkatachalam (vivekv2@gmail.com)

    properties (SetAccess = immutable)

        % Memory Mapped object on disk
        View

    end

    methods

        function obj = MemoryMappedArray(filename, dtype, shape)
            % x = MEMORYMAPPEDARRAY(filename, type, shape)
            %
            %   Creates read-only reference to data in a binary file. The
            %   number of time slices is determined by the size of the
            %   file, so only the first N-1 dimensions should be specified
            %   in shape.
            
            obj.View = memmapfile(filename, ...
                'Format', {dtype, shape, 'data'});
                
            obj.Size = [obj.View.Format{2} length(obj.View.Data)];

            obj.ElementClass = dtype;

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

            data = this.View.Data(t).data;

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
