classdef ColorArray < ArrayBase
    %COLORARRAY Creates an array with different color channels.
    %
    % Author: Vivek Venkatachalam (vivekv2@gmail.com)

    properties (SetAccess = immutable)

        % Set of color arrays
        Colors

        PreimageSize

        % Determines the color dimension. This is typically 3 for a 2D
        % array and 4 for a 3D array.
        CDim = 3;

        % Convert images to RGB or not
        ConvertToRGB = false;

        % Number of color channels
        SizeC

        RawSize

        NSlices

        SlicedDimension

    end

    properties

        % RGB values for the colors, with each column corresponding to a
        % separate color channel.
        LUT = eye(3);

    end

    methods

        function obj = ColorArray(colors, cdim, rgb, use_rgb, sliced_dim)
            % obj = COLORARRAY(colors, cdim, rgb, use_rgb)
            %
            %   Creates an array from the color channels in colors.

            if ~iscell(colors)
                colors = {colors};
            end

            N = length(colors);

            obj.Colors = colors;
            obj.CDim = cdim;
            obj.SizeC = N;

            if nargin > 2
                obj.LUT = rgb;
            end
            obj.LUT = obj.LUT(:, 1:N);

            if nargin > 3
                obj.ConvertToRGB = use_rgb;
            end


            p = size(colors{1});

            image_size = [p(1:cdim-1) N p(cdim:end)];

            obj.PreimageSize = p;
            obj.RawSize = image_size;

            if nargin > 4
                obj.SlicedDimension = sliced_dim;
            else
                obj.SlicedDimension = length(obj.RawSize);
            end


            if obj.SlicedDimension == cdim
                obj.SlicedDimension = cdim + 1;
            end

            obj.ElementClass = element_class(colors{1});

            obj.Size = image_size;
            if obj.ConvertToRGB
                obj.Size(obj.CDim) = 3;
            end

        end

        function [varargout] = subsref(this, S)

            % Determine which slices we will need to transform
            requested = S.subs{this.SlicedDimension};

            data = zeros([this.Size(1:end-1) length(requested)], ...
                this.ElementClass);

            idx = num2cell(repmat(':', 1, length(this.Size)));
            for i = 1:length(requested)

                idx{end} = i;
                data(idx{:}) = this.get_slice(requested(i));

            end

            new_S = S;
            new_S.subs{this.SlicedDimension} = ':';
            varargout{1} = subsref(data, new_S);

        end

        function data = get_slice(this, t)

            if this.ConvertToRGB
                data = this.get_rgb_slice(t);
            else
                data = this.get_raw_slice(t);
            end

        end

        function data = get_raw_slice(this, t)

            assert(numel(t)==1, ...
                'get_slice can only be called on single slices');

            idx = num2cell(repmat(':', 1, this.SlicedDimension));
            idx{end} = 1;

            if length(this.Size) < this.SlicedDimension
                data_size = [this.RawSize 1];
                get_color = @(x) get_array_data(this.Colors{x});
            else
                data_size = [this.RawSize(1:end-1) 1];
                get_color = @(x) get_slice(this.Colors{x}, t);
            end

            data = zeros(data_size, this.ElementClass);

            slice_idx = idx;
            for c = 1:this.SizeC
                slice_idx{this.CDim} = c;
                data(slice_idx{:}) = get_color(c);
            end

        end

        function data = get_rgb_slice(this, t)

            cdata = get_raw_slice(this,t);
            data = multiply_nd(this.LUT, cdata, this.CDim);

        end

        function data = get_mipz_slice(this, t)

            data = max_intensity_z(get_rgb_slice(this, t));

        end

        function colors = get_colors(this)

            colors = this.Colors;

        end

        function data = get_array_data(this)

            if length(this.Size) < this.SlicedDimension
                data = get_slice(this, 1);
            else
                data = get_array_data(this);
            end

        end
        
        function obj = lazily_transform(this, fn)
            
            A = LazyArray(this, fn);
            
            colors = {};
            for i = 1:this.SizeC
                colors{i} = LazyArray(A, @(x) get_slice(x, i, this.CDim));
            end
            
            obj = ColorArray(colors, this.CDim, this.LUT, ...
                this.ConvertToRGB, this.SlicedDimension);
            
        end


    end

    methods (Static)


    end

end
