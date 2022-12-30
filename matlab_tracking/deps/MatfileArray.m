classdef MatfileArray < ArrayBase
    % MATFILEARRAY Immutable array that indexes into a matfile on disk.
    %
    % Author: Vivek Venkatachalam (vivekv2@gmail.com)

    properties (SetAccess = immutable)

        % matfile object that corresponds to where the data is located on
        % disk
        Matfile

        % Field in matfile this array will access.
        Field

    end

    methods

        function obj = MatfileArray(filename, field)
            % x = MATFILEARRAY(filename, field)
            %
            %   Creates read-only reference to the array contained in a
            %   specified field of a .mat file.

            if isa(filename, 'matlab.io.MatFile')
                obj.Matfile = filename;
            elseif isa(filename, 'char')
                obj.Matfile = matfile(filename);
            else
                error('Invalid argument sent to Matfile constructor.');
            end
            
            if nargin < 2
                field = 'data';
            end

            obj.Field = field;
            obj.Size = size(obj.Matfile, field);

            % Determine the type and size of the array.
            mfile_info = whos(obj.Matfile);

            for i = 1:length(mfile_info)
                if strcmp(mfile_info(i).name, field)
                    array_info = mfile_info(i);     
                end
            end

            obj.Size = array_info.size;
            obj.ElementClass = array_info.class;

        end

        function [varargout] = subsref(this, S)

            for i = 1:length(S.subs)
                if isnumeric(S.subs{i})
                    idx{i} = [min(S.subs{i}):max(S.subs{i})];
                    s{i} = S.subs{i} - min(S.subs{i}) + 1;
                else
                    idx{i} = S.subs{i};
                    s{i} = ':';
                end
            end
            
            A = this.Matfile.(this.Field)(idx{:});
            
            varargout{1} = A(s{:});

        end

    end

end
