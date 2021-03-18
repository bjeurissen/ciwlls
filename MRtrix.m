classdef MRtrix
    % Class Wrapper written by Ben Jeurissen, Visionlab, University of Antwerp
    % original read and write functions by the MRtrix developers (www.mrtrix.org)
    %
    
    methods (Access=public,Static=true)
        function image = read(filename)
            image.comments = {};
            
            f = fopen (filename, 'r');
            if (f<1)
                disp (['error opening ' filename ]);
                return
            end
            L = fgetl(f);
            if ~strncmp(L, 'mrtrix image', 12)
                fclose(f);
                disp ([filename ' is not in MRtrix format']);
                return
            end
            
            transform = [];
            dw_scheme = [];
            
            while 1
                L = fgetl(f);
                if ~ischar(L), break, end
                L = strtrim(L);
                if strcmp(L, 'END'), break, end
                d = strfind (L,':');
                if isempty(d)
                    disp (['invalid line in header: ''' L ''' - ignored']);
                else
                    key = lower(strtrim(L(1:d(1)-1)));
                    value = strtrim(L(d(1)+1:end));
                    if strcmp(key, 'dim')
                        image.dim = str2num(char(MRtrix.split_strings (value, ',')))';
                    elseif strcmp(key, 'vox')
                        image.vox = str2num(char(MRtrix.split_strings (value, ',')))';
                    elseif strcmp(key, 'layout')
                        image.layout = value;
                    elseif strcmp(key, 'datatype')
                        image.datatype = value;
                    elseif strcmp(key, 'labels')
                        image.labels = MRtrix.split_strings (value, '\');
                    elseif strcmp(key, 'units')
                        image.units = MRtrix.split_strings (value, '\');
                    elseif strcmp(key, 'transform')
                        transform(end+1,:) = str2num(char(MRtrix.split_strings (value, ',')))';
                    elseif strcmp(key, 'comments')
                        image.comments{end+1} = value;
                    elseif strcmp(key, 'file')
                        file = value;
                    elseif strcmp(key, 'dw_scheme')
                        dw_scheme(end+1,:) = str2num(char(MRtrix.split_strings (value, ',')))';
                        %else
                        %     disp (['unknown key ''' key ''' - ignored']);
                    end
                end
            end
            fclose(f);
            
            
            if ~isempty(transform)
                image.transform = transform;
                image.transform(4,:) = [ 0 0 0 1 ];
            end
            
            if ~isempty(dw_scheme)
                image.dw_scheme = dw_scheme;
            end
            
            if ~isfield (image, 'dim') || ~exist ('file') || ...
                    ~isfield (image, 'layout') || ~isfield (image, 'datatype')
            disp ('critical entries missing in header - not reading data')
            return
            end
            
            layout = MRtrix.split_strings(image.layout, ',');
            order = (abs(str2num (char(layout)))+1)';
            
            [ file, offset ] = strtok(file);
            if isempty(offset), offset = 0; else; offset = str2num(char(offset)); end
            [f,g] = fileparts(filename);
            if strcmp(file,'.'), file = filename; else; file = fullfile (f, file); end
            
            datatype = lower(image.datatype);
            byteorder = datatype(end-1:end);
            
            if strcmp(byteorder, 'le')
                f = fopen (file, 'r', 'l');
                datatype = datatype(1:end-2);
            elseif strcmp(byteorder, 'be')
                f = fopen (file, 'r', 'b');
                datatype = datatype(1:end-2);
            else
                if strcmp(datatype, 'bit')
                    datatype = 'bit1';
                    f = fopen(file, 'r', 'b');
                else
                    f = fopen(file, 'r');
                end
            end
            
            if (f<1)
                disp (['error opening ' filename ]);
                return
            end
            
            fseek (f, offset, -1);
            if strcmp(datatype,'bit1')
                datatype__ = datatype;
            else
                datatype__ = [datatype '=>' datatype];
            end
            
            image.data = fread (f, inf, datatype__);
            
            if strcmp(datatype,'bit1')
                image.data = logical(image.data);
            end
            
            fclose (f);
            
            order(order)= 1:size(order,2);
            
            image.data = reshape (image.data, image.dim(order));
            image.data = ipermute (image.data, order);
            for i=1:size(order,2)
                if layout{i}(1) == '-'
                    image.data = flipdim(image.data, i);
                end
            end
        end
        
        function write(image, filename)
            fid = fopen (filename, 'w');
            fprintf (fid, 'mrtrix image\ndim: ');
            
            if isstruct(image)
                dim = size(image.data);
            else
                dim = size(image);
            end
            fprintf (fid, '%d', dim(1));
            fprintf (fid, ',%d', dim(2:end));
            
            fprintf (fid, '\nvox: ');
            if isstruct (image) && isfield (image, 'vox')
                fprintf (fid, '%f', image.vox(1));
                fprintf (fid, ',%f', image.vox(2:end));
            else
                fprintf(fid, '2');
                fprintf(fid, ',%d', 2*ones(1,size(dim,2)-1));
            end
            
            fprintf (fid, '\nlayout: +0');
            fprintf (fid, ',+%d', 1:(size(dim,2)-1));
            
            [computerType, maxSize, endian] = computer;
            if isstruct (image) && isfield (image, 'datatype')
                datatype = lower(image.datatype);
                byteorder = datatype(end-1:end);
                
                if strcmp (byteorder, 'le')
                    precision = datatype(1:end-2);
                    byteorder = 'l';
                elseif strcmp(byteorder, 'be')
                    precision = datatype(1:end-2);
                    byteorder = 'b';
                else
                    if strcmp(datatype, 'bit')
                        precision = 'bit1';
                        byteorder = 'n';
                    elseif strcmp (datatype, 'int8') || strcmp (datatype, 'uint8')
                        precision = datatype;
                        byteorder = 'n';
                        if endian == 'L'
                            datatype(end+1:end+3) = 'le';
                        else
                            datatype(end+1:end+3) = 'be';
                        end
                    end
                end
            else
                if endian == 'L'
                    datatype = 'float32le';
                else
                    datatype = 'float32be';
                end
                precision = 'float32';
                byteorder = 'n';
            end
            fprintf (fid, [ '\ndatatype: ' datatype ]);
            
            if isstruct (image) && isfield (image, 'comments') && ~iscell(image.comments)
                fprintf (fid, '\ncomments: %s', image.comments);
            end
            
            if isstruct (image) && isfield (image, 'transform')
                fprintf (fid, '\ntransform: %d', image.transform(1,1));
                fprintf (fid, ',%d', image.transform(1,2:4));
                fprintf (fid, '\ntransform: %d', image.transform(2,1));
                fprintf (fid, ',%d', image.transform(2,2:4));
                fprintf (fid, '\ntransform: %d', image.transform(3,1));
                fprintf (fid, ',%d', image.transform(3,2:4));
            end
            
            if isstruct (image) && isfield (image, 'dw_scheme')
                for i=1:size(image.dw_scheme,1)
                    fprintf (fid, '\ndw_scheme: %d', image.dw_scheme(i,1));
                    fprintf (fid, ',%d', image.dw_scheme(i,2:end));
                end
            end
            
            if filename(end-3:end) == '.mif'
                datafile = filename;
                dataoffset = ftell (fid) + 24;
                fprintf (fid, '\nfile: . %d\nEND\n                         ', dataoffset);
            elseif filename(end-3:end) == '.mih'
                datafile = [ filename(end-3:end) '.dat' ];
                dataoffset = 0;
                fprintf (fid, '\nfile: %s %d\nEND\n', datafile, dataoffset);
            else
                disp ('unknown file suffix - aborting');
                return
            end
            
            fclose(fid);
            
            fid = fopen (datafile, 'r+', byteorder);
            fseek (fid, dataoffset, -1);
            
            if isstruct(image)
                fwrite (fid, image.data, precision);
            else
                fwrite (fid, image, precision);
            end
            fclose (fid);
        end
        
        function mif = fromMatlab(data,vox,v2w,grad)
            mif.data = data;
            mif.dim = size(data);
            mif.vox = vox;
            mif.transform = v2w;
            v2v_matlab2mrtrix = eye(4,4); v2v_matlab2mrtrix(1:3,4) = 1;
            mif.transform = mif.transform*v2v_matlab2mrtrix;
            mif.transform = bsxfun(@rdivide,mif.transform,[vox 1]);
            if exist('grad','var')
                mif.dw_scheme = grad;
            end
        end
        
        function [data,dim,vox,v2w,grad] = toMatlab(mif)
            data = mif.data;
            dim = size(data);
            vox = mif.vox(1:3);
            v2w = mif.transform;
            v2w = bsxfun(@times,v2w,[vox 1]);
            v2v_matlab2mrtrix = eye(4,4); v2v_matlab2mrtrix(1:3,4) = -1;
            v2w = v2w*v2v_matlab2mrtrix;
            if nargout > 4
                if isfield(mif,'dw_scheme')
                    grad = mif.dw_scheme;
                else
                    grad = [];
                end
            end
        end
        
        function [data,vox,v2w,grad] = readToMatlab(filename)
            [data,~,vox,v2w,grad] = MRtrix.toMatlab(MRtrix.read(filename));
        end
        
        function writeFromMatlab(data,vox,v2w,grad,filename)
            MRtrix.write(MRtrix.fromMatlab(data,vox,v2w,grad),filename);
        end
        
    end
    
    methods (Access=private,Static=true)
        function S = split_strings (V, delim)
            S = {};
            while size(V,2) > 0
                [R, V] = strtok(V,delim);
                S{end+1} = R;
            end
        end
    end
end


