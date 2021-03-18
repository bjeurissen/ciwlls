classdef Volumes
    % Authors: Ben Jeurissen (ben.jeurissen@uantwerpen.be)
    %
    %
    %    Copyright (c) 2020 University of Antwerp
    %
    %    Permission is hereby granted, free of charge, to any non-commercial
    %    entity ('Recipient') obtaining a copy of this software and associated
    %    documentation files (the 'Software'), to the Software solely for
    %    non-commercial research, including the rights to use, copy and modify
    %    the Software, subject to the following conditions:
    %
    %    1. The above copyright notice and this permission notice shall be
    %    included by Recipient in all copies or substantial portions of the
    %    Software.
    %
    %    2. The Software shall not be distributed to any third parties
    %    without written approval of the authors.
    %
    %    3. The Software is provided 'as is', without warranty of any kind,
    %    express or implied, including but not limited to the warranties of
    %    merchantability, fitness for a particular purpose and noninfringement.
    %    In no event shall the authors or copyright holders be liable for any
    %    claim, damages or other liability, whether in an action of contract,
    %    tort or otherwise, arising from, out of or in connection with the
    %    Software or the use or other dealings in the Software.
    %
    %    4. The Software may only be used for non-commercial research and may
    %    not be used for clinical care.
    %
    %    5. Prior to publication of research involving the Software, the
    %    Recipient shall inform the Authors listed above.
    %
    
    methods (Access=public,Static=true)
        function [y, mask] = vec(x,mask)
            if ~exist('mask','var')
                mask = ~isnan(x(:,:,:,1));
            end
            if islogical(x)
                y = false([size(x,4) sum(mask(:))]);
            else
                y = zeros([size(x,4) sum(mask(:))], class(x));
            end
            for k = 1:size(x,4)
                Dummy = x(:,:,:,k);
                y(k,:) = Dummy(mask(:));
            end
        end
        
        function y = unvec(x,mask)
            dims = [size(mask,1) size(mask,2) size(mask,3)];
            if isfloat(x)
                y = NaN([dims(1) dims(2) dims(3) size(x,1)], class(x));
            elseif islogical(x)
                y = false([dims(1) dims(2) dims(3) size(x,1)]);
            else
                y = zeros([dims(1) dims(2) dims(3) size(x,1)], class(x));
            end
            
            for k = 1:size(x,1)
                if isfloat(x)
                    Dummy = NaN(dims, class(x));
                elseif islogical(x)
                    Dummy = false(dims);
                else
                    Dummy = zeros(dims, class(x));
                end
                Dummy(mask) = x(k,:);
                y(:,:,:,k) = Dummy;
            end
        end
        
        function y = unvec0(x,mask)
            dims = [size(mask,1) size(mask,2) size(mask,3)];
            if islogical(x)
                y = false([dims(1) dims(2) dims(3) size(x,1)]);
            else
                y = zeros([dims(1) dims(2) dims(3) size(x,1)], class(x));
            end
            
            for k = 1:size(x,1)
                if islogical(x)
                    Dummy = false(dims, class(x));
                else
                    Dummy = zeros(dims, class(x));
                end
                Dummy(mask) = x(k,:);
                y(:,:,:,k) = Dummy;
            end
        end
        
        function y = mask(x,mask)
            y = x;
            nanmask = ~mask & true([1 1 1 size(y,4)]);
            y(nanmask) = NaN;
        end
    end
end
