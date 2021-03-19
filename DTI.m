classdef DTI < LogLinear
    % Authors: Ben Jeurissen (ben.jeurissen@uantwerpen.be), Jan Morez (jan.morez@uantwerpen.be)
    %
    % Basic usage:
    %
    % model = DTI(grad);
    % y = Volumes.mask(y,mask);
    % x = model.solve(y);
    % m = model.metrics(x);
    %
    % with
    %
    % grad: n_w × 4 (gradient direction + b-value, preferrably expressed in ms/um^2)
    % y: n_x × n_y × n_z × n_w (weighted image series)
    % mask : n_x × n_y × n_z (boolean processing mask)
    % x: n_x × n_y × n_z × 7 (model parameters)
    % m: struct with scalar metrics
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
    %    included by the Recipient in all copies or substantial portions of the
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
    
    methods (Access = public, Static = false)
        function obj = DTI(grad, varargin)
            fprintf(1, 'Setting up DTI model ...\n');
            
            % parse DTI specific options
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('constr', 1);
            p.addOptional('constr_dirs', 100);
            p.parse(varargin{:});
            
            % set up problem matrix
            grad = double(grad);
            grad(:, 1:3) = bsxfun(@rdivide, grad(:, 1:3), sqrt(sum(grad(:, 1:3).^2, 2))); grad(isnan(grad)) = 0;
            A = [ones([size(grad, 1) 1], class(grad)) DTI.grad2A(grad)];
            
            % set up constraint matrix
            constr = p.Results.constr;
            n = p.Results.constr_dirs;
            Aneq = [];
            if exist('constr', 'var') && any(constr)
                dirs = Directions.get(n);
                if constr(1)
                    Aneq = [Aneq; [zeros(n, 1) DTI.grad2A(dirs)]];
                end
            end
            
            % set up constraint vector
            bneq = [];
            if size(Aneq, 1) > 0
                bneq = zeros(size(Aneq, 1), 1);
            end
            
            % set up generic y = exp(A*x) problem
            obj = obj@LogLinear(A,Aneq,bneq,varargin{:});
        end
    end
    
    methods (Access = private, Static = true)
        function b0 = b0(x)
            b0 = exp(x(1,:));
        end
        
        function fa = fa(eigval)
            l1 = eigval(1,:); l2 = eigval(2,:); l3 = eigval(3,:);
            fa = sqrt(1/2).*sqrt((l1-l2).^2+(l2-l3).^2+(l3-l1).^2)./sqrt(l1.^2+l2.^2+l3.^2);
        end
        function colfa = colfa(eigval,eigvec)
            colfa = abs(eigvec(1:3,:)).*DTI.fa(eigval);
        end
        
        function ad = ad(eigval)
            ad = eigval(1,:);
        end
        
        function rd = rd(eigval)
            rd = mean(eigval(2:3,:),1);
        end
        
        function md = md(eigval)
            md = mean(eigval,1);
        end
    end
    
    methods (Access = public, Static = true)
        function v = ind()
            v = [1 1; 1 2; 1 3; 2 2; 2 3; 3 3];
        end
        
        function v = cnt()
            v = [1 2 2 1 2 1];
        end
        
        function A = grad2A(grad)
            if size(grad,2) < 4
                grad(:,4) = 1;
            end
            A = -(grad(:, 4)).*prod(reshape(grad(:,DTI.ind()),[],6,2),3)*diag(DTI.cnt());
        end
        
        function [eigval, eigvec] = eig(x,k)
            if nargin < 2
                k = 3;
            end
            t = x(2:7,:);
            eigval = zeros(k,size(t,2));
            eigvec = zeros(3*k,size(t,2));
            for i = 1:size(t,2)
                ti = reshape(t([1 2 3 2 4 5 3 5 6],i), [3 3]);
                [vec, val] = eigs(ti,k);
                eigval(:,i) = diag(val);
                eigvec(:,i) = vec(:);
            end
        end
        
        function adc = adc(x, dir)
            adc = -DTI.grad2A(dir)*x(2:7,:);
        end
        
        function metrics = metrics(x)
            fprintf(1, 'Calculating DTI metrics ...\n');
            if ndims(x) ~= 2; [x, mask] = Volumes.vec(x); end
            metrics.b0 = DTI.b0(x);
            [metrics.eigval, metrics.eigvec] = DTI.eig(x);
            metrics.ad = DTI.ad(metrics.eigval);
            metrics.rd = DTI.rd(metrics.eigval);
            metrics.md = DTI.md(metrics.eigval);
            metrics.fa = DTI.fa(metrics.eigval);
            metrics.colfa = DTI.colfa(metrics.eigval,metrics.eigvec);
            if exist('mask','var'); metrics = Volumes.unvec_struct(metrics); end
        end
    end
end
