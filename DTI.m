classdef DTI
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
    % grad: n_g × 4 (gradient direction + b-value, preferrably expressed in ms/um^2)
    % y: n_x × n_y × n_z × n_g (DW image series)
    % mask : n_x × n_y × n_z (boolean processing mask)
    % x: n_x × n_y × n_z × 7 (DTI parameters)
    % m: struct with scalar DTI metrics
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
    %    2. The Softwre shall not be distributed to any third parties
    %    without written approval of the authors.
    %   
    %    3. The software is provided 'as is', without warranty of any kind,
    %    express or implied, including but not limited to the warranties of
    %    merchantability, fitness for a particular purpose and noninfringement.
    %    In no event shall the authors or copyright holders be liable for any
    %    claim, damages or other liability, whether in an action of contract,
    %    tort or otherwise, arising from, out of or in connection with the
    %    software or the use or other dealings in the software. 
    %    
    %    4. The Software may only be used for non-commercial research and may
    %    not be used for clinical care. 
    %
    %    5. Before publication by Recipient of research involving the Software
    %    shall contact the authors listed above.
    %
    
    properties (Access = public, Constant = false)
        A
        iter
        estimatorname
        estimator
        init_estimator
    end
    
    methods (Access = public, Static = false)
        function obj = DTI(grad, varargin)
            if nargin > 0
                fprintf(1, 'Setting up DTI model ...\n');
                p = inputParser;
                p.addOptional('estimator', 'wlls');
                p.addOptional('iter', 2);
                p.addOptional('constr', 1);
                p.addOptional('constr_dirs', 100);
                p.parse(varargin{:});
                
                obj.estimatorname = p.Results.estimator;
                if strcmp(obj.estimatorname, 'wlls')
                    obj.iter = p.Results.iter;
                end
                constr = p.Results.constr;
                n = p.Results.constr_dirs;
                
                grad = double(grad);
                grad(:, 1:3) = bsxfun(@rdivide, grad(:, 1:3), sqrt(sum(grad(:, 1:3).^2, 2))); grad(isnan(grad)) = 0;
                obj.A = [ones([size(grad, 1) 1], class(grad)) DTI.grad2A(grad)];
                
                Aneq = [];
                bneq = [];
                if exist('constr', 'var') && any(constr)
                    dirs = Directions.get(n);
                    if constr(1)
                        Aneq = [Aneq; -[zeros(n, 1), prod(reshape(dirs(:,DTI.ind()),[],6,2),3)*diag(DTI.cnt())]];
                    end
                end
                
                if size(Aneq, 1) > 0
                    bneq = zeros(size(Aneq, 1), 1);
                end
                
                switch obj.estimatorname
                    case 'lls'
                        obj.estimator = LLS(obj.A ,Aneq, bneq, [], []);
                    case 'wlls'
                        obj.estimator = LLS(obj.A, Aneq, bneq, [], []);
                    case 'nls'
                        obj.init_estimator = LLS(obj.A, [], [], [], []);
                        obj.estimator = NLS(@obj.ssd, Aneq, bneq, [], []);
                    otherwise
                        error('estimator not supported!')
                end
            end
        end
        
        function x = solve(obj, y, x0)
            if ndims(y) ~= 2; [y, mask] = Volumes.vec(y); if nargin > 2; x0 = Volumes.vec(x0,mask); end; end %#ok<*ISMAT>
            f = 1000/median(y(:));
            y = y.*f;
            y(y<eps) = eps;
            switch obj.estimatorname
                case 'lls'
                    fprintf(1, 'Performing LLS fitting ...\n');
                    x = obj.estimator.solve(log(y));
                case 'wlls'
                    fprintf(1, 'Performing WLLS fitting ...\n');
                    logy = log(y);
                    fprintf(1, 'Initial weighting...\n');
                    x = obj.estimator.solve(logy, y);
                    for it = 1:obj.iter
                        fprintf(1, 'Iterative reweighting #%i...\n', it);
                        x = obj.estimator.solve(logy, obj.predict(x));
                    end
                case 'nls'
                    fprintf(1, 'Performing NLS fitting ...\n');
                    if nargin > 2
                        x0(1, :) = x0(1, :) + log(f);
                        x = obj.estimator.solve(y, x0);
                    else
                        fprintf(1, 'Initial LLS fitting ...\n');
                        x0 = obj.init_estimator.solve(log(y));
                        fprintf(1, 'Final NLS fitting ...\n');
                        warning('off','MATLAB:nearlySingularMatrix');
                        x = obj.estimator.solve(y, x0);
                        warning('on','MATLAB:nearlySingularMatrix');
                    end
            end
            x(1, :) = x(1, :) - log(f);
            if exist('mask','var'); x = Volumes.unvec(x, mask); end
        end
        
        function y = predict(obj, x)
            if ndims(x) ~= 2; [x, mask] = Volumes.vec(x); end
            y = exp(obj.A*x);
            if exist('mask','var'); y = Volumes.unvec(y, mask); end
        end
                
        function [f, g] = ssd(obj, x, y)
            if nargin == 1
                f = size(obj.A,2);
            else
                y_hat = obj.predict(x);
                d = y_hat-y;
                f = sum(d.^2, 1);
                if nargout > 1
                    g = sum(2*obj.A.*(d.*y_hat), 1);
                end
            end
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
            adc = prod(reshape(dir(:,DTI.ind()),[],6,2),3)*diag(DTI.cnt())*x(2:7,:);
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
            if exist('mask','var')
                f = fieldnames(metrics);
                for i = 1:size(f,1)
                    fn = f{i};
                    metrics.(fn) = Volumes.unvec(metrics.(fn), mask);
                end
            end
        end
    end
end