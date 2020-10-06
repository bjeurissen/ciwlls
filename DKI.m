classdef DKI < DTI
    % Authors: Ben Jeurissen (ben.jeurissen@uantwerpen.be), Jan Morez (jan.morez@uantwerpen.be)
    %
    % Basic usage:
    %
    % model = DKI(grad);
    % y = Volumes.mask(y,mask);
    % x = model.solve(y);
    % m = model.metrics(x);
    %
    % with
    %
    % grad: n_g × 4 (gradient direction + b-value, preferrably expressed in ms/um^2)
    % y: n_x × n_y × n_z × n_g (DW image series)
    % mask : n_x × n_y × n_z (boolean processing mask)
    % x: n_x × n_y × n_z × 22 (DKI parameters)
    % m: struct with scalar DKI metrics
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
    
    methods (Access = public, Static = false)
        function obj = DKI(grad, varargin)
            if nargin > 0
                fprintf(1, 'Setting up DKI model...\n');
                p = inputParser;
                p.addOptional('estimator', 'wlls');
                p.addOptional('iter', 2);
                p.addOptional('constr', [1 0 1]);
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
                obj.A = [ones([size(grad, 1) 1], class(grad)) DTI.grad2A(grad) DKI.grad2A(grad)];
                
                Aneq = [];
                bneq = [];
                if exist('constr', 'var') && any(constr)
                    dirs = Directions.get(n);
                    if constr(1)
                        Aneq = [Aneq; [zeros(n, 1) DTI.grad2A(dirs) zeros(n, 15)]];
                    end
                    if constr(2)
                        Aneq = [Aneq; [zeros(n, 7) -DKI.grad2A(dirs)]];
                    end
                    if constr(3)
                        Aneq = [Aneq; [zeros(n, 1) DTI.grad2A(dirs) 2*max(grad(:, 4))*DKI.grad2A(dirs)]];
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
    end
    
    methods (Access = private, Static = true)
        function [ak,rk,mk] = armk(x)
            ak = zeros(1,size(x,2));
            rk = zeros(1,size(x,2));
            mk = zeros(1,size(x,2));
            [~, pdv] = DTI.eig(x,1);
            for i = 1:size(x,2)
                ak(1,i) = DKI.akc(x(:,i),pdv(:,i)');
                rk(1,i) = mean(DKI.akc(x(:,i),Directions.equatorpoints(pdv(:,i),300)'),1);
                mk(1,i) = mean(DKI.akc(x(:,i),Directions.get(300)),1);
            end
        end
        
        function akc = akc(x, dir)
            adc = DTI.adc(x, dir);
            akc = 6*DKI.grad2A(dir)*x(8:22, :);
            akc = akc./adc.^2;
        end
    end
    
    methods (Access = public, Static = true)
        function v = ind()
            v =[1 1 1 1; 1 1 1 2; 1 1 1 3; 1 1 2 2; 1 1 2 3; 1 1 3 3; 1 2 2 2; 1 2 2 3; 1 2 3 3; 1 3 3 3; 2 2 2 2; 2 2 2 3; 2 2 3 3; 2 3 3 3; 3 3 3 3];
        end
        
        function v = cnt()
            v = [1 4 4 6 12 6 4 12 12 4 1 4 6 4 1];
        end
        
        function A = grad2A(grad)
            if size(grad,2) < 4
                grad(:,4) = 1;
            end
            bsqd6 = (grad(:, 4).^2)/6;
            A = (bsqd6 .* prod(reshape(grad(:, DKI.ind()), [], 15, 4), 3))*diag(DKI.cnt());
        end
        
        function metrics = metrics(x)
            if ndims(x) ~= 2; [x, mask] = Volumes.vec(x); end
            metrics = DTI.metrics(x);
            fprintf(1, 'Calculating DKI metrics ...\n');
            [metrics.ak, metrics.rk, metrics.mk] = DKI.armk(x);
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
