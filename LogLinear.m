classdef (Abstract) LogLinear
    % Authors: Ben Jeurissen (ben.jeurissen@uantwerpen.be), Jan Morez (jan.morez@uantwerpen.be)
    %
    % Basic usage:
    %
    % model = LogLinear(A,Aneq,bneq,params);
    % y = Volumes.mask(y,mask);
    % x = model.solve(y);
    % m = model.metrics(x);
    %
    % with
    %
    % A: n_w × n_p (relates x to y as y = exp(A*x))
    % Aneq: n_c × n_p (inequality constraint matrix)
    % bneq: n_c × 1 (unequality constraint vector)
    % y: n_x × n_y × n_z × n_w (weighted image series)
    % mask : n_x × n_y × n_z (boolean processing mask)
    % x: n_x × n_y × n_z × n_p (model parameters)
    % m: struct with scalar model metrics
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
        init_weight
        init_estimator
    end
    
    methods (Access = public, Static = false)
        function obj = LogLinear(A, Aneq, bneq, varargin)
            fprintf(1, 'Setting up generic exp(A*x) model ...\n');
            obj.A = A;
            p = inputParser;
            p.addOptional('estimator', 'wlls');
            p.addOptional('init_estimator', 'lls');
            p.addOptional('init_weight', 'data');
            p.addOptional('iter', 2);
            p.parse(varargin{:});
            
            obj.estimatorname = p.Results.estimator;
            obj.init_estimator = p.Results.init_estimator;
            
            if strcmp(obj.estimatorname, 'wlls')
                obj.init_weight = p.Results.init_weight;
                obj.iter = p.Results.iter;
            end
            
            switch obj.estimatorname
                case 'lls'
                    obj.estimator = LLS(obj.A, Aneq, bneq, [], []);
                case 'wlls'
                    obj.estimator = LLS(obj.A, Aneq, bneq, [], []);
                case 'nls'
                    obj.init_estimator = LLS(obj.A, [], [], [], []);
                    obj.estimator = NLS(@obj.ssd, Aneq, bneq, [], []);
                otherwise
                    error('estimator not supported!')
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
                    switch obj.init_weight
                        case 'data'
                            fprintf(1, 'Initial weighting using data...\n');
                            x = obj.estimator.solve(logy, y);
                        case 'ones'
                            fprintf(1, 'Initial weighting using ones...\n');
                            x = obj.estimator.solve(logy, ones(size(y)));
                    end
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
                        x = obj.estimator.solve(y, x0);
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
end
