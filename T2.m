classdef T2 < LogLinear
    % Authors: Ben Jeurissen (ben.jeurissen@uantwerpen.be), Jan Morez (jan.morez@uantwerpen.be)
    %
    % Basic usage:
    %
    % model = T2(te);
    % y = Volumes.mask(y,mask);
    % x = model.solve(y);
    % m = model.metrics(x);
    %
    % with
    %
    % te: n_w × 1 (echo times, preferrably expressed in units that have order of magnitude of 1)
    % y: n_x × n_y × n_z × n_w (weighted image series)
    % mask : n_x × n_y × n_z (boolean processing mask)
    % x: n_x × n_y × n_z × 2 (model parameters)
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
    
    methods (Access = public, Static = false)
        function obj = T2(te, varargin)
            fprintf(1, 'Setting up T2 model ...\n');
            
            % parse T2 specific options
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('constr', 1);
            p.parse(varargin{:});
            
            % set up problem matrix
            te = double(te);
            A = [ones([size(te, 1) 1], class(te)) -te];
            
            % set up constraint matrix
            constr = p.Results.constr;
            Aneq = [];
            if exist('constr', 'var') && any(constr)
                if constr(1)
                    Aneq = [Aneq; 0 -1];
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
        function rho = rho(x)
            rho = exp(x(1,:));
        end
        
        function t2 = t2(x)
            t2 = 1./x(2,:);
        end
    end
    
    methods (Access = public, Static = true)
        function metrics = metrics(x)
            fprintf(1, 'Calculating T2 metrics ...\n');
            if ndims(x) ~= 2; [x, mask] = Volumes.vec(x); end
            metrics.rho = T2.rho(x);
            metrics.t2 = T2.t2(x);
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
