classdef NLS
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
    
    properties
        fun
        Aneq
        bneq
        Aeq
        beq
        optimopt
    end
    
    methods
        function obj = NLS(fun,Aneq,bneq,Aeq,beq)
            obj.fun = fun;
            obj.Aneq = Aneq;
            obj.bneq = bneq;
            obj.Aeq = Aeq;
            obj.beq = beq;
            if isempty(Aneq) && isempty(Aeq)
                obj.optimopt = optimoptions('fminunc','Display','none','SpecifyObjectiveGradient',true,'MaxIterations',400,'MaxFunctionEvaluations',inf,'OptimalityTolerance',1e-8,'StepTolerance',1e-12);
            else
                obj.optimopt = optimoptions('fmincon','Display','none','SpecifyObjectiveGradient',true,'MaxIterations',400,'MaxFunctionEvaluations',inf,'OptimalityTolerance',1e-8,'StepTolerance',1e-12,'ConstraintTolerance',1e-8);
            end
        end
        
        function x = solve(obj,y,x0)
            warning('off','MATLAB:nearlySingularMatrix');
            nvox = size(y,2);
            p = 0;
            D = parallel.pool.DataQueue;
            afterEach(D, @nUpdateWaitbar);
            if nargin > 2
                x = x0;
            else
                x = zeros([obj.fun() nvox]);
            end
            fprintf(1,'Progress: %3d%%\n',0);
            if isempty(obj.Aneq) && isempty(obj.Aeq)
                parfor i = 1:nvox
                    try
                        x(:,i) = fminunc(@(x) obj.fun(x,y(:,i)),x(:,i),obj.optimopt);
                    catch
                        x(:,i) = NaN;
                        warning('NLS:solve could not recover from NaN or Inf');
                    end
                    send(D,i)
                end
            else
                parfor i = 1:nvox
                    x(:,i) = fmincon(@(x) obj.fun(x,y(:,i)),x(:,i),obj.Aneq,obj.bneq,obj.Aeq,obj.beq,[],[],[],obj.optimopt);
                    send(D,i)
                end
            end
            fprintf('\n');
            function nUpdateWaitbar(~)
                p = p + 1;
                if ~mod(p,round(nvox/100))
                    fprintf(1,'\b\b\b\b%3.0f%%',p/nvox*100);
                end
            end
            warning('on','MATLAB:nearlySingularMatrix');
        end
    end
end
