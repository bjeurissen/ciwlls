classdef LLS
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
    
    properties
        A
        Aneq
        bneq
        Aeq
        beq
        optimopt
        H
    end
    
    methods
        function obj = LLS(A,Aneq,bneq,Aeq,beq)
            obj.A = A;
            obj.Aneq = Aneq;
            obj.bneq = bneq;
            obj.Aeq = Aeq;
            obj.beq = beq;
            if ~isempty(Aneq) || ~isempty(Aeq)
                obj.H = obj.A'*obj.A;
                obj.optimopt = optimoptions('quadprog','Display','none','LinearSolver','sparse','MaxIterations',400,'OptimalityTolerance',1e-8,'StepTolerance',1e-12,'ConstraintTolerance',1e-8);
            end
        end
        
        function x = solve(obj,y,w)
            nvox = size(y,2);
            p = 0;
            D = parallel.pool.DataQueue;
            afterEach(D, @nUpdateWaitbar);
            if nargin < 3
                if isempty(obj.Aneq) && isempty(obj.Aeq)
                    x = obj.A\y;
                else
                    x = zeros([size(obj.A,2) nvox]);
                    fprintf(1,'Progress: %3d%%\n',0);
                    parfor i = 1:nvox
                        x(:,i) = quadprog(obj.H,-obj.A'*y(:,i),obj.Aneq,obj.bneq,obj.Aeq,obj.beq,[],[],[],obj.optimopt); %#ok<PFBNS>
                        send(D,i)
                    end
                    fprintf('\n');
                end
            else
                wy = w.*y;
                x = zeros([size(obj.A,2) nvox]);
                fprintf(1,'Progress: %3d%%\n',0);
                parfor i = 1:nvox
                    wA = w(:,i).*obj.A; %#ok<PFBNS>
                    if isempty(obj.Aneq) && isempty(obj.Aeq)
                        x(:,i) = wA\wy(:,i);
                    else
                        x(:,i) = quadprog(wA'*wA,-wA'*wy(:,i),obj.Aneq,obj.bneq,obj.Aeq,obj.beq,[],[],[],obj.optimopt);
                    end
                    send(D,i)
                end
                fprintf('\n');
            end
            function nUpdateWaitbar(~)
                p = p + 1;
                if ~mod(p,round(nvox/100))
                    fprintf(1,'\b\b\b\b%3.0f%%',p/size(y,2)*100);
                end
            end
        end
    end
end
