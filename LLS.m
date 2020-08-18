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
                obj.optimopt = optimoptions('quadprog','Display','none','LinearSolver','sparse','MaxIterations',100,'OptimalityTolerance',1e-10,'StepTolerance',1e-10,'ConstraintTolerance',1e-10);
            end
        end
        
        function x = solve(obj,y,w)
            if nargin < 3
                if isempty(obj.Aneq) && isempty(obj.Aeq)
                    x = obj.A\y;
                else
                    x = zeros([size(obj.A,2) size(y,2)]);
                    fprintf(1,'Progress: %3d%%\n',0);
                    for i = 1:size(y,2)
                        if ~mod(i,round(size(y,2)/100))
                            fprintf(1,'\b\b\b\b%3.0f%%',i/size(y,2)*100);
                        end
                        x(:,i) = quadprog(obj.H,-obj.A'*y(:,i),obj.Aneq,obj.bneq,obj.Aeq,obj.beq,[],[],[],obj.optimopt);
                    end
                    fprintf('\n');
                end
            else
                wy = w.*y;
                x = zeros([size(obj.A,2) size(y,2)]);
                fprintf(1,'Progress: %3d%%\n',0);
                for i = 1:size(y,2)
                    if ~mod(i,round(size(y,2)/100))
                        fprintf(1,'\b\b\b\b%3.0f%%',i/size(y,2)*100);
                    end
                    wA = w(:,i).*obj.A;
                    if isempty(obj.Aneq) && isempty(obj.Aeq)
                        x(:,i) = wA\wy(:,i);
                    else
                        x(:,i) = quadprog(wA'*wA,-wA'*wy(:,i),obj.Aneq,obj.bneq,obj.Aeq,obj.beq,[],[],[],obj.optimopt);
                    end
                end
                fprintf('\n');
            end
        end
    end
end