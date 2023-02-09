classdef QTI < LogLinear
    % Authors: Ben Jeurissen (ben.jeurissen@uantwerpen.be), Jan Morez (jan.morez@uantwerpen.be)
    %
    % Basic usage:
    %
    % model = QTI(grad);
    % y = Volumes.mask(y,mask);
    % x = model.solve(y);
    % m = model.metrics(x);
    %
    % with
    %
    % grad: n_w × 5 (gradient direction + b-value + b-delta, preferrably expressed in ms/um^2)
    % y: n_x × n_y × n_z × n_w (weighted image series)
    % mask : n_x × n_y × n_z (boolean processing mask)
    % x: n_x × n_y × n_z × 28 (model parameters)
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
        function obj = QTI(grad, varargin)
            fprintf(1, 'Setting up QTI model...\n');

            % parse QTI specifc options
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('constr', [0 0 1 1 1]);
            p.addOptional('constr_dirs', 100);
            p.addOptional('rank23', false);
            p.parse(varargin{:});

            % set up problem matrix
            grad = double(grad);
            grad(:, 1:3) = bsxfun(@rdivide, grad(:, 1:3), sqrt(sum(grad(:, 1:3).^2, 2))); grad(isnan(grad)) = 0;
            A = [ones([size(grad, 1) 1], class(grad)) -Tensor.tpars_to_1x6(grad(:,4), grad(:,5), grad(:,1:3)) 0.5*Tensor.t_1x6_to_1x21(Tensor.tpars_to_1x6(grad(:,4), grad(:,5), grad(:,1:3)))];

            % set up constraint matrix
            constr = p.Results.constr;
            n = p.Results.constr_dirs;
            Aneq = [];
            if exist('constr', 'var') && any(constr)
                dirs = Directions.get(n);
                c1 = Tensor.tpars_to_1x6(ones(n,1), ones(n,1), dirs); % linear
                c2 = Tensor.t_1x6_to_1x21(c1);
                E_bulk = Tensor.iso_1x21();
                if constr(1)
                    disp('Constraining to non-negative diffusivity')
                    Aneq = [Aneq; -[zeros(n, 1) c1 zeros(n, 21)]];
                end
                if constr(2)
                    disp('Constraining to non-negative total kurtosis')
                    Aneq = [Aneq; -[zeros(n, 7) c2]];
                end
                if constr(3)
                    disp('Constraining to non-negative isotropic kurtosis')
                    Aneq = [Aneq; -[zeros(1, 7) E_bulk]];
                end
                if constr(4)
                    disp('Constraining to non-negative anisotropic kurtosis')
                    Aneq = [Aneq; -[zeros(n, 7) c2-E_bulk]];
                end
                if constr(5)
                    disp('Constraining to monotonic signal decay')
                    Aneq = [Aneq; [zeros(n, 1) -c1 max(grad(:, 4))*c2]];
                end
            end

            % set up constraint vector
            bneq = [];
            if size(Aneq, 1) > 0
                bneq = zeros(size(Aneq, 1), 1);
            end

            % add equality constraints if rank23 is true (needed to support LTE+STE only data)
            Aeq = []; beq = [];
            if p.Results.rank23
                disp('Constraining parameters 12,13,14,15, and 16 to zero.');
                error('Currently the rank23 option is not supported... Please wait for an update....')
                Aeq = zeros(5,28);
                Aeq(1,12) = 1;
                Aeq(2,13) = 1;
                Aeq(3,14) = 1;
                Aeq(4,15) = 1;
                Aeq(5,16) = 1;
                beq = zeros(5,1);
            end

            r = rank(A,1e-10);
            if r < size(A,2)
                if r >= 23
                    if ~p.Results.rank23
                        error('A in exp(A*x) is not full rank. Use "QTI(..., ''rank23'', true)" to support LTE+STE only data.')
                    end
                else
                    error('A in exp(A*x) is not full rank.');
                end
            else
                if p.Results.rank23
                    disp('WARNING: Constraining parameters 12,13,14,15, and 16 to zero, despite full rank A.');
                end
            end

            % set up generic y = exp(A*x) problem
            obj = obj@LogLinear(A,Aneq,bneq,Aeq,beq,varargin{:});
        end
    end

    methods (Access = public, Static = true)
        function metrics = metrics(x)
            if ndims(x) ~= 2; [x, mask] = Volumes.vec(x); end %#ok<ISMAT>
            fprintf(1, 'Calculating QTI metrics ...\n');
            dt_1x6 = x(2:7,:)';
            B0 = exp(x(1,:));
            FA = Tensor.fa(dt_1x6);
            MD = Tensor.md(dt_1x6);
            L = Tensor.eigval(dt_1x6);
            AD = real(L(:,1));
            RD = mean(real(L(:,2:3)),2);

            dt2_1x21 = Tensor.t_1x6_to_1x21(dt_1x6);
            ct_1x21 = x(8:28,:)';
            [E_bulk, E_shear, E_iso] = Tensor.iso_1x21();

            V_MD2    = Tensor.inner(dt2_1x21, E_bulk);
            V_shear2 = Tensor.inner(dt2_1x21, E_shear);
            V_iso2   = Tensor.inner(dt2_1x21, E_iso);

            V_MD     = Tensor.inner(ct_1x21 , E_bulk);
            V_shear  = Tensor.inner(ct_1x21 , E_shear);
            V_iso    = Tensor.inner(ct_1x21 , E_iso);

            V_MD1    = V_MD    + V_MD2;
            V_shear1 = V_shear + V_shear2;
            V_iso1   = V_iso   + V_iso2;

            C_MD =       V_MD     ./ V_MD1;
            C_mu = 1.5 * V_shear1 ./ V_iso1;
            C_M  = 1.5 * V_shear2 ./ V_iso2;
            C_c  =            C_M ./ C_mu;

            MKi  = 3 * V_MD ./ V_MD2;
            MKa  = (6/5) * V_shear1 ./ V_MD2;
            MKad = (6/5) * V_shear  ./ V_MD2;
            MKt  = MKi + MKa;
            MK   = MKad + MKi;
            MKd  = MKa - MKad;
            uFA  = sqrt(C_mu);

            S_I = sqrt(V_MD.*(V_MD > 0));
            S_A = sqrt(V_shear1.*(V_shear1 > 0));

            metrics.b0 = B0;
            metrics.fa = FA';
            metrics.md = MD';
            metrics.ad = AD';
            metrics.rd = RD';

            metrics.v_md2 = V_MD2';
            metrics.v_shear2 = V_shear2';
            metrics.v_iso2 = V_iso2';

            metrics.v_md = V_MD';
            metrics.v_shear = V_shear';
            metrics.v_iso = V_iso';

            metrics.v_md1 = V_MD1';
            metrics.v_shear1 = V_shear1';
            metrics.v_iso1 = V_iso1';

            metrics.c_md = C_MD';
            metrics.c_mu = C_mu';
            metrics.c_m = C_M';
            metrics.c_c = C_c';

            metrics.mki = MKi';
            metrics.mka = MKa';
            metrics.mkad = MKad';
            metrics.mkt = MKt';
            metrics.mk = MK';
            metrics.mkd = MKd';
            metrics.ufa = uFA';

            metrics.s_i = S_I';
            metrics.s_a = S_A';
            if exist('mask','var'); metrics = Volumes.unvec_struct(metrics,mask); end
        end
    end
end
