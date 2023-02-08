classdef Directions
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

    methods (Access = public, Static = true)
        function S = c2s(C)
            r = sqrt(sum(C.^2, 2));
            S(:,1) = atan2(C(:,2), C(:,1));
            S(:,2) = acos(C(:,3)./r);
        end

        function C = s2c(S)
            C = [ sin(S(:,2)).*cos(S(:,1)), ...
                sin(S(:,2)).*sin(S(:,1)), ...
                cos(S(:,2)) ];
        end

        function p = equatorpoints(d, k)
            dt = pi/k;
            theta = 0:dt:(pi-dt);
            p = [cos(theta); sin(theta); zeros(size(theta))];
            R = vrrotvec2mat(vrrotvec([0 0 1],d));
            p = R*p;
        end

        function dirs = get(n)
            fname = fullfile(mfilename('fullpath'),sprintf('dir%04i.txt',n));
            try
                dirs = TextFile.read(fname);
            catch
                error('Number of directions not supported')
            end
        end
    end
end
