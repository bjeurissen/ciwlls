classdef Tensor
    methods(Static)
        function alpha = rot_6x6(R)
            alpha = [
                R(1,1).^2             R(1,2).^2             R(1,3).^2             sqrt(2)*R(1,1)*R(1,2)       sqrt(2)*R(1,3)*R(1,1)       sqrt(2)*R(1,2)*R(1,3)      ;
                R(2,1).^2             R(2,2).^2             R(2,3).^2             sqrt(2)*R(2,1)*R(2,2)       sqrt(2)*R(2,3)*R(2,1)       sqrt(2)*R(2,2)*R(2,3)      ;
                R(3,1).^2             R(3,2).^2             R(3,3).^2             sqrt(2)*R(3,1)*R(3,2)       sqrt(2)*R(3,3)*R(3,1)       sqrt(2)*R(3,2)*R(3,3)      ;
                sqrt(2)*R(1,1)*R(2,1) sqrt(2)*R(1,2)*R(2,2) sqrt(2)*R(1,3)*R(2,3) R(1,2)*R(2,1)+R(2,2)*R(1,1) R(1,3)*R(2,1)+R(2,3)*R(1,1) R(1,2)*R(2,3)+R(2,2)*R(1,3);
                sqrt(2)*R(3,1)*R(1,1) sqrt(2)*R(3,2)*R(1,2) sqrt(2)*R(3,3)*R(1,3) R(3,2)*R(1,1)+R(1,2)*R(3,1) R(3,3)*R(1,1)+R(1,3)*R(3,1) R(3,2)*R(1,3)+R(1,2)*R(3,3);
                sqrt(2)*R(2,1)*R(3,1) sqrt(2)*R(2,2)*R(3,2) sqrt(2)*R(2,3)*R(3,3) R(2,2)*R(3,1)+R(3,2)*R(2,1) R(2,3)*R(3,1)+R(3,3)*R(2,1) R(2,2)*R(3,3)+R(3,2)*R(2,3);
                ];
        end
        
        function E_iso = iso_3x3()
            E_iso = eye(3)/3;
        end
        
        function [E_bulk, E_shear, E2_iso] = iso_6x6()
            E_iso = Tensor.t_3x3_to_1x6(Tensor.iso_3x3());
            E_bulk  = Tensor.outer(E_iso, E_iso);
            E_shear = eye(6)/3 - E_bulk;
            E2_iso  = eye(6)/3;
        end
        
        
        function [E_bulk, E_shear, E2_iso] = iso_1x21()
            [E_bulk, E_shear, E2_iso] = Tensor.iso_6x6();
            E_bulk  = Tensor.t_6x6_to_1x21(E_bulk);
            E_shear = Tensor.t_6x6_to_1x21(E_shear);
            E2_iso  = Tensor.t_6x6_to_1x21(E2_iso);
        end
        
        function t_1x6 = t_3x3_to_1x6(t_3x3)
            t_1x6 = t_3x3([1 5 9 2 3 6]) .* [1 1 1 sqrt(2) sqrt(2) sqrt(2)];
        end
        
        function t_3x3 = t_1x6_to_3x3(t_1x6)
            assert(size(t_1x6,2) == 6)
            t = t_1x6 .* [1 1 1 sqrt(1/2) sqrt(1/2) sqrt(1/2)];
            t_3x3 = t([1 4 5; 4 2 6; 5 6 3]);
        end
        
        function t_1x6 = v_1x3_to_1x6(v_1x3,ad,rd)
            assert(size(v_1x3,2) == 3)
            if nargin < 3
                rd = 0;
            end
            if nargin < 2
                ad = 1;
            end
            x = v_1x3(:,1);
            y = v_1x3(:,2);
            z = v_1x3(:,3);
            o = zeros([1 3]);
            e = ones([1 3]);
            c = sqrt(2);
            t_1x6 = [x.*x y.*y z.*z c.*x.*y c.*x.*z c.*y.*z] .* (ad - rd) + [e o] .* rd;
        end
        
        function t_1x15 = v_1x3_to_1x15(v_1x3)
            assert(size(v_1x3,2) == 3)
            t_1x15 = [...
                v_1x3(:,1) .* v_1x3(:,1) .* v_1x3(:,1) .* v_1x3(:,1) * sqrt(1) ...
                v_1x3(:,1) .* v_1x3(:,1) .* v_1x3(:,1) .* v_1x3(:,2) * sqrt(4) ...
                v_1x3(:,1) .* v_1x3(:,1) .* v_1x3(:,1) .* v_1x3(:,3) * sqrt(4) ...
                v_1x3(:,1) .* v_1x3(:,1) .* v_1x3(:,2) .* v_1x3(:,2) * sqrt(6) ...
                v_1x3(:,1) .* v_1x3(:,1) .* v_1x3(:,2) .* v_1x3(:,3) * sqrt(12) ...
                v_1x3(:,1) .* v_1x3(:,1) .* v_1x3(:,3) .* v_1x3(:,3) * sqrt(6) ...
                v_1x3(:,1) .* v_1x3(:,2) .* v_1x3(:,2) .* v_1x3(:,2) * sqrt(4) ...
                v_1x3(:,1) .* v_1x3(:,2) .* v_1x3(:,2) .* v_1x3(:,3) * sqrt(12) ...
                v_1x3(:,1) .* v_1x3(:,2) .* v_1x3(:,3) .* v_1x3(:,3) * sqrt(12) ...
                v_1x3(:,1) .* v_1x3(:,3) .* v_1x3(:,3) .* v_1x3(:,3) * sqrt(4) ...
                v_1x3(:,2) .* v_1x3(:,2) .* v_1x3(:,2) .* v_1x3(:,2) * sqrt(1) ...
                v_1x3(:,2) .* v_1x3(:,2) .* v_1x3(:,2) .* v_1x3(:,3) * sqrt(4) ...
                v_1x3(:,2) .* v_1x3(:,2) .* v_1x3(:,3) .* v_1x3(:,3) * sqrt(6) ...
                v_1x3(:,2) .* v_1x3(:,3) .* v_1x3(:,3) .* v_1x3(:,3) * sqrt(4) ...
                v_1x3(:,3) .* v_1x3(:,3) .* v_1x3(:,3) .* v_1x3(:,3) * sqrt(1) ];
        end
        
        
        function t = t_1x15_to_1x21(N)
            N = N ./ repmat(sqrt([1 4 4 6 12 6 4 12 12 4 1 4 6 4 1]), size(N,1), 1);
            
            t = [...
                [N(:,1)  N(:,11) N(:,15)] * sqrt(1) ...
                [N(:,4)  N(:,6)  N(:,13)] * sqrt(2) ...
                [N(:,2)  N(:,8)  N(:,14)] * sqrt(4) ...
                [N(:,5)  N(:,3)  N(:,12)] * sqrt(4) ...
                [N(:,7)  N(:,10) N(:,9)]  * sqrt(4) ...
                [N(:,13) N(:,6)  N(:,4)]  * sqrt(4) ...
                [N(:,9)  N(:,8)  N(:,5)]  * sqrt(8) ];
            
        end
        
        function d_1x6 = tpars_to_1x6(trace, delta, uvec)
            assert(size(uvec,2) == 3)
            t_stick  = Tensor.v_1x3_to_1x6(uvec, 1, 0);
            t_sphere = Tensor.t_3x3_to_1x6(Tensor.iso_3x3());
            c2 = trace .* delta;
            c1 = trace - c2;
            d_1x6 = c1 .* t_sphere + c2 .* t_stick;
        end
        
        
        function t_1x21 = t_1x6_to_1x21(t_1x6)
            assert(size(t_1x6,2) == 6)
            xx = t_1x6(:,1);
            yy = t_1x6(:,2);
            zz = t_1x6(:,3);
            xy = t_1x6(:,4);
            xz = t_1x6(:,5);
            yz = t_1x6(:,6);
            t_1x21 = [...
                [xx.*xx yy.*yy zz.*zz] * sqrt(1) ...
                [xx.*yy xx.*zz yy.*zz] * sqrt(2) ...
                [xx.*yz yy.*xz zz.*xy] * sqrt(2) ...
                [xx.*xy xx.*xz]        * sqrt(2) ...
                [yy.*xy yy.*yz]        * sqrt(2) ...
                [zz.*xz zz.*yz]        * sqrt(2) ...
                [xy.*xy xz.*xz yz.*yz] * sqrt(1) ...
                [xy.*xz xy.*yz xz.*yz] * sqrt(2) ];
            %             % equivalent:
            %             t_1x21 = zeros([size(t_1x6,1) 21]);
            %             for i = 1:size(t_1x6,1)
            %                 t_1x21(i,:) = Tensor.t_6x6_to_1x21(Tensor.outer(t_1x6(i,:),t_1x6(i,:)));
            %             end
        end
        
        function t_1x21 = t_6x6_to_1x21(t_6x6)
            assert(all(size(t_6x6) == [6 6]))
            t_1x21 = [...
                [t_6x6(1,1) t_6x6(2,2) t_6x6(3,3)] * sqrt(1) ...
                [t_6x6(1,2) t_6x6(1,3) t_6x6(2,3)] * sqrt(2) ...
                [t_6x6(1,6) t_6x6(5,2) t_6x6(4,3)] * sqrt(2) ... % xxyz yyxz zzxy
                [t_6x6(1,4) t_6x6(1,5) t_6x6(2,4)] * sqrt(2) ...
                [t_6x6(2,6) t_6x6(3,5) t_6x6(3,6)] * sqrt(2) ...
                [t_6x6(4,4) t_6x6(5,5) t_6x6(6,6)] * sqrt(1) ...
                [t_6x6(4,5) t_6x6(4,6) t_6x6(5,6)] * sqrt(2)];
        end
        
        function t_6x6 = t_1x21_to_6x6(t_1x21)
            % function t = tm_1x21_to_6x6(t)
            %
            % Convert fourth-order tensor in 1x21 format to 6x6 format
            assert(size(t_1x21,2) == 21)
            xxxx = t_1x21(1); % all same
            yyyy = t_1x21(2);
            zzzz = t_1x21(3);
            
            c = sqrt( 2 );
            xxyy = t_1x21(4) / c;
            xxzz = t_1x21(5) / c;
            yyzz = t_1x21(6) / c;
            
            c = sqrt( 2 );
            xxyz = t_1x21(7) / c;
            yyxz = t_1x21(8) / c;
            zzxy = t_1x21(9) / c;
            
            c = sqrt( 2 );
            xxxy = t_1x21(10) / c;
            xxxz = t_1x21(11) / c;
            yyxy = t_1x21(12) / c;
            yyyz = t_1x21(13) / c;
            zzxz = t_1x21(14) / c;
            zzyz = t_1x21(15) / c;
            
            c = 1;
            xyxy = t_1x21(16) / c;
            xzxz = t_1x21(17) / c;
            yzyz = t_1x21(18) / c;
            
            c = sqrt(2);
            xyxz = t_1x21(19) / c;
            xyyz = t_1x21(20) / c;
            xzyz = t_1x21(21) / c;
            
            
            A = [...
                xxxx xxyy xxzz;
                xxyy yyyy yyzz;
                xxzz yyzz zzzz]; % 6 unique
            
            B = [...
                xxxy xxxz xxyz;
                yyxy yyxz yyyz;
                zzxy zzxz zzyz]; % 9 unique
            
            C = [...
                xyxy xyxz xyyz;
                xyxz xzxz xzyz;
                xyyz xzyz yzyz];
            
            % this is fishy, but probably correct given the assumptions on cross terms
            % w = sqrt(ones(3,3)*2 + eye(3)*2);
            % w = ones(3,3);
            
            t_6x6 = [A B; B' C];
        end
        
        function t_3x3x3x3 = t_1x21_to_3x3x3x3(t_1x21)
            assert(size(t_1x21,2)==21)
            xxxx = t_1x21(1); % all same
            yyyy = t_1x21(2);
            zzzz = t_1x21(3);
            
            c = sqrt(2);
            xxyy = t_1x21(4) / c;
            xxzz = t_1x21(5) / c;
            yyzz = t_1x21(6) / c;
            
            c = sqrt(4);
            xxyz = t_1x21(7) / c;
            yyxz = t_1x21(8) / c;
            zzxy = t_1x21(9) / c;
            
            c = sqrt(4);
            xxxy = t_1x21(10) / c;
            xxxz = t_1x21(11) / c;
            yyxy = t_1x21(12) / c;
            yyyz = t_1x21(13) / c;
            zzxz = t_1x21(14) / c;
            zzyz = t_1x21(15) / c;
            
            c = sqrt(4);
            xyxy = t_1x21(16) / c;
            xzxz = t_1x21(17) / c;
            yzyz = t_1x21(18) / c;
            
            c = sqrt(8);
            xyxz = t_1x21(19) / c;
            xyyz = t_1x21(20) / c;
            xzyz = t_1x21(21) / c;
            
            t_3x3x3x3 = zeros(3,3,3,3);
            
            % 9
            t_3x3x3x3(1,1,1,1) = xxxx;
            t_3x3x3x3(1,1,1,2) = xxxy;
            t_3x3x3x3(1,1,1,3) = xxxz;
            
            t_3x3x3x3(1,1,2,1) = xxxy;
            t_3x3x3x3(1,1,2,2) = xxyy;
            t_3x3x3x3(1,1,2,3) = xxyz;
            
            t_3x3x3x3(1,1,3,1) = xxxz;
            t_3x3x3x3(1,1,3,2) = xxyz;
            t_3x3x3x3(1,1,3,3) = xxzz;
            
            % 9+9
            t_3x3x3x3(1,2,1,1) = xxxy;
            t_3x3x3x3(1,2,1,2) = xyxy;
            t_3x3x3x3(1,2,1,3) = xyxz;
            
            t_3x3x3x3(1,2,2,1) = xyxy;
            t_3x3x3x3(1,2,2,2) = yyxy;
            t_3x3x3x3(1,2,2,3) = xyyz;
            
            t_3x3x3x3(1,2,3,1) = xyxz;
            t_3x3x3x3(1,2,3,2) = xyyz;
            t_3x3x3x3(1,2,3,3) = zzxy;
            
            % 27
            t_3x3x3x3(1,3,1,1) = xxxz;
            t_3x3x3x3(1,3,1,2) = xyxz;
            t_3x3x3x3(1,3,1,3) = xzxz;
            
            t_3x3x3x3(1,3,2,1) = xyxz;
            t_3x3x3x3(1,3,2,2) = yyxz;
            t_3x3x3x3(1,3,2,3) = xzyz;
            
            t_3x3x3x3(1,3,3,1) = xzxz;
            t_3x3x3x3(1,3,3,2) = xzyz;
            t_3x3x3x3(1,3,3,3) = zzxz;
            
            % 36
            t_3x3x3x3(2,1,1,1) = xxxy;
            t_3x3x3x3(2,1,1,2) = xyxy;
            t_3x3x3x3(2,1,1,3) = xyxz;
            
            t_3x3x3x3(2,1,2,1) = xyxy;
            t_3x3x3x3(2,1,2,2) = yyxy;
            t_3x3x3x3(2,1,2,3) = xyyz;
            
            t_3x3x3x3(2,1,3,1) = xyxz;
            t_3x3x3x3(2,1,3,2) = xyyz;
            t_3x3x3x3(2,1,3,3) = zzxy;
            
            % 45
            t_3x3x3x3(2,2,1,1) = xxyy;
            t_3x3x3x3(2,2,1,2) = yyxy;
            t_3x3x3x3(2,2,1,3) = yyxz;
            
            t_3x3x3x3(2,2,2,1) = yyxy;
            t_3x3x3x3(2,2,2,2) = yyyy;
            t_3x3x3x3(2,2,2,3) = yyyz;
            
            t_3x3x3x3(2,2,3,1) = yyxz;
            t_3x3x3x3(2,2,3,2) = yyyz;
            t_3x3x3x3(2,2,3,3) = yyzz;
            
            % 54
            t_3x3x3x3(2,3,1,1) = xxyz;
            t_3x3x3x3(2,3,1,2) = xyyz;
            t_3x3x3x3(2,3,1,3) = xzyz;
            
            t_3x3x3x3(2,3,2,1) = xyyz;
            t_3x3x3x3(2,3,2,2) = yyyz;
            t_3x3x3x3(2,3,2,3) = yzyz;
            
            t_3x3x3x3(2,3,3,1) = xzyz;
            t_3x3x3x3(2,3,3,2) = yzyz;
            t_3x3x3x3(2,3,3,3) = zzyz;
            
            % 63
            t_3x3x3x3(3,1,1,1) = xxxz;
            t_3x3x3x3(3,1,1,2) = xyxz;
            t_3x3x3x3(3,1,1,3) = xzxz;
            
            t_3x3x3x3(3,1,2,1) = xyxz;
            t_3x3x3x3(3,1,2,2) = yyxz;
            t_3x3x3x3(3,1,2,3) = xzyz;
            
            t_3x3x3x3(3,1,3,1) = xzxz;
            t_3x3x3x3(3,1,3,2) = xzyz;
            t_3x3x3x3(3,1,3,3) = zzxz;
            
            % 72
            t_3x3x3x3(3,2,1,1) = xxyz;
            t_3x3x3x3(3,2,1,2) = xyyz;
            t_3x3x3x3(3,2,1,3) = xzyz;
            
            t_3x3x3x3(3,2,2,1) = xyyz;
            t_3x3x3x3(3,2,2,2) = yyyz;
            t_3x3x3x3(3,2,2,3) = yzyz;
            
            t_3x3x3x3(3,2,3,1) = xzyz;
            t_3x3x3x3(3,2,3,2) = yzyz;
            t_3x3x3x3(3,2,3,3) = zzyz;
            
            % 81
            t_3x3x3x3(3,3,1,1) = xxzz;
            t_3x3x3x3(3,3,1,2) = zzxy;
            t_3x3x3x3(3,3,1,3) = zzxz;
            
            t_3x3x3x3(3,3,2,1) = zzxy;
            t_3x3x3x3(3,3,2,2) = yyzz;
            t_3x3x3x3(3,3,2,3) = zzyz;
            
            t_3x3x3x3(3,3,3,1) = zzxz;
            t_3x3x3x3(3,3,3,2) = zzyz;
            t_3x3x3x3(3,3,3,3) = zzzz;
            
        end
        
        
        function x = inner(t1,t2)
            x = t1 * t2';
        end
        
        function x = outer(t1,t2)
            x = t1' * t2;
        end
        
        function L = eigval(t_1x6)
            assert(size(t_1x6,2) == 6)
            L = zeros(size(t_1x6,1), 3);
            for i = 1:size(t_1x6,1)
                t_3x3 = Tensor.t_1x6_to_3x3(t_1x6(i,:));
                L(i,:) = eigs(t_3x3);
            end
        end
        
        
        function md = md(t_1x6)
            assert(size(t_1x6,2) == 6)
            md = Tensor.inner(t_1x6, Tensor.t_3x3_to_1x6(Tensor.iso_3x3()));
        end
        
        function v_lambda = v_lambda(t_1x6)
            assert(size(t_1x6,2) == 6)
            [~,E_shear] = Tensor.iso_6x6();
            v_lambda = Tensor.inner(Tensor.t_1x6_to_1x21(t_1x6), Tensor.t_6x6_to_1x21(E_shear));
        end
        
        function fa = fa(a,b)
            if nargin > 1
                md = a;
                vl = b;
            else
                assert(size(a,2) == 6)
                md = Tensor.md(a);
                vl = Tensor.v_lambda(a);
            end
            fa = sqrt((3/2) * (1 + md.^2./vl).^(-1));
            fa = real(fa);
        end
    end
end
