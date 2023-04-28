classdef RigidBodySolver < handle
    
    properties
       Mass, Inertia, invInertia, dt, g , currentWr2omegaDot
    end
    
    methods
        function obj = RigidBodySolver(Mass, Inertia, dt,g)
            obj.Mass = Mass; 
            obj.Inertia = Inertia;
            obj.invInertia = inv(Inertia);
            obj.dt = dt;
            obj.g = g;
        end
        
        function state = nextStep(RBS, currentState, Force, Moments)
            K = zeros(12, 4);
    
            K(:, 1) = RBS.dt*DOF6(RBS, currentState ,Force, Moments);
            K(:, 2) = RBS.dt*DOF6(RBS, currentState+0.5*K(:, 1) ,Force, Moments);
            K(:, 3) = RBS.dt*DOF6(RBS, currentState+0.5*K(:, 2) ,Force, Moments);
            K(:, 4) = RBS.dt*DOF6(RBS, currentState+K(:, 3) ,Force, Moments);

            state = currentState + (...
                K(:, 1)+...
                2*K(:, 2)+...
                2*K(:, 3)+...
                K(:, 4))/6;
        end
        
        function F = DOF6(RBS, currentState, forces, Moments)
            
            % (Sin, Cos, Tan) of (phi, theta, epsi)
            [S, C, T] = SCT(RBS, currentState(7:9));
            s_theta = S.theta;
            c_theta = C.theta;
            t_theta = T.theta;
            s_epsi = S.epsi;
            c_epsi = C.epsi;
            s_phi = S.phi;
            c_phi = C.phi;

            Forces = forces + inertia2Body(RBS, S, C)*[
                0;
                0;
                RBS.Mass*RBS.g
            ];

            % (u, v, w) dot
            u_v_w_dot = (1/RBS.Mass)*Forces - cross(...
                currentState(4:6, 1), currentState(1:3, 1)...
            );

            % (p, q, r) dot
            p_q_r_dot = RBS.invInertia *(Moments - cross(...
                currentState(4:6, 1), RBS.Inertia * currentState(4:6, 1)...
            ));

            % (phi, theta, epsi) dot
            RBS.currentWr2omegaDot = [
                1, s_phi*t_theta, c_phi*t_theta;
                0, c_phi, -s_phi;
                0, s_phi/c_theta, c_phi/c_theta;
            ];
            phi_theta_epsi_dot = RBS.currentWr2omegaDot * currentState(4:6, 1);

            % (x, y, z) dot
            x_y_z_dot = [
                c_theta*c_epsi, (s_phi*s_theta*c_epsi - c_phi*s_epsi), (c_phi*s_theta*c_epsi + s_phi*s_epsi);
                c_theta*s_epsi, (s_phi*s_theta*s_epsi + c_phi*c_epsi), (c_phi*s_theta*s_epsi - s_phi*c_epsi);
                -s_theta, s_phi*c_theta, c_phi*c_theta
            ] * currentState(1:3, 1);

            F = [u_v_w_dot; p_q_r_dot; phi_theta_epsi_dot; x_y_z_dot];
            
        end
        
        function [S, C, T] = SCT(RBD, s)
            
            sct = [...
                sin(s(1)), sin(s(2)), sin(s(3));
                cos(s(1)), cos(s(2)), cos(s(3));
                tan(s(1)), tan(s(2)), tan(s(3));
            ];
            S = struct(...
                'phi', sct(1,1),...
                'theta', sct(1,2),...
                'epsi', sct(1,3)...
            );
            C = struct(...
                'phi', sct(2,1),...
                'theta', sct(2,2),...
                'epsi', sct(2,3)...
            );
            T = struct(...
                'phi', sct(3,1),...
                'theta', sct(3,2),...
                'epsi', sct(3,3)...
            );
            
        end
        
        function R_i2b = inertia2Body(~, S, C)
            R1 = [...
                1       0     0;
                0   C.phi S.phi;
                0  -S.phi C.phi;
            ];
            R2 = [...
                C.theta   0 -S.theta;
                      0   1        0;
                S.theta   0  C.theta;
            ];
            R3 = [...
                C.epsi S.epsi 0;
               -S.epsi C.epsi 0;
                     0      0 1;
            ];
            R_i2b = R1*R2*R3;
            
        end
        
    end
end

