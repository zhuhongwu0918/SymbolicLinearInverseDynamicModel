function robot = CreateScara()
%% robot = CreateScara()
% Creates a robot model of the Scara that is compatible with the
% spatial_v2 dynamics library.

a1 = .1;
a2 = .2;

aa1 = sym('a1');
aa2 = sym('a2');

% Modified Denavit Hartenberg parameters alpha_{i-1}, a_{i-1}, and d_i 
% as in J.J.Craig's classic text
dhTable = [0 0 0; 
           0 a1 0 ; 
           0 a2 0 ; 
           0 0 0 ;  
           ];
       
       
dhTable_sym = [0 0 0; 
           0 aa1 0 ; 
           0 aa2 0 ; 
           0 0 0 ;  
           ];

% Constraints on motor parameters to ensure rotational symmetry (about z)
motor_symmetry_constraint = zeros(6,10);
motor_symmetry_constraint(1,2) = 1;
motor_symmetry_constraint(2,3) = 1;
motor_symmetry_constraint(3,5) = 1; motor_symmetry_constraint(3,6) = -1;
motor_symmetry_constraint(4,8) = 1;
motor_symmetry_constraint(5,9) = 1;
motor_symmetry_constraint(6,10)= 1; 


B = RangeBasis( null(motor_symmetry_constraint) ); 
P = B*B'; % Projector onto nullspace of symmetry constraint

robot.NB = 4;
robot.parent = [0 1 2 3];
robot.jtype = {'Rz', 'Rz', 'Pz','Rz'};
robot.jtype_motor = {'Rz', 'Rz', 'Rz','Rz'};

robot.gr = {5, 5,5,5};
for i = 1:4
    robot.gr{i} = 5;
end

robot.Xtree = {  };
for i = 1:4
    alpha = dhTable(i,1);
    a = dhTable(i,2);
    d = dhTable(i,3);
    
    a_s = dhTable_sym(i,2);
    robot.Xtree{i} = xlt([0 0 d]) * xlt([a 0 0]) * rotx(alpha);
    %xlt空间坐标平移rotx空间坐标绕旋转(X-axis rotation)
    robot.Xtree_sym{i} = xlt([a_s 0 0]) * round(rotx(alpha));
    robot.Xtree_motor{i} = robot.Xtree{i};
    robot.Xtree_motor_sym{i} = robot.Xtree_sym{i};
    
    robot.motor_constraint{i} = P;
    
end

robot.gravity = [0 0 -9.81];		
robot.I = { eye(6),eye(6),eye(6),eye(6)};
robot.I_motor = { eye(6),eye(6),eye(6),eye(6)};