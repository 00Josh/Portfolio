% Matlab matrix truss solver (2D)
clear;
clc;

% To use: make sure that the beam legnths(L), NodeAxisLbl (..AxisLabels)
%         and axis/beam angles are all in the same beam-order

% material properties
E = 200e9;                              % Pa
A = (0.100^2-0.090^2);                  % m^2 

% force/displacment matrixes (also imput into solve function below)
syms 'R%d' 'x%d' 'y%d' [1 10]% declered as symbilic variables
Q = [ 9e3; -12e3; R3; 0; 0; -15e3; R7; R8; 0; -10e3];      % N
D = [ x1; y2; 0; y4; x5; y6; 0; 0; x9; y10];             % m

% beam lengths (nombers correspond to beam lables)
L = [ 3, 4, 3, 4, 5, 5, 3, 3*(2)^.5]; % m

% axis labes at each beam (grouped by beam nomber, order of pars subject to beam direction)
NodeAxisLbl = {[7,8,1,2];
               [1,2,3,4];
               [3,4,5,6];
               [5,6,7,8];
               [3,4,7,8];
               [5,6,1,2];
               [7,8,9,10];
               [9,10,1,2]};

% axis angles (from gloabl 'x' axis), to mach 'axisLbl', each row is a beam
% and each column is each end of beam
axisAgl = [0, 0;  % degrees
           0,310;
         310, 0;
           0, 0;
         310, 0;
           0, 0;
           0, 0;
           0, 0;];

% beam angles from global x axis, direction should piont in the same
% directin as the arrys in 'axisLbl'
beamAgl = [ 90; % degrees
           180;
           270;
             0;
360-atand(3/4);
    atand(3/4);
             0;
         90+45];




% building individual stiffnes matrixes
k = cell(1, length(L));
for Beam = 1:length(L)
    k{Beam} = stiff(beamAgl(Beam), axisAgl(Beam,:),E,A,L(Beam));
    %k{Beam} = stiff(lamnax(theta(Beam)),lamnay(theta(Beam)),E,A,L(Beam));
end

% Merging stiffnes matrixes into large 'K' matrix
K_size = max(size(Q));
K = zeros(K_size);
for Beam = 1:length(L)
    K = K + k_resort(k{Beam}, NodeAxisLbl{Beam}, K_size);
end

% solve the matrix equation varialbles (gives node displacment and cotraint reactions)
AnsInStruct = solve(Q==K*D);

% % solving manualy:
% Q_relavant_rows = [1,2,4,5,6,9,10];
% D_relavant_rows = [1,2,4,5,6,9,10];
% [Q_ans, D_ans] = Gau(Q,K,D, Q_relavant_rows, D_relavant_rows)

% reformat output
AnsInCell = struct2cell(AnsInStruct);
for i=1:max(size(AnsInCell))
    AnsDecimal(i) = double(AnsInCell{i}); %#ok<SAGROW>
end

% find beam axial forces
D_des = symToDecimal(D,AnsInStruct, AnsDecimal);
F_beam = cell(length(L),1);
F_beamAx = nan(length(L),1);
for beam=1:length(L)
    F_beam{beam} = k{beam}*D_des(NodeAxisLbl{beam});
    F_beamAx(beam) = (F_beam{beam}(1)^2 + F_beam{beam}(2)^2)^0.5;
        
    % check if ends of same beam have the same axial force (+-)1N
    if abs(F_beamAx(beam) - (F_beam{beam}(3)^2 + F_beam{beam}(4)^2)^0.5) > 1
        fprintf('Error: beam angle dose not match axis and beam axial force direction\n')
    end
    
    % chech whether in compression or tension (by cheching if rsultant pionts towards or away from beam)
    a = atan2d(F_beam{beam}(2),F_beam{beam}(1)) + axisAgl(beam,1) - beamAgl(beam);
    A_dif = abs(mod360(a));
    if A_dif < 10
        F_beamAx(beam) = -1*F_beamAx(beam);
    elseif abs(A_dif-180) < 10
    else
        error('axial force resultant angle does not piont along beam');
    end
end

% Print out answers
Ans_names = fieldnames(AnsInStruct);
AnsInCell = struct2cell(AnsInStruct);
fprintf('\n Force in (N), Displacments (m)\n');
for i=1:length(Ans_names)
    fprintf('    %s:  %0.4d\n',Ans_names{i}, double(AnsInCell{i}));%#ok<SAGROW>
end
fprintf('\nBeam Axil Forces (N)\n');
for i=1:length(F_beamAx)
    fprintf('    Beam %0.0f: %0.5d\n', i, F_beamAx(i));
end

% clean up workspace
clearvars -except AnsDecimal Q D k K E A L beamAgl axisAgl axisLbl NodeAxisLbl F_beamAx





function k = stiff(beamAngle, axisAngle, E, A, L)

lamx = [ cosd(beamAngle - axisAngle(1)),  cosd(beamAngle - axisAngle(2))];
lamy = [ cosd(90 - beamAngle + axisAngle(1)),  cosd(90 - beamAngle + axisAngle(2))];

a = [      lamx(1)^2,  lamx(1)*lamy(1),    -lamx(1)*lamx(2),    -lamx(1)*lamy(2);
     lamx(1)*lamy(1),        lamy(1)^2,    -lamx(2)*lamy(1),    -lamy(1)*lamy(2);
    -lamx(1)*lamx(2), -lamx(2)*lamy(1),           lamx(2)^2,     lamx(2)*lamy(2);
    -lamx(1)*lamy(2), -lamy(1)*lamy(2),     lamx(2)*lamy(2),           lamy(2)^2  ];

k = (A*E/L)*a;
end

function K = k_resort(k, axisLbl, size)
K = zeros(size);
for i=1:size
    for j=1:size
        if RowColMatch(axisLbl,i,j)
            K(i,j) = K(i,j) + k(find(axisLbl==i), find(axisLbl==j)); %#ok<FNDSB>
        end
    end
end
end

function Return = RowColMatch(nodeAxisLbls,row,col)
RC = [row,col];
Return = 1;
for demention=1:2
    if Return==0
        break;
    end
    for l=1:length(nodeAxisLbls)
        if RC(demention)==nodeAxisLbls(l)
            break;
        end
        if l==length(nodeAxisLbls)
            Return=0;
        end
    end
end
end

function a = mod360(A)
while A>360 || A<0
    A = A-(A/abs(A))*360;
end
a=A;
end

function D_des = symToDecimal(D,AnsInStruct, AnsDecimal)
AnsLbles = fieldnames(AnsInStruct);
D_des = nan(length(D), 1);
for i=1:length(D)
    for j=1:length(AnsLbles)
        if symvar(D(i))==AnsLbles(j)
            D_des(i) = AnsDecimal(j);
            break
        elseif j==length(AnsLbles)
            D_des(i) = double(D(i));
        end
    end
end
end

% Gousinan solving (so that you can see what is goning on)
% - when solving manualy this is where you exclude rows and colums in the Q=K*D equation 
% E.G: [Q_ans, D_ans] = Gau(D,K,Q, [3,7,8],[1,2,4] )
function [Q_out, D_out] = Gau(Q,K,D, Q_relevant_rows, D_relevant_rows)
Q_out = nan(size(Q));
D_out = nan(size(D));

% eliminate irrelevant rows and culomns
D = D(D_relevant_rows);
Q = Q(Q_relevant_rows);
K_r = K(D_relevant_rows, Q_relevant_rows);

% test if all dementions match (all sould have the same rows, as well as be detruminant, K(Rows==Columns))
if size(Q,1)~=size(D,1) || size(D,1)~=size(K_r,1) || size(K_r,1)~=size(K_r,2)
    fprintf("Error: can't solve Gausian matrix as system is indterminant\n    Try ajusting the input\n\n");
else
    
% join   
M = cat(2,double(K_r), double(Q));
rows = size(M,1);
% solve with gausian elimination
for d=1:rows
    M(d,:) = M(d,:)/M(d,d);
    for i=d+1:rows
        M(i,:) = M(i,:) - M(i,d)*M(d,:);
    end
end

% sort answers
D_ans = nan(size(Q));
D_ans(end)=M(end,end);
for r=rows-1:-1:1
    D_ans(r) = M(r,end) - M(r, r+1:end-1) * D_ans(r+1:end);
end

Q_ans_rows = rmmissing(standardizeMissing(1:length(Q_out), Q_relevant_rows));

Q_ans = K( Q_ans_rows, D_relevant_rows) * D_ans;

Q_out(Q_ans_rows) = Q_ans;
D_out(D_relevant_rows) = D_ans;
end
end

% plain stiffness matrix not including loacl axis's
%function k = stiff(lamx, lamy, E, A, L)
% 
% a = [   lamx^2, lamx*lamy,    -lamx^2, -lamx*lamy;
%      lamx*lamy,    lamy^2, -lamx*lamy,    -lamy^2;
%        -lamx^2, -lamx*lamy,    lamx^2,  lamx*lamy;
%     -lamx*lamy,    -lamy^2, lamx*lamy,     lamy^2  ];
% 
% k = (A*E/L)*a;
% end

