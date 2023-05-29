% Matlab matrix beam solver (2D)
clear;
clc;
% force/displacment matrixes (also imput into solve function below)
% declered as symbilic variables
syms 'V%d' 'M%d' 'v%d' 'a%d' [1 10]

% material properties
E = 200e9;                              % Pa
I = 100e-6;                  % m^4 

% input force and displacement matrixes as column matrixes
Q = [ V1;  M1; 0; 0];      % N, N/m
Q_dash = [  -28; 9.33e3; 34.96e3; 58.6e3];      % N, N/m
D = [ 0;0;v3;a4];            % m, rad

% beam lengths (nombers correspond to beam lables)
L = [ 2, 5, 6]; % m

% axis labes at each beam (grouped by beam nomber, order of pars subject to beam direction)
NodeAxisLbl = {[1,2,3,4];
               [3,4,5,6];
               [5,6,7,8]};


k = cell(1, length(L));
K_size = max(size(Q));
K = zeros(K_size);
for Beam = 1:length(L)
    % building individual stiffnes matrixes
    k{Beam} = Mstiff(E,I,L(Beam));
    
    % Merging stiffnes matrixes into large 'K' matrix
    K = K + k_resort(k{Beam}, NodeAxisLbl{Beam}, K_size);
end

% solve the matrix equation varialbles (gives node displacment and cotraint reactions)
AnsInStruct = solve( Q==K*D - Q_dash);

% Print out answers
Ans_names = fieldnames(AnsInStruct);
AnsInCell = struct2cell(AnsInStruct);
fprintf('Force in (N), Moments in (N/m), angles in (rad)\n');
for i=1:length(Ans_names)
    fprintf('    %s:  %0.4d\n',Ans_names{i}, double(AnsInCell{i}));%#ok<SAGROW>
end

% % clean up workspace
clearvars -except Q D k K E I L axisLbl NodeAxisLbl

% solving manualy:
% S1 = cat(2,double(K([3,7,8],[1,2,4,5,6,9,10])), Q([3,7,8]))
% s1=double(S1);
% Gau(s1);

function M = Gau(m)
M=m;
rows = size(m,1);
for d=1:rows
    M(d,:) = M(d,:)/M(d,d);
    for i=d+1:rows
        M(i,:) = M(i,:) - M(i,d)*M(d,:);
    end
end
end

function k = Mstiff(E, I, L)
k = (I*E/L^3)*  [12,   6*L,  -12,    6*L;
                6*L, 4*L^2, -6*L, -2*L^2;
                -12,  -6*L,   12,   -6*L;
                6*L, 2*L^2, -6*L,  4*L^2];
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

