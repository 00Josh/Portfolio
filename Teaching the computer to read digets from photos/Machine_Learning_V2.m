%% Machine Learning V2

%% First go at Machine Learing
% This program is the first attemt at a nural network
% based on the info at http://neuralnetworksanddeeplearning.com/
% currently learns binery, but te layers and sizes of each layer can be
% edited
% this uses for loops and a lot of cell structures wich are very slow

% clear;
clc;
close all;

%% import database

[img, lbl] = import_training_data('train-images.idx3-ubyte', 'train-labels.idx1-ubyte');
% to open: image(images(:,:,x))

% reshape each image to be a column vector X
X = reshape( img, size(img,1)*size(img,2), [], 1)/255;

% % debug version (learn binery 0-3)
% X = [0,1,0,1;0,0,1,1];lbl = [0,1,2,3];
% X = join(X, ex);
% lbl = join(lbl, ex);


%% set up layer structure and batch size and initual vaues
% declar size of each layer in the nural network(can change length and nomber of layers)
layer_str = [ size(X,1), 80, 40, 10];

% log every enth entry
log_per = 8;  %traing images

% declare batch size and laring rate
sz_bch = 14;
learnR = 95;
loops = 1;

% structure the weights andf biasses in the correct dementions with normal
% distribution
[W, b] = init_val(layer_str);


%% Prealocate stuff

% declar lyer info
lys = length(layer_str);

% Prealocate space
a = cell(lys-1, 1);
z = cell(lys-1, 1);
b_dash = cell(sz_bch, size(b,2));
W_dash = cell(sz_bch, size(W,2));
Bcnt = 1;
sigma = cell(lys-1,1);

% set up loging system for progress
thr = log_per;
log_count = 1;
log_W = cell(1, length(lbl)/log_per);
log_b = log_W;
log_Av_C = nan(1, length(lbl)/log_per);
log_Acc = nan(1, length(lbl)/log_per);
Av_acc = 0;
Av_C = 0;


%% start training by baches

% Itterate by ten immages at a time
Img_cnt = 1;
% tic;
for restarts = 1:loops
    while  Img_cnt+sz_bch<size(X,2)%toc<120 &&
        for Im = Img_cnt:Img_cnt+sz_bch-1
            Img_cnt=Img_cnt+1;

            % find the Z
            a{1} = double(X(:,Im));
            for i = 1:lys - 1
                z{i} = W{i}*a{i} + b{i};
                a{i+1} = sig(z{i});
            end

            % find the sigma function
            sigma{end} = (a{end} - logic(lbl(Im),size(a{end},1))).*sig_dash(z{end});
            for l = lys-2:-1:1
                sigma{l} = (W{l+1}'*sigma{l+1}).*sig_dash(z{l});
            end

            % compute sigma
            for L = 1:lys-1
                % updata wiegths and bises
                b_dash{Bcnt, L} = b{L} - sigma{L}*learnR/sz_bch;
                W_dash{Bcnt, L} = W{L} - sigma{L}*a{L}'*learnR/sz_bch;
            end
            Bcnt = Bcnt + 1;


            % record mooving average accuracy and cost
            out = find(max(a{end})==a{end})-1;
            if out==lbl(Im)
                crct = 1;
            else
                crct = 0;
            end
            Av_acc = mvavg(Av_acc, crct, 200);
            Av_C = mvavg(Av_C, cost(lbl(Im), a{end}), 200);

            % output info
    %         fprintf('Individual Cost: %0.3f\n', (lbl(Im) - a{end}).^2);
    %         fprintf('Cost: %0.3f\n', cost(lbl(Im), a{end}));         
    %         fprintf('Ans: %d, Image: %d, Correct: %d, Acc %0.4f\n\n', find(max(a{end})==a{end}), lbl(Im), crct, Av_acc);

            % output the image that is currently in use
    %         Q = figure(1);
    %         image(img(:,:,Im));
    %         Q.Position = [-998.2000 259.4000 560 420.0000];
    %         pause(0.01);

            % log the wieghts and biases every 2000 trainig images
            if Img_cnt > thr
                log_W{log_count} = W;
                log_b{log_count} = b;
                log_Av_C(log_count) = Av_C;
                log_Acc(log_count) = Av_acc;
                thr = thr + log_per;
                log_count = log_count + 1;
            end
        end

        % take the average of the b_dash and W_dash and apply
        Bcnt = 1;
        b = avg(b_dash);
        W = avg(W_dash);
    end
    
end

% output info
% plot progres
f1 = figure;
plot(log_Av_C);
% f1.Position =  [-1483, 273, 560, 420];
title('Recorded Average Cost');
f2 = figure;
plot(log_Acc);
% f2.Position =  [-876.6, 279.4, 560, 420];
title('Recorded Average Accuracy');
fprintf('Layer Structure\n');
fprintf('%d\n', layer_str);
fprintf('learning rate %d\nBatch Size %d\nLoops %d\nLoging every %d\n\n', learnR, sz_bch, loops, log_per);


%% Functions

function y = join(x, n)
for i=1:n
    x = cat(2, x, x);
end
y = x;
end

function av = mvavg(cur, new, lgt)
av = ((lgt-1)*cur+new)/lgt;
end

function B = avg(A)
% take average verticaly along each column (cell imput)
B = cell(1, size(A,2));
for j = 1:size(A,2)
    a = zeros(size(A{1,j}));
    for i = 1:size(A,1)
        a = a + A{i,j};
    end
    B{j} = a/i;
end
end

function [W, b] = init_val(sz)
layN = length(sz)-1;
W = cell(1,layN);
b = cell(1,layN);
for i = 1:layN
    W{i} = randn(sz(i+1),sz(i));
    b{i} = randn(sz(i+1),1);
end
end

function out = sig(x)
out = (1+exp(-x)).^-1;
end

function out = sig_dash(x)
out = exp(-x).*(1+exp(-x)).^-2;
end

function C = cost(An, a)
y = logic(An, size(a,1));
C = sum((y - a).^2);
end

function A = logic(y, sz)
A = zeros(sz,1);
A(y+1) = 1;
end

function [images, labels] = import_training_data(path_to_digits, path_to_labels)

% Open files
fid1 = fopen(path_to_digits, 'r');

% The labels file
fid2 = fopen(path_to_labels, 'r');

% Read in magic numbers for both files
A = fread(fid1, 1, 'uint32');
magicNumber1 = swapbytes(uint32(A)); % Should be 2051
fprintf('Magic Number - Images: %d\n', magicNumber1);

A = fread(fid2, 1, 'uint32');
magicNumber2 = swapbytes(uint32(A)); % Should be 2049
fprintf('Magic Number - Labels: %d\n', magicNumber2);

% Read in total number of images
% Ensure that this number matches with the labels file
A = fread(fid1, 1, 'uint32');
totalImages = swapbytes(uint32(A));
A = fread(fid2, 1, 'uint32');
if totalImages ~= swapbytes(uint32(A))
    error('Total number of images read from images and labels files are not the same');
end
fprintf('Total number of images: %d\n', totalImages);

% Read in number of rows
A = fread(fid1, 1, 'uint32');
numRows = swapbytes(uint32(A));

% Read in number of columns
A = fread(fid1, 1, 'uint32');
numCols = swapbytes(uint32(A));

fprintf('Dimensions of each digit: %d x %d\n', numRows, numCols);

% For each image, store into an individual slice
images = zeros(numRows, numCols, totalImages, 'uint8');
for k = 1 : totalImages
    % Read in numRows*numCols pixels at a time
    A = fread(fid1, numRows*numCols, 'uint8');

    % Reshape so that it becomes a matrix
    % We are actually reading this in column major format
    % so we need to transpose this at the end
    images(:,:,k) = reshape(uint8(A), numCols, numRows).';
end

% Read in the labels
labels = fread(fid2, totalImages, 'uint8');

% Close the files
fclose(fid1);
fclose(fid2);

end

