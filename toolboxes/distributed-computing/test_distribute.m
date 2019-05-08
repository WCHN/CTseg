function test_distribute

% -------
% OPTIONS
% -------

dist = struct;

dist.server.ip      = 'holly';
dist.server.login   = 'ybalba';
dist.server.folder  = '/data/ybalba/distribute';
dist.client.folder  = '/Users/balbasty/Desktop/FIL/data/distribute';

dist.matlab.bin     = '/share/apps/matlab';
dist.matlab.add     = {'/data/ybalba/matlab/shape-toolbox' ...
                       '/data/ybalba/matlab/utility-functions'};
dist.matlab.addsub  = '/data/ybalba/matlab/viewers-toolbox';
dist.spm.path       = '/data/ybalba/matlab/spm-trunk';
dist.spm.toolboxes  = {'Shoot'};

dist.translate      = {'/Users/balbasty/Desktop/FIL/data' '/data/ybalba'};
dist.restrict       = 'file_array';
dist.clean          = true;

dist.job.batch = true;

dist = distribute_default(dist);

%% ------
%  Test 1
%  ------

% t = {60 60 60 60 60};
t = num2cell(10:10:100);
% t = {1 1 1 1 1};
[dist,foo] = distribute(dist, @pause, 'iter', t);

%% ------
%  Test 2
%  ------

N = 10;
DIM = [5 5];
a = cell(1,N);
b = cell(1,N);
for i=1:N
    a{i} = randn(DIM);
    b{i} = randn(DIM);
end

true_c = cellfun(@plus, a, b, 'UniformOutput', false);
dist_c = distribute(dist, 'plus', 'iter', a, 'iter', b);

%% ------
%  Test 3
%  ------

N = 10;
f = cell(1,N);
a = struct('f', f);
a = distribute(dist, @setfield, 'inplace', a, 'f', 3);


