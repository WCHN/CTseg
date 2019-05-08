function TestPushPull

win     = [2 2 2];
lat_ref = [172 221 168];
vs_ref  = [1 1 1];
lat_sub = [28 221 168];
vs_sub  = [6 1 1];

M_sub = [diag(vs_sub) [0;0;0]; 0 0 0 1];
M_ref = [diag(vs_ref) [0;0;0]; 0 0 0 1];
% R = [cos(pi/3) -sin(pi/3) 0; sin(pi/3) cos(pi/3) 0; 0 0 1];
R = eye(3);
R = [R [0;0;0]; 0 0 0 1];

A = M_sub\(R*M_ref);
J = inv(A(1:3,1:3));
J = single(reshape(J, [1 1 1 3 3]));

y = spm_warps('affine', inv(A), lat_sub);

ref = zeros(lat_ref, 'single');
[X,Y,Z] = ndgrid((1:lat_ref(1))-floor(lat_ref(1)/2),...
                 (1:lat_ref(2))-floor(lat_ref(2)/2),...
                 (1:lat_ref(3))-floor(lat_ref(3)/2));
ref(sqrt(X.^2+Y.^2+Z.^2)<50) = ref(sqrt(X.^2+Y.^2+Z.^2)<50) + 3;
sub = zeros(lat_sub, 'single');
sub(14,111,84) = 1;

kernel = blur_fun(lat_ref, eye(3), [6 1 1]);

%%
fprintf('Time spm_diffeo pull/push\n');
tic, ref2sub = spm_diffeo('pull', ref, y); toc
tic, sub2ref = spm_diffeo('push', sub, y, lat_ref); toc

%%
fprintf('Time mtv pull/push (0=dirac)\n');
tic, ref2sub = pushpull('pull', ref, y, J, [0 0 0]); toc
tic, sub2ref = pushpull('push', sub, y, J, [0 0 0], lat_ref); toc

%%
fprintf('Time mtv pull/push (1=gauss)\n');
tic, ref2sub = pushpull('pull', ref, y, J, [1 1 1]); toc
tic, sub2ref = pushpull('push', sub, y, J, [1 1 1], lat_ref); toc

%%
fprintf('Time mtv pull/push (2=rect)\n');
tic, ref2sub = pushpull('pull', ref, y, J, [2 2 2]); toc
tic, sub2ref = pushpull('push', sub, y, J, [2 2 2], lat_ref); toc

%%
fprintf('Time fft pull/push\n');
tic, ref2sub = spm_diffeo('pull',single(real(ifftn(fftn(ref).*kernel))),single(y)); toc
tic, sub2ref = single(real(ifftn(fftn(spm_diffeo('push',sub,single(y),lat_ref)).*kernel))); toc

%%
fprintf('Check adjointness mtv pull/push (1=gauss)\n');
win = [1 1 1];
u = randn(lat_ref, 'single');
v = randn(lat_sub, 'single');
Au = pushpull('pull', u, y, J, win);
Atv = pushpull('push', v, y, J, win, lat_ref);
v_dot_Au = v(:)' * Au(:);
Atv_dot_v = Atv(:)' * u(:);
v_dot_Au - Atv_dot_v

%%
fprintf('Check adjointness mtv pull/push (2=rect)\n');
win = [2 2 2];
u = randn(lat_ref, 'single');
v = randn(lat_sub, 'single');
Au = pushpull('pull', u, y, J, win);
Atv = pushpull('push', v, y, J, win, lat_ref);
v_dot_Au = v(:)' * Au(:);
Atv_dot_v = Atv(:)' * u(:);
v_dot_Au - Atv_dot_v

foo = 0;