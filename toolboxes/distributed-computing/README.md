# Distributed Computing in Matlab

This project contains a set of functions to help distribute the execution of Matlab functions on a cluster.

## Prerequisites

- The cluster (i.e. "the server") and the main workstation (i.e., "the client") should to be able to access some shared disk space. It is usually possible to do so by mounting a network folder on both stations. If the server and client do not share the same network, mounting can be performed over ssh with sshfs ([Linux](https://doc.ubuntu-fr.org/sshfs), [MacOS](https://osxfuse.github.io), [Windows](https://github.com/Foreveryone-cz/win-sshfs))

- The client should have an SSH software installed. Both [ssh](https://doc.ubuntu-fr.org/ssh), for linux, and [PuTTY](http://www.putty.org), for Windows, are currently handled.

- The client should be able to connect to the cluster without having to type a password. This is usually managed by registering RSA keys on the cluster. On linux, this is done by

    1) Generate a set of public and private keys
    ```shell
    ssh-keygen -t rsa
    ```

    2) Register the key with the cluster
    ```shell
    ssh-copy-id login@cluster
    ```

- Matlab should be installed on the server, with an infinite (or large) number of licenses.

### Limitations

- For now, we only support clusters managed with Sun Grid Engine. It should also work with its forks ([OGS](http://gridscheduler.sourceforge.net), ...) but they were not tested. We plan on supporting other queuing systems, such as PBS, in the future?

## Usage

### Configuration

First an option structure must be generated with the `distribute_default` function. Beforehand, a number of mandatory options should be manually set (else, jobs wil be run locally).

Here is an example of a typical configuration:
```matlab
opt = struct;

opt.server.ip      = 'cluster.university.ac.uk';
opt.server.login   = 'me';
opt.server.folder  = '/home/me/distribute';
opt.client.folder  = '/Users/me/distribute';

opt.matlab.bin     = '/share/apps/matlab';

opt.translate      = {'/Users/me/mydata' '/home/me/data'};

opt = distribute_default(opt);
```

### Options

Several additional options can be set in order to specify the cluster configuration more precisely, or to load matlab packages at runtime. Here is the complete list:

### GENERAL
```
mode          - Parallelisation mode: 'qsub'/'parfor'/'for'                 - ['for']
verbose       - Speak during processing                                     - [false]
```

#### CLUSTER
```
server.ip     - IP adress (or alias name) of the cluster                    - ['' = (server == client)]
server.login  - Login with which to connect                                 - ['']
server.source - Files to source on server side                              - [try to find bashrc and/or bash_profile]
server.folder - Shared folder for writing data, scripts, etc.               - ['~/.distribute']
```

#### LOCAL
```
client.source  - Files to source on server side                             - [auto]
client.workers - Number of local workers                                    - [auto]
client.folder  - Shared folder for writing data, scripts, etc.              - ['~/.distribute']
```

#### SUBMIT JOBS
```
ssh.type      - SSH software to use 'ssh'/'putty'                           - [try to detect]
ssh.bin       - Path to the ssh binary                                      - [try to detect]
ssh.opt       - SSH options                                                 - ['-x']
sched.sub     - Path to the submit binary                                   - [try to detect]
sched.stat    - Path to the stat binary                                     - [try to detect]
sched.acct    - Path to the acct binary                                     - [try to detect]
sched.type    - Type of scheduler 'sge'/'pbs'                               - [try to detect]
job.batch     - Submit jobs as a batch (force same mem for all)             - [true]
job.mem       - (Initial) Max memory usage by a single job                  - ['2G']
job.est_mem   - Estimate max memory usage from previous runs                - [true]
job.sd        - Amount of extra memory to add to estimated max memory       - [0.1]
job.use_dummy - Uses a dummy job to decide when job have finished           - [false]
```

#### MATLAB
```
matlab.bin    - Path to matlab binary                                       - [try to detect]
matlab.add    - Paths to add to Matlab path                                 - [{}]
matlab.addsub - Paths to add to Matlab path, with subdirectories            - [{}]
matlab.opt    - Commandline options to pass to matlab                       - ['-nojvm -nodesktop -nosplash -singleCompThread']
spm.path      - Path to SPM                                                 - []
spm.toolboxes - List of SPM toolboxes to add to Matlab path                 - [{}]
```

#### DATA
```
translate  - Cell array of size 2xN with translation between client and     - [{client.folder server.folder}].
             server paths. Example:
                  {'/home/me/'     '/mnt/users/me/' ;
                   '/shared/data/' '/mnt/shared/data'}
restrict   - Restrict translation to a class: 'char'/'file_array'           - ['']
clean      - Clean tmp data when finished                                   - [true]
clean_init - Initially clean tmp data                                       - [false]
```

### Run

The main function is `distribute`. Its syntax, is quite straightforward: it takes the option structure, a function name or handle, and the list of arguments to pass to the function. Arguments that should be sliced (i.e., iterated over), should be preceded by `'iter'`. Arguments that should be sliced *and* are both inputs and outputs (in particular, structure arrays) should be preceded by `'inplace'`. The option structure is returned because some functionalities (RAM usage estimation) need this structure to be updated.

```
FORMAT [opt, out1, ...] = distribute(opt, func, ('iter'/'inplace'), arg1, ...)

opt  - Option structure. See 'help distribute_default'.
func - Matlab function to apply (string or function handle)
arg  - Arguments of the function
       > Arguments that should be iterated over should be preceded
         with 'iter'
       > Arguments that should be iterated in place (i.e. the output
         replaces the input) should be preceded with 'inplace'
out  - Output of the function. Each one is a cell array of outputs,
       unless some arguments were 'inplace'. In this case, the
       function should return inplace arguments first, in the same
       order as their input order.
```

#### Compact versions

To parallelise locally, one can also replace the option structure with the number of workers. In this case, no option structure is returned:
```
FORMAT [out1, ...] = distribute(nworkers, func, ('iter'/'inplace'), arg1, ...)
```

To loop without parallelisation, one can also remove entirely the option argument. In this case, no option structure is returned:
```
FORMAT [out1, ...] = distribute(func, ('iter'/'inplace'), arg1, ...)
```

### Examples

Let us give a first use case, were we have a list of pairs of arrays that should be summed:
```
% Initialise arrays
N   = 10;
DIM = [5 5];
a   = cell(1,N);
b   = cell(1,N);
for i=1:N
    a{i} = randn(DIM);
    b{i} = randn(DIM);
end

% Local processing
true_c = cellfun(@plus, a, b, 'UniformOutput', false);

% Distributed processing
[opt, dist_c] = distribute(opt, 'plus', 'iter', a, 'iter', b);
```

Let us now apply a distributed process to a structure array, which is both input and output
```
% Initialise structure array
N = 10;
f = cell(1,N);
a = struct('f', f);

% Set the value 3 in all fields 'f'
[opt, a] = distribute(dist, @setfield, 'inplace', a, 'f', 3);
```

## Future developments

We intend to allow:
- job batching, where a single job processes several "subjects".
- optimising cluster use by choosing between local and distributed processing based on the cluster load.
- distributing scripts/binaries on top of Matlab functions. This can be helpful for working with compiled Matlab scripts.

## Contributors

This software was developed by Mikael Brudfors and Yaël Balbastre in John Ashburner's [Computational Anatomy Group](http://www.fil.ion.ucl.ac.uk/Ashburner/) at the [Wellcome Centre for Human Neuroimaging](http://www.fil.ion.ucl.ac.uk/) in UCL.

If you encounter any difficulty, please send an email to `y.balbastre` or `mikael.brudfors.15` _at_ `ucl.ac.uk`

## License

This software is released under the [GNU General Public License version 3](LICENSE) (GPL v3). As a result, you may copy, distribute and modify the software as long as you track changes/dates in source files. Any modifications to or software including (via compiler) GPL-licensed code must also be made available under the GPL along with build & install instructions.


[TL;DR: GPL v3](https://tldrlegal.com/license/gnu-general-public-license-v3-(gpl-3))
