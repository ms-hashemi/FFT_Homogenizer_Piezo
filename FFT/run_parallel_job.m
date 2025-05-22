% run_parallel_job
%   Run a parallel Matlab job on a 1 cluster node with a pool of 8 Matlab workers.

% First initialize a cluster object based on an existing cluster profile,
% in this case we are using the 'condo R2019b' profile.
c = parcluster('nova R2020b')

% Use the AdditionalProperties property of the cluster object to set job specific details:
c.AdditionalProperties.NumNodes = 1;                           % Number of nodes requested. 
c.AdditionalProperties.EmailAddress = 'mhashemi@iastate.edu';  % Your Email address (please modify).
c.AdditionalProperties.ProcsPerNode = 10;                       % 1 more than number of Matlab workers.
c.AdditionalProperties.WallTime = '2:00:00';                   % The max wall time for the job.
c.AdditionalProperties.QueueName = 'short_1node192';
c.AdditionalProperties.AdditionalSubmitArgs = '--mem=0 --chdir=/work/sheidaei/mhashemi';

% Examples of other properties that you might need to set:
%    To set a specific queue name, in this case to use the 'freecompute' free tier:
%        c.AdditionalProperties.QueueName = 'freecompute';
% 
%    To set the Slurm job name.  (if not set, Matlab will use "JobN" where N is determined by Matlab):
%        c.AdditionalProperties.AdditionalSubmitArgs = '--job-name=xxx';
%    NOTE: The value of AdditionalProperties.AdditionalSubmitArgs is simply added on to the sbatch command
%          so this can be used supply any additional options to sbatch. 
% c.AdditionalProperties.AdditionalSubmitArgs = '-p=amd';
% c.AdditionalProperties.AdditionalSubmitArgs = '--mem=500000';
% Start a job timer for recording the elapsed time for the job:
tic

% The batch command below creates a job object called 'myjob' that runs a
% Matlab job with 8 parallel pool workers.
% NOTE: Matlab will add an additional worker to the pool for its own use so be 
% sure that the number of processors requested from Slurm (NumNodes X ProcsPerNode)
% is greater than the total number of workers needed by Matlab.
% We also set the parameter AutoAddClientPath to false so that Matlab won't complain when paths on 
% your desktop don't exist on the compute node (this is typical and can be ignored).

% myjob = batch(c, 'main3D_newton_energy', 'pool', 9, 'AutoAddClientPath', false, 'AttachedFiles', "000009.zip")
myjob = batch(c, 'main3D_newton_energy_sym_efficient', 'pool', 9, 'AutoAddClientPath', false, 'CurrentFolder', '/work/sheidaei/mhashemi');

% see https://www.mathworks.com/help/parallel-computing/batch.html for additional tips and examples.

% Wait for the job to finish before continuing. 
% wait(myjob)
% diary(myjob)

% load the 'A' array from the job results. (The values for 'A' are calculated in parallel_mywave.m):

% load(myjob,'A');

%-- plot the results --%
% plot(A);

% print the elapsed time for the job:
toc

