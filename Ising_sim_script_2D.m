%% simulation_bluehive

pc = parcluster('local')
JOB_ID = getenv('SLURM_JOBID')
CPUS = getenv('SLURM_CPUS_PER_TASK')
pc.JobStorageLocation = strcat('/local_scratch/',JOB_ID)
parpool(pc,str2num(CPUS))

addpath('classes');
%addpath('.');
Ls_ = [16,32,64,128];
numruns = 100;

k = 1;
for i = str2num(getenv('SLURM_ARRAY_TASK_ID'))
	model = Ising_2D_sim(Ls_,i);
    model.runsim_hybrid();
	k = k + 1;
end
delete(gcp('nocreate'))
