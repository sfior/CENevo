#!/bin/bash

#SBATCH --job-name=angsd_gl
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4000
#SBATCH --tmp=30000
#SBATCH --time=24:00:00
#SBATCH --output=./logs/angsd_gl.%A.%a.log
#SBATCH --error=./errs/angsd_gl.%A.%a.err

###############################################################
#### Notes ###
# For deciding between -c, --cpus-per-task vs -n, --ntasks; see: https://stackoverflow.com/questions/51139711/hpc-cluster-select-the-number-of-cpus-and-threads-in-slurm-sbatch and https://support.ceci-hpc.be/doc/_contents/SubmittingJobs/SlurmFAQ.html#Q05
# -c, --cpus-per-task vs -n, --ntasks. Note that generally: -n cores, --ntasks=cores should be used for MPI jobs/processes (distributed memory) and --ntasks=1 --cpus-per-task=cores for OpenMP jobs/process (shared memory). See: https://scicomp.ethz.ch/wiki/LSF_to_Slurm_quick_reference#Parallel_job
# The --ntasks value refers to the number of tasks to be launched by SLURM only. This usually equates to the number of MPI tasks launched. Reduce this from nodes*76 if demanded by memory requirements, or if OMP_NUM_THREADS>1.
# --array=1-1%1

###############################################################
#### Define job/run parameters and variables ###

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
mem_gb=$(( $SLURM_MEM_PER_NODE / 1000 ))

# Define job index for job arrays
#IDX=$SLURM_ARRAY_TASK_ID

# Load modules
module load angsd/0.940

#! Full path to application executable:
#application=

#! Run options for the application:
#options=

#! Work directory (i.e. where the job will run):
workdir="$SLURM_SUBMIT_DIR"

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this safe value to no more than 76:
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#! Number of MPI tasks to be started by the application per node and in total (do not change):
np=$[${numnodes}*${mpi_tasks_per_node}]

#! The following variables define a sensible pinning strategy for Intel MPI tasks -
#! this should be suitable for both pure MPI and hybrid MPI/OpenMP jobs:
export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets
#! Notes:
#! 1. These variables influence Intel MPI only.
#! 2. Domains are non-overlapping sets of cores which map 1-1 to MPI tasks.
#! 3. I_MPI_PIN_PROCESSOR_LIST is ignored if I_MPI_PIN_DOMAIN is set.
#! 4. If MPI tasks perform better when sharing caches/sockets, try I_MPI_PIN_ORDER=compact.

#! Uncomment one choice for CMD below (add mpirun/mpiexec options if necessary):

#! Choose this for a MPI code (possibly using OpenMP) using Intel MPI.
#CMD="mpirun -ppn $mpi_tasks_per_node -np $np $application $options"

#! Choose this for a pure shared-memory OpenMP parallel program on a single node:
#! (OMP_NUM_THREADS threads will be created):
#CMD="$application $options"

#! Choose this for a MPI code (possibly using OpenMP) using OpenMPI:
#CMD="mpirun -npernode $mpi_tasks_per_node -np $np $application $options"

###############################################################
#### RUN ###

# Change working directory
cd $workdir
echo -e "Changed directory to `pwd`.\n"

# Set parameters
REF=/cluster/work/gdc/shared/p936/secondNGSdataset/Assemblies/Dsylvestris/assembly_homozygous.fa

# Print job info
echo -e "JobID: ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}\n======"
echo "Time: $(date)"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"
if [ "$SLURM_JOB_NODELIST" ]; then
        echo -e "\nNodes allocated:\n================"
fi
echo -e "\nnumtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)"
echo -e "\nExecuting command:\n==================\n$CMD\n"

# Run command
# Calculate genotype likelihoods of focal populations
# Results are same whether including  -doPost 1 or not.

# Increase limit on number of open files
#ulimit -n 2048

# Change directory
out_dir="/cluster/work/gdc/shared/p936/DSYL_ms_revisions/GLF_files"
cd Recalibrated_and_clippedOverlapped_BAMs

# CEN (scaffold4_size532381:179265-204502)
angsd -GL 2 -out ${out_dir}/GL_196indsAlps_below1400m_CEN_region_w25kbFlanks_pVal10maf25depth2 -ref ${REF} -nThreads $SLURM_CPUS_PER_TASK -doGlf 2 -doMajorMinor 1 -minMaf 0.25 -SNP_pval 1e-10 -doMaf 1 -only_proper_pairs 0 -minMapQ 1 -minQ 1 -C 50 -remove_bads 1 -doCounts 1 -setMinDepth 392 -setMaxDepth 2940 -bam bam_196inds_centralAlps_below1400m.filelist -r scaffold4_size532381:154265-229502 # CEN region with 25kb flanking regions, all pops

# TCF1 (scaffold1_size1318325_1_614794:211891-222933)
#angsd -GL 2 -out ${out_dir}/GL_196indsAlps_below1400m_TCF1_region_w25kbFlanks_pVal10maf25depth2 -ref ${REF} -nThreads $SLURM_CPUS_PER_TASK -doGlf 2 -doMajorMinor 1 -minMaf 0.25 -SNP_pval 1e-10 -doMaf 1 -only_proper_pairs 0 -minMapQ 1 -minQ 1 -C 50 -remove_bads 1 -doCounts 1 -setMinDepth 392 -setMaxDepth 2940 -bam bam_196inds_centralAlps_below1400m.filelist -r scaffold1_size1318325_1_614794:186891-247933 # TCF1 region with 25kb flanking regions

#eval $CMD

###############################################################
##Get a summary of the job
myjobs -j ${SLURM_JOB_ID}

####job controlling, grep "Job:" *.out
echo "Job: ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} successfully finished $(date)"
# happy end
exit 0
###############################################################
