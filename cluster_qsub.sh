#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=4G,h_rt=6:00:00,tmem=4G
# tscratch is per rule, from cluster.yaml via submit.sh -- NOT fixed here.
# It is per slot like tmem, so the old flat 50G meant a 4-slot job needed
# 200G of local scratch, more than any node has (max ~190G), and such a job
# queues behind the one or two nodes that could ever run it.

# join stdout and stderr output
#$ -j y
#$ -sync y
#$ -R y

echo $JOB_ID

{exec_job}
