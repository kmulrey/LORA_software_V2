#! /bin/bash
#SBATCH --time=3:00:00


source /vol/astro7/lofar/hpandya/root/my_build/bin/thisroot.sh




RUNNR=$(printf %06d $SLURM_ARRAY_TASK_ID)

cd /vol/astro3/lofar/lora/LORA_software_V2/LORA_software_V2

python LORA_by_event_debug.py -i $RUNNR
