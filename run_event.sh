#! /bin/bash
#SBATCH --time=3:00:00


use root

RUNNR=$(printf %06d $SLURM_ARRAY_TASK_ID)

cd /vol/astro3/lofar/lora/LORA_software_V2

python LORA_by_event.py -i $RUNNR
