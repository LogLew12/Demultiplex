#!/usr/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mail-user='llewis3@uoregon.edu'
#SBATCH --mail-type=END,FAIL

/usr/bin/time -v python demux.py \
-r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz \
-r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz \
-i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz \
-i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz \
-b /projects/bgmp/shared/2017_sequencing/indexes.txt \
-c 35 \
# -r1 /projects/bgmp/llewis3/bioinformatics/Bi622/demux/unit_test/test_R1.fastq \
# -r2 /projects/bgmp/llewis3/bioinformatics/Bi622/demux/unit_test/test_R2.fastq \
# -i1 /projects/bgmp/llewis3/bioinformatics/Bi622/demux/unit_test/test_I1.fastq \
# -i2 /projects/bgmp/llewis3/bioinformatics/Bi622/demux/unit_test/test_I2.fastq \