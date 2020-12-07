#alignment
import os
import random
fastq_folder = '../PM01/65hrs_7_2-40591052'
os.chdir('/data/home/llong35/data/CB4856_genome')
files = os.listdir(fastq_folder)
for file in files:
    if not file.endswith('.fastq.gz'):
        continue
    fastq_file = os.path.join(fastq_folder,file)
    sam_file = os.path.join(fastq_folder,file.split('.')[0]+'.sam')
    command = 'bwa mem CB4856_genome.fa '+fastq_file + ' > '+sam_file
    print(command)
    os.system(command)