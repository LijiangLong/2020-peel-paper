# sam to bam in each folder
import os
import random
sam_folder = '../PM01/65hrs_7_2-40591052'
os.chdir('/data/home/llong35/data/WS230')
files = os.listdir(sam_folder)
for file in files:
    if not file.endswith('sam'):
        continue
    sam_file = os.path.join(sam_folder,file)
    bam_file = os.path.join(sam_folder,file.split('.')[0]+'.bam')
    if os.path.exists(bam_file):
        continue
    random_number = random.randrange(9999999999999)
    command = 'samtools sort -O bam -T temp'+str(random_number)
    command += ' '+sam_file + ' >| ' + bam_file
    print(command)
    os.system(command)
    
