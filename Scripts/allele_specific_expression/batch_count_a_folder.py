#batch count a folder
import os


dir = '/data/home/llong35/data/PM01/65hrs_7_2-40591052'
for i in [1]:
    end = 'L00'+str(i)+'.bam'
    for file in os.listdir(dir):
        if not file.endswith(end):
            continue
        file_path = os.path.join(dir,file)
        CB4856_file = file.split('.')[0]+'_R1_001.bam'
        CB4856_path = os.path.join(dir,CB4856_file)
        output_file = file.split('.')[0]+'_results.csv'
        output_path = os.path.join(dir,output_file)
        script_file='assign_each_read_to_each_allele_of_each_gene.py'
        command = 'python '+ script_file + ' -b_N2 '+file_path+' -b_CB4856 '+CB4856_path + ' -o_1 '+output_path
        print(command)
        os.system(command)
