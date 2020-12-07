import subprocess,argparse,os
from Bio import SeqIO
import pysam
import pdb
import scipy.stats


N2_CB4856_chrom = {'CHROMOSOME_I': 'gi|809001836|gb|CM003206.1|',
              'CHROMOSOME_II': 'gi|809001828|gb|CM003207.1|',
              'CHROMOSOME_III': 'gi|809001817|gb|CM003208.1|',
              'CHROMOSOME_IV': 'gi|809001806|gb|CM003209.1|',
              'CHROMOSOME_V': 'gi|809001797|gb|CM003210.1|',
              'CHROMOSOME_X': 'gi|809001786|gb|CM003211.1|',
              'CHROMOSOME_MtDNA': 'gi|809001771|gb|CM003212.1|'}

N2_CB4856_chrom = {'CHROMOSOME_I': 'I',
              'CHROMOSOME_II': 'II',
              'CHROMOSOME_III': 'III',
              'CHROMOSOME_IV': 'IV',
              'CHROMOSOME_V': 'V',
              'CHROMOSOME_X': 'X',
              'CHROMOSOME_MtDNA': 'MtDNA'}
nucleotides=['A','T','C','G']         

class gene():
    def __init__(self, gene_name, public_name):
        self.gene_name = gene_name
        self.public_name = public_name
        self.SNV_list = []
        self.N2_read_list = []
        self.CB4856_read_list = []
     
    def add_SNV(self, SNV):
        self.SNV_list.append(SNV)
        
    def merge_read_list(self):
        read_list = {}
        for SNV in self.SNV_list:
            for read in SNV.N2_read_list:
                try:
                    read_list[read]+=1
                except:
                    read_list[read]=1
        for SNV in self.SNV_list:
            for read in SNV.CB4856_read_list:
                try:
                    read_list[read]-=1
                except:
                    read_list[read]=-1
        return read_list
    
    def merge_N2_read_list(self):
        for SNV in self.SNV_list:
            self.N2_read_list += SNV.N2_read_list
        self.N2_read_list = list(set(self.N2_read_list))
    
    def merge_CB4856_read_list(self):
        for SNV in self.SNV_list:
            self.CB4856_read_list += SNV.CB4856_read_list
        self.CB4856_read_list = list(set(self.CB4856_read_list))
    
    def count_N2_CB4856_reads(self):
        read_list = self.merge_read_list()
        N2 = 0
        CB4856 = 0
        for k,v in read_list.items():
            if v > 0:
                N2 += 1
            elif v < 0:
                CB4856 += 1
        return N2, CB4856
    
    def cal_Binom_p(self):
        self.merge_N2_read_list()
        self.merge_CB4856_read_list()
        N2_count = len(self.N2_read_list)
        CB4856_count = len(self.CB4856_read_list)
        binom_p = scipy.stats.binom_test(N2_count,N2_count+CB4856_count)
        self.binom_p = binom_p
        return binom_p
        
            
        

class SNV():
    def __init__(self,chrom, N2_position,CB4856_position,N2_base,CB4856_base,feature):
        self.chrom = chrom
        self.N2_position = N2_position
        self.CB4856_position = CB4856_position
        self.N2_base = N2_base
        self.CB4856_base = CB4856_base
        self.feature = feature
        self.N2_read_list = []
        self.CB4856_read_list = []
    
    def add_N2_read(self, read_name):
        if read_name not in self.N2_read_list:
            self.N2_read_list.append(read_name)

    
    def add_CB4856_read(self, read_name):
        if read_name not in self.CB4856_read_list:
            self.CB4856_read_list.append(read_name)

    
    
parser = argparse.ArgumentParser()
HOME_DIR='/data/home/llong35/data/PM01/65hrs_7_2-40591052/'
parser.add_argument('-b_N2','--bam_file_mapped_to_N2', type = str,help='bam_file mapped to N2',default=HOME_DIR+'65hrs-7-2_S1_L001.bam')
parser.add_argument('-b_CB4856','--bam_file_mapped_to_CB4856', type = str,help='bam_file mapped to CB4856',default=HOME_DIR+'65hrs-7-2_S1_L001_R1_001.bam')
parser.add_argument('-a','--snp_annotation_file_between_CB4856_WS230', type = str,help='snp_annotation_file_between_CB4856_WS230',default='/data/home/llong35/data/CB4856_genome/SNV_N2_CB4856')
parser.add_argument('-o_1','--output_file', help='output file of gene summary',default='/data/home/llong35/data/AE_output/7_2_1')
# parser.add_argument('-r','--reference_files', type = str,help='reference file the bam file mapped to')
args = parser.parse_args()
# pdb.set_trace()
try:
    output_f = open(args.output_file,'w')
except:
    output_f = open('snp_info_verification_output','w')

f = open(args.snp_annotation_file_between_CB4856_WS230,'r')
bamfile_N2 = pysam.AlignmentFile(args.bam_file_mapped_to_N2)
bamfile_CB4856 = pysam.AlignmentFile(args.bam_file_mapped_to_CB4856)
i=0
gene_name_list = []
gene_list = []
for line in f:
        if line.startswith('#'):
            continue
        i+=1
        if i%500==0:
            print(str(i))
        chrom, N2_position, N2_nucleotide, CB4856_position, CB4856_nucleotide = line.split()[0:5]
        feature = line.split()[6]
        SNV_object = SNV(chrom,N2_position,CB4856_position,N2_nucleotide,CB4856_nucleotide,feature)
        gene_name, public_gene_name = line.split()[8:10]
        if gene_name == 'NA':
            continue
        if gene_name not in gene_name_list:
            gene_object = gene(gene_name, public_gene_name)
            gene_name_list.append(gene_name)
            gene_list.append(gene_object)
        chrom = 'CHROMOSOME_'+chrom
        N2_position = int(N2_position)-1
        CB4856_position = int(CB4856_position)-1
        pileups = bamfile_N2.pileup(chrom,N2_position,N2_position+1)
        try:
            for column in pileups:
                if column.pos == N2_position:
                    break
            for read in column.pileups:
                if not read.is_del and not read.is_refskip and read.alignment.mapq != 0:
                    alignment = read.alignment
                    base = alignment.seq[read.query_position]
                    if base == N2_nucleotide:
                        SNV_object.add_N2_read(alignment.qname)
        except:
            pass
        pileups = bamfile_CB4856.pileup(N2_CB4856_chrom[chrom],CB4856_position,CB4856_position+1)
        try:
            for column in pileups:
                if column.pos == CB4856_position:
                    break
            for read in column.pileups:
                if not read.is_del and not read.is_refskip and read.alignment.mapq != 0:
                    alignment = read.alignment
                    base = alignment.seq[read.query_position]
                    if base == CB4856_nucleotide:
                        SNV_object.add_CB4856_read(alignment.qname)
        except:
            pass
        
        for gene_object in gene_list:
            if gene_object.gene_name == gene_name:
                break
        gene_object.add_SNV(SNV_object)

# output_f.write('gene_name,public_name,N2_specific_read_count, CB4856_specific_read_count,binom_p\n')
# for gene in gene_list:
#     file_name = '/Volumes/Lijiang_data/datasets/Tajima\'s_D_region/Mar_14_data/' + gene.gene_name+'_'+ gene.public_name+'.csv'
#     output_file = open(file_name,'w')
#     binom_p = gene.cal_Binom_p()
#     for SNV in gene.SNV_list:
#         SNV_info = ''
#         SNV_info += SNV.chrom+','+SNV.N2_position+','+SNV.N2_base+','+SNV.CB4856_position+','+SNV.CB4856_base+','
# #         +':feature:'+SNV.feature
#         SNV_info += str(len(SNV.N2_read_list))+','
# #        ':N2_base:'+SNV.N2_base+':N2_specific_read_count:'+
#         SNV_info += str(len(SNV.CB4856_read_list))+','
# #         'CB4856_base:'+SNV.CB4856_base+':CB4856_specific_read_count:'+
#         SNV_info += SNV.feature +'\n'
#         output_file.write(SNV_info)
#     output_file.close()
#     output_f.write(gene.gene_name+','+gene.public_name+','+str(len(gene.N2_read_list))+','+str(len(gene.CB4856_read_list))+','+str(gene.binom_p)+'\n')
output_f.write('gene_name,public_name,N2_specific_read_count, CB4856_specific_read_count\n')


for gene in gene_list:
    N2_count,CB4856_count = gene.count_N2_CB4856_reads()
    output_f.write(gene.gene_name+','+gene.public_name+','+str(N2_count)+','+str(CB4856_count)+'\n')


output_f.close()
f.close()

        
        
        
                
        

