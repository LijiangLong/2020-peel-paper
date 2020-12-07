import pdb
import pandas as pd
import os


def main():
    read_count_file = '/Users/lijiang/Downloads/Reads_count.xlsx'
    ae_folder = '/Users/lijiang/Desktop/results'
    output_file = '/Users/lijiang/Desktop/summary_results.csv'
    df = pd.read_excel(read_count_file)
    sample_count = {}
    
    for i in range(len(df)):
        sample,replicate,count = df.iloc[i,:]
        try:
            sample_count[sample] += count
        except KeyError:
            sample_count[sample] = count
    record = {}
    pdb.set_trace()
    for file in os.listdir(ae_folder):
        if not file.endswith('.csv'):
            continue
        key = file[6:9]
        if not key in record:
            record[key] = {}
        df = pd.read_csv(os.path.join(ae_folder,file))
        for i in range(len(df)):
            gene_name,public_name,n2_read,cb_read = df.iloc[i,:]
            gene_key = ','.join([gene_name,public_name])
            if not gene_key in record[key]:
                record[key][gene_key]  = {'n2':0,'cb':0}
            record[key][gene_key]['n2'] += n2_read
            record[key][gene_key]['cb'] += cb_read
    with open(output_file,'w') as output_f:
        output_f.write('name,public_name,sample,n2_read,cb4856_read\n')
        for sample in record:
            for gene in record[sample]:
                n2_read = record[sample][gene]['n2']
                n2_read = n2_read/sample_count[sample]*1000000
                cb4856_read = record[sample][gene]['cb']
                cb4856_read = cb4856_read/sample_count[sample]*1000000
                output_f.write(','.join([gene,sample,'{:2f}'.format(n2_read),'{:2f}'.format(cb4856_read)]))
                output_f.write('\n')
        


if __name__ == "__main__":
    main()