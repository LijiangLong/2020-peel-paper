import numpy as np
import pdb
import pandas as pd
import matplotlib.pyplot as plt

def selection_coefficient_calculator(genotype_freq_file):
    output_file = 'temp.csv'
    df = pd.read_csv(genotype_freq_file)
    with open(output_file,'w') as output_f:
        output_f.write('p,A,a\n')
        for i in range(len(df)-1):
            this_generation = df.iloc[i,:]
            p_this_gen = this_generation[0]+this_generation[1]/2
            this_gen = np.array([p_this_gen,1-p_this_gen])
            
            next_generation = df.iloc[i+1,:]
            p_next_gen = next_generation[0]+next_generation[1]/2
            next_gen = np.array([p_next_gen,1-p_next_gen])
            relative_freq = next_gen/this_gen
#             max_outcome = max(relative_freq)
            max_outcome = relative_freq[1]
            fitness = relative_freq/max_outcome
            if p_this_gen > 0.9999999:
                break
            output_f.write(str(p_this_gen))
            for f in fitness:
                output_f.write(','+str(f))
            output_f.write('\n')
    return output_file
            
def plot_selection_coefficient(allele_fitness_file):
    df = pd.read_csv(allele_fitness_file)
    x = df.p
    y = 1-df.A
#     _ = plt.plot(x,y,color='black')
#     _ = plt.ylim(0,1)
#     _ = plt.xlim(0,1)
#     _ = plt.savefig('temp.svg', dpi=300)
    return x,y
def main():
    _ = plt.rc('axes.spines', **{'bottom': True, 'left': True, 'right': False, 'top': False})
    genotype_freq_files = []
    genotype_freq_files.append('beneficial.csv')

#     genotype_freq_files.append('s_05_start_high.csv')
#     genotype_freq_files.append('s_05_start_low.csv')
#     genotype_freq_files.append('s_1_start_high.csv')
#     genotype_freq_files.append('s_1_start_low.csv')
#     genotype_freq_files.append('fitness_cost_0_trajectory.csv')
#     genotype_freq_files.append('fitness_cost_0.5_trajectory.csv')
#     genotype_freq_files.append('fitness_cost_1_trajectory.csv')
    
    for genotype_freq_file in genotype_freq_files:
        allele_fitness_file = selection_coefficient_calculator(genotype_freq_file)
        x,y = plot_selection_coefficient(allele_fitness_file)
        _ = plt.plot(x,y,label=genotype_freq_file)
    _ = plt.legend(bbox_to_anchor=(1.05,1.05))
    _ = plt.xlim(0,1)
#     _ = plt.ylim(0,1)
    plt.show()
    pdb.set_trace()
    
    

#     _ = plt.savefig('beneficial_allele_selection_coefficient.png', dpi=300)


if __name__ == "__main__":
    main()
