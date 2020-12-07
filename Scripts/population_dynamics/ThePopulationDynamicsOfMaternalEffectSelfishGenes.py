import numpy as np
import pdb
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
_ = plt.rc('axes.spines', **{'bottom': True, 'left': True, 'right': False, 'top': False})

def three_genetypes_onegeneration(mm,mo,oo,s=0.1,h=1,t=1):
    new_mm,new_mo,new_oo =  0,0,0
    # family 1
    new_mm += (1-s)*mm*mm
    # family 2
    new_mm += 0.5 * (1-s) * mo * mm
    new_mo += 0.5 * (1-s) * mo * mm
    # family 3
    new_mo += (1-s) * oo * mm
    # family 4
    new_mm += 0.5 * (1-h*s) * mm * mo
    new_mo += 0.5 * (1-h*s) * mm * mo
    # family 5
    new_mm += 0.25 * (1-h*s) * mo * mo
    new_mo += 0.5 * (1-h*s) * mo * mo
    new_oo += 0.25 * (1-h*s) * mo * mo * (1-t)
    # family 6
    new_mo += 0.5 * (1-h*s) * oo * mo
    new_oo += 0.5 * (1-h*s) * oo * mo * (1-t)
    # family 7
    new_mo += mm * oo
    # family 8
    new_mo += 0.5 * mo * oo
    new_oo += 0.5 * mo * oo
    # family 9
    new_oo += oo * oo

    new_total = new_mm + new_mo + new_oo
    return new_mm/new_total,new_mo/new_total,new_oo/new_total

def partial_selfing_maternal(mm,mo,oo,k=0.5,s=0.1,h=1,t=1):
    new_mm, new_mo, new_oo = 0, 0, 0
    # mm selfing
    new_mm += (1-k) * (1-s) * mm
    # mo selfing
    new_mm += (1-k) * (1-h*s) * mo * 0.25
    new_mo += (1-k) * (1-h*s) * mo * 0.5
    # oo selfing
    new_oo += (1-k) * oo
    # family 1
    new_mm += (1-s)*mm*mm * k
    # family 2
    new_mm += 0.5 * (1-s) * mo * mm* k
    new_mo += 0.5 * (1-s) * mo * mm* k
    # family 3
    new_mo += (1-s) * oo * mm* k
    # family 4
    new_mm += 0.5 * (1-h*s) * mm * mo* k
    new_mo += 0.5 * (1-h*s) * mm * mo* k
    # family 5
    new_mm += 0.25 * (1-h*s) * mo * mo* k
    new_mo += 0.5 * (1-h*s) * mo * mo* k
    new_oo += 0.25 * (1-h*s) * mo * mo * (1-t)* k
    # family 6
    new_mo += 0.5 * (1-h*s) * oo * mo* k
    new_oo += 0.5 * (1-h*s) * oo * mo * (1-t)* k
    # family 7
    new_mo += mm * oo* k
    # family 8
    new_mo += 0.5 * mo * oo* k
    new_oo += 0.5 * mo * oo* k
    # family 9
    new_oo += oo * oo* k

    new_total = new_mm + new_mo + new_oo
    return new_mm/new_total,new_mo/new_total,new_oo/new_total


def partial_selfing_paternal(mm,mo,oo,k=0.5,s=0.1,h=1,t=1):
    new_mm, new_mo, new_oo = 0, 0, 0
    # mm selfing
    new_mm += (1-k) * (1-s) * mm
    # mo selfing
    new_mm += (1-k) * (1-h*s) * mo * 0.25
    new_mo += (1-k) * (1-h*s) * mo * 0.5
    # oo selfing
    new_oo += (1-k) * oo
    # family 1
    new_mm += (1-s)*mm*mm * k
    # family 2
    new_mm += 0.5 * (1-s) * mo * mm* k
    new_mo += 0.5 * (1-s) * mo * mm* k
    # family 3
    new_mo += (1-s) * oo * mm* k
    # family 4
    new_mm += 0.5 * (1-h*s) * mm * mo* k
    new_mo += 0.5 * (1-h*s) * mm * mo* k
    # family 5
    new_mm += 0.25 * (1-h*s) * mo * mo* k
    new_mo += 0.5 * (1-h*s) * mo * mo* k
    new_oo += 0.25 * (1-h*s) * mo * mo * (1-t)* k
    # family 6
    new_mo += 0.5 * (1-h*s) * oo * mo* k
    new_oo += 0.5 * (1-h*s) * oo * mo * k
    # family 7
    new_mo += mm * oo* k
    # family 8
    new_mo += 0.5 * mo * oo* k
    new_oo += 0.5 * mo * oo* k* (1-t)
    # family 9
    new_oo += oo * oo* k

    new_total = new_mm + new_mo + new_oo
    return new_mm/new_total,new_mo/new_total,new_oo/new_total

def partial_selfing_paternal_fitness(mm,mo,oo,k=0.5,s=0.1,h=1,t=1):
    new_mm, new_mo, new_oo = 0, 0, 0
    # mm selfing
    new_mm += (1-k) * (1-s) * mm
    # mo selfing
    new_mm += (1-k) * (1-h*s) * mo * 0.25
    new_mo += (1-k) * (1-h*s) * mo * 0.5
    # oo selfing
    new_oo += (1-k) * oo
    # family 1
    new_mm += (1-s)*mm*mm * k
    # family 2
    new_mm += 0.5 * (1-h*s) * mo * mm* k
    new_mo += 0.5 * (1-h*s) * mo * mm* k
    # family 3
    new_mo += oo * mm* k
    # family 4
    new_mm += 0.5 * (1-s) * mm * mo* k
    new_mo += 0.5 * (1-s) * mm * mo* k
    # family 5
    new_mm += 0.25 * (1-h*s) * mo * mo* k
    new_mo += 0.5 * (1-h*s) * mo * mo* k
    new_oo += 0.25 * (1-h*s) * mo * mo * (1-t)* k
    # family 6
    new_mo += 0.5  * oo * mo* k
    new_oo += 0.5  * oo * mo * k
    # family 7
    new_mo += (1-s)*mm * oo* k
    # family 8
    new_mo += (1-h*s)*0.5 * mo * oo* k
    new_oo += (1-h*s)*0.5 * mo * oo* k* (1-t)
    # family 9
    new_oo += oo * oo* k

    new_total = new_mm + new_mo + new_oo
    return new_mm/new_total,new_mo/new_total,new_oo/new_total

def partial_selfing_paternal_both_fitness(mm,mo,oo,k=0.5,s=0.1,h=1,t=1):
    new_mm, new_mo, new_oo = 0, 0, 0
    # mm selfing
    new_mm += (1-k) * (1-s) * mm* (1-s)
    # mo selfing
    new_mm += (1-k) * (1-h*s) * mo * 0.25* (1-h*s)
    new_mo += (1-k) * (1-h*s) * mo * 0.5* (1-h*s)
    # oo selfing
    new_oo += (1-k) * oo
    # family 1
    new_mm += (1-s)*mm*mm * k * (1-s)
    # family 2
    new_mm += 0.5 * (1-h*s) * mo * mm* k* (1-s)
    new_mo += 0.5 * (1-h*s) * mo * mm* k* (1-s)
    # family 3
    new_mo += oo * mm* k* (1-s)
    # family 4
    new_mm += 0.5 * (1-s) * mm * mo* k* (1-h*s)
    new_mo += 0.5 * (1-s) * mm * mo* k* (1-h*s)
    # family 5
    new_mm += 0.25 * (1-h*s) * mo * mo* k* (1-h*s)
    new_mo += 0.5 * (1-h*s) * mo * mo* k* (1-h*s)
    new_oo += 0.25 * (1-h*s) * mo * mo * (1-t)* k* (1-h*s)
    # family 6
    new_mo += 0.5  * oo * mo* k* (1-h*s)
    new_oo += 0.5  * oo * mo * k* (1-h*s)
    # family 7
    new_mo += (1-s)*mm * oo* k
    # family 8
    new_mo += (1-h*s)*0.5 * mo * oo* k
    new_oo += (1-h*s)*0.5 * mo * oo* k* (1-t)
    # family 9
    new_oo += oo * oo* k

    new_total = new_mm + new_mo + new_oo
    return new_mm/new_total,new_mo/new_total,new_oo/new_total

def beneficial_allele(mm,mo,oo,s=0.1,h=0.5):
    new_mm = 1*mm
    new_mo = (1-h*s)*mo
    new_oo = (1-s)*oo

    new_total = new_mm + new_mo + new_oo
    return new_mm/new_total,new_mo/new_total,new_oo/new_total


def allele_freq(mm,mo,oo):
    return mm+mo/2

def one_simulation(mm=0,mo=0.1,oo=0.9,k=0.9, s=0.75, h=0, t=1,gen=100):
    freq = []

    for i in range(gen):
        mm,mo,oo = partial_selfing_paternal(mm, mo, oo, k, s, h, t)
        af = allele_freq(mm,mo,oo)
        freq.append(af)
#         if af == 1 or af == 0:
#             return freq
#         try:
#             if freq[-1] == freq[-2] and freq[-2] == freq[-3] and freq[-3] == freq[-4]:
#                 return freq
#         except:
#             pass
    return freq
    
    
def genotype_tracked_simulation(mm=0,mo=0.1,oo=0.9,k=0.9, s=0.75, h=0, t=1,gen=100):
    freq = []
    for i in range(gen):
        mm,mo,oo = partial_selfing_paternal(mm, mo, oo, k, s, h, t)
        freq.append([mm,mo,oo])
        if mm > 0.999:
            break
    return freq




def selection_coefficient_curve(k=1, s=0, h=0.5, t=1,effect='peel'):
    ms = np.arange(0,1,0.01)
    ys=  []
    for m in ms:
        o=1-m
        if k == 1:
            mm=m*m
            mo=2*m*o
            oo=o*o
        else:
            sigma = 1-k
            mm = m*m+sigma*m*o/2/(1-sigma/2)
            mo = 2*m*o-2*(sigma*m*o/2/(1-sigma/2))
            oo = o*o + sigma*m*o/2/(1-sigma/2)
        mm,mo,oo = partial_selfing_paternal(mm, mo, oo, k, s, h, t)
        if effect == 'peel': 
            mm,mo,oo = partial_selfing_paternal(mm, mo, oo, k, s, h, t)
        elif effect == 'beneficial':
            mm,mo,oo = beneficial_allele(mm, mo, oo, s, h)
        af = allele_freq(mm,mo,oo)
#         ys.append(af-m)


        previous = np.array([m,1-m])
        now = np.array([af,1-af])
        relative = now/previous
        relative = relative/relative[0]
        ys.append(1-relative[1])
        
    return ms,ys

def allele_freq_change_curve(k=1, s=0, h=0.5, t=1,effect='peel'):
    iterations = 1
    ms = np.arange(0.0001,1,0.01)
    ys=  []
    for m in ms:
        o=1-m
        if k == 1:
            mm=m*m
            mo=2*m*o
            oo=o*o
        else:
            sigma = 1-k
            mm = m*m+sigma*m*o/2/(1-sigma/2)
            mo = 2*m*o-2*(sigma*m*o/2/(1-sigma/2))
            oo = o*o + sigma*m*o/2/(1-sigma/2)
        if effect == 'peel': 
            for _ in range(iterations):
                mm,mo,oo = partial_selfing_paternal(mm, mo, oo, k, s, h, t)
        elif effect == 'beneficial':
            for _ in range(iterations):
                mm,mo,oo = beneficial_allele(mm, mo, oo, s, h)
        af = allele_freq(mm,mo,oo)
        ys.append((af-m)/iterations)
    return ms,ys

def allele_freq_change_from_file(input_file):
    df = pd.read_csv(input_file)
    xs = []
    ys = []
    for i in range(len(df)-1):
        this_generation = df.iloc[i,:]
        p_this_gen = this_generation[0]+this_generation[1]/2
        this_gen = np.array([p_this_gen,1-p_this_gen])
            
        next_generation = df.iloc[i+1,:]
        p_next_gen = next_generation[0]+next_generation[1]/2
        next_gen = np.array([p_next_gen,1-p_next_gen])
        change = p_next_gen - p_this_gen
        xs.append(p_this_gen)
        ys.append(change)
    return xs, ys


def h_and_s_curve(k=1, s=0, h=0.5, t=1,effect='peel'):
    ms = np.arange(0.001,1,0.001)
    hs=  []
    ss = []
    
    for m in ms:
        o = 1 - m
        if k == 1: 
            mm=m*m
            mo=2*m*o
            oo=o*o
        else:
            sigma = 1-k
            mm = m*m+sigma*m*o/2/(1-sigma/2)
            mo = 2*m*o-2*(sigma*m*o/2/(1-sigma/2))
            oo = o*o + sigma*m*o/2/(1-sigma/2)
        if effect == 'peel': 
            new_mm,new_mo,new_oo = partial_selfing_paternal(mm, mo, oo, k, s, h, t)
        elif effect == 'beneficial':
            new_mm,new_mo,new_oo = beneficial_allele(mm, mo, oo, s, h)
        previous = np.array([mm,mo,oo])
        now = np.array([new_mm,new_mo,new_oo])
        relative = now/previous
        relative = relative/relative[0]
        s = 1-relative[2]
        h = (1-relative[1])/s
        ss.append(s)
        hs.append(h)
        
    return ms,hs,ss

def main():
    
    for s in [0,-0.18]:
        for k in [1,0.5]:
            freq = one_simulation(mm=0,mo=0.05,oo=0.95,k=k, s=s, h=0.5, t=1,gen=100)
            plt.plot(range(100),freq,label=str(s)+' '+str(k))
    _ = plt.ylim(0,1)
    _ = plt.xlim(0,100)
    _ = plt.legend()
#     _ = plt.show()
#     
#     pdb.set_trace()
    _ = plt.savefig('added_finness_effect.svg', dpi=300)
    
    
    
    
    
    
# def main():
# 
#     output_file = 'beneficial_new.csv'
#     gen = 1000
#     m=0.01
#     o=1-m
#     mm=m*m
#     mo=2*m*o
#     oo=o*o
#     
#     with open(output_file,'w') as output_f:
#         output_f.write('mm,mo,oo\n')
#         for i in range(gen):
#             mm,mo,oo = beneficial_allele(mm,mo,oo,s=0.5,h=0.5)
#             output_f.write(','.join([str(x) for x in [mm,mo,oo]]))
#             output_f.write('\n')


# def main():
#     input_file = 'beneficial_new.csv'
#     
#     xs,ys = allele_freq_change_from_file(input_file)
#     _  = plt.plot(xs,ys)
#     plt.show()
#     pdb.set_trace()
            
# def main():
#     ms,ys = allele_freq_change_curve(k=0.1, s=0, h=0.5, t=1)
#     _ = plt.plot(ms,ys,label='k=0.9')
#     ms,ys = allele_freq_change_curve(k=1, s=0.3, h=0.5, t=1)
#     _ = plt.plot(ms,ys,label='s=0.3')
#     
#     ms,ys = allele_freq_change_curve(k=1, s=0.6, h=0.5, t=1)
#     _ = plt.plot(ms,ys,label='s=0.6')
#     ms,ys = selection_coefficient_curve(k=1, s=0.69, h=0.5, t=0)
#     ys = [-1*y for y in ys]
#     _ = plt.plot(ms,ys,label='s=0.69,t=0')
#     ms,ys = selection_coefficient_curve(k=1, s=1, h=0.5, t=0)
#     ys = [-1*y for y in ys]
#     _ = plt.plot(ms,ys,label='s=1,t=0')
#     ms,ys = selection_coefficient_curve(k=1, s=0.2, h=0.5, t=0)
#     ys = [-1*y for y in ys]
#     _ = plt.plot(ms,ys,label='s=0.2,t=0')
#     for k in np.arange(0.1,1.1,0.1):
#     ks = [1,0.1,0.001,0.0001]
#     max_ys = []
#     pdb.set_trace()
#     for k in ks:
#     for s in np.arange(0,1.01,0.5):
#     for s in np.logspace(-4,0,50):
#         ms,ys = allele_freq_change_curve(k=1, s=s, h=0.5, t=0,effect='beneficial')
#         max_ys.append(max(ys))
# #         _ = plt.plot(ms,ys,label='k={k:.2f}'.format(k=k))
#         
#     _ = plt.plot(np.logspace(-4,0,50),max_ys,c='black')
#     ms,ys = allele_freq_change_curve(k=1, s=0.44, h=0.5, t=0,effect='beneficial')
#     _ = plt.plot(ms,ys,label='beneficial allele')
#     _ = plt.legend(bbox_to_anchor=(1.01,1.05))
#     _ = plt.plot((0,1),(0,0),linestyle=':')
#     _ = plt.ylim(0.00001,0.1)
#     _ = plt.xticks(np.arange(1,5,1))
#     _ = plt.xlim(0,1)
#     _ = plt.yscale("log")
#     _ = plt.xscale("log")
#     _ = plt.show()
#     pdb.set_trace()
    
    
#     _ = plt.savefig('beneficial_allele_max_allele_freq_change.svg', dpi=300)

# def main():
#     
#     ms,hs,ss = h_and_s_curve(k=0.1, s=0, h=0.5, t=1,effect='peel')
#     _ = plt.plot(ms,hs,label='h')
#     _ = plt.plot(ms,ss,label='s')
#     _ = plt.plot((0,1),(0,0),linestyle=':')
#     _ = plt.ylim(-2,1)
#     plt.legend()
# #     plt.show()
# #     pdb.set_trace()
# 
#     
#     _ = plt.savefig('h_s_change_10%_outcrossing.png', dpi=300)


if __name__ == "__main__":
    main()



