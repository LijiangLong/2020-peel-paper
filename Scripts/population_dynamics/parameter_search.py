from thepopulationdynamicsofmaternaleffectselfishgenes import *
import numpy as np
import seaborn as sns
import pdb
# 
# def search():
#     gen = 1000
#     _ = plt.rc('axes.spines', **{'bottom': True, 'left': True, 'right': False, 'top': False})
#     fig = plt.figure(figsize=(10, 8))
# 
#     axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#     xs = []
#     ys = []
#     cs = []
#     for s in np.arange(-0.02,1.0001,0.03):
#         for k in np.arange(0.01,1.0001,0.03):
#             freq = one_simulation(mm=0, mo=0.1, oo=0.9, k=k, s=s, h=0.4, t=1, gen=gen)
#             # if len(freq) < gen and freq[-1] > 0 and freq[-1] < 1:
#             #     cs.append('blue')
#             if freq[-1] > freq[0] and freq[-1] < 0.9999:
#                 cs.append('blue')
#             else:
#                 cs.append('red')
#             xs.append(s)
#             ys.append(k)
# 
#     _ = plt.scatter(xs,ys,c=cs,s=15)
#     _ = plt.xlabel('fitness cost s')
#     _ = plt.ylabel('outcrossing rate k')
#     # plt.legend()
#     plt.ylim(0,1.02)
#     # plt.show()
#     # pdb.set_trace()
#     plt.savefig('h0.4.svg', dpi=300)

# 
def search():
    gen = 10000
    _ = plt.rc('axes.spines', **{'bottom': True, 'left': True, 'right': False, 'top': False})
    fig = plt.figure(figsize=(10, 8))

    axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    matrix= []
    m=0.5
    o=1-m
    mm=m*m
    mo=2*m*o
    oo=o*o
    for k in np.arange(0.02,-0.00001,-0.0002):
#     for k in np.logspace(-4,-2,num=10):
        finals = []
        for s in np.arange(0,0.020001,0.0002):
#         for s in np.logspace(-4,-2,num=10):
            freq = one_simulation(mm=mm, mo=mo, oo=oo, k=k, s=s, h=0, t=1, gen=gen)
            finals.append(freq[-1])
        matrix.append(finals)

    
#     _ = plt.scatter(xs,ys,c=cs,s=15)
#     _ = plt.xlabel('fitness cost s')
#     _ = plt.ylabel('outcrossing rate k')
    # plt.legend()
#     plt.ylim(0,1.02)
    _ = sns.heatmap(matrix)
    _ = plt.xticks(np.arange(0,101,20))
    _ = plt.yticks(np.arange(0,101,20))
#     _ = plt.show()

    plt.savefig('parameter_space_heatmap_starting50_enlarged_2.svg', dpi=300)
#     pdb.set_trace()
search()