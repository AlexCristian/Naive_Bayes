year = '2010'
import numpy as np
import glob
from sklearn.metrics import accuracy_score
pref = np.zeros((5, 3))
for ix in range(1, 6):
    file = 'results/' + year + '_j' + str(ix) + '.txt'
    pair = []
    with open(file) as f:
        for line in f:
            if line[:6] == 'Genome':
                label = int(line.split(' ')[3][:-1])
                pred = int(line.split(' ')[-1])
                pair.append((label, pred))

    tmp = glob.glob('/lustre/scratch/zz374/NBC/process/temporal_data/' + year + '/save_' + str(ix) + '/*.dat')
    labeled_genomes = [int(item.split('/')[9].split('-')[0]) for item in tmp]
    labeled_genomes = set(labeled_genomes)
    
    data = np.zeros((len(pair), 3))
    for i in range(len(pair)):
        data[i, 0] = pair[i][0]
        data[i, 1] = pair[i][1]
        if pair[i][0] in labeled_genomes: 
            data[i, 2] = 1
        else:
            data[i, 2] = 0
    pref[ix - 1, 0] = accuracy_score(data[:, 0], data[:, 1])
    pref[ix - 1, 1]	= accuracy_score(data[data[:, 2] == 1, 0], data[data[:, 2] == 1, 1])
    pref[ix - 1, 2]	= data[data[:, 2] == 1].shape[0]/data.shape[0]
    # print('overall_accuracy:', accuracy_score(data[:, 0], data[:, 1]))
    # print('labeled_accuracy:', accuracy_score(data[data[:, 2] == 1, 0], data[data[:, 2] == 1, 1]))
    # print('labeled__percent:', data[data[:, 2] == 1].shape[0]/data.shape[0])
print('overall_accuracy:')
for i in range(5):
    print(pref[i, 0])
print('labeled__percent:')
for i in range(5):
    print(pref[i, 2])
print('labeled_accuracy:')
for i in range(5):
    print(pref[i, 1])

