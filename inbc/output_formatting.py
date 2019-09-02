import sys
file_name = sys.argv[1]
c = 0
fw = open('Outputfile.txt', 'a+')
with open(file_name) as f:
    for line in f:
        taxid, seq = line.strip().split(',')
        taxid = taxid.split(' ')[1]
        readid = seq.split('/')[-1].split('.')[0]
        try:
            output = '{}\t{}\n'.format(readid, taxid)
            fw.write(output)
        except:
            c += 1
fw.close()
print('wrong number: {}'.format(c), flush = True)
