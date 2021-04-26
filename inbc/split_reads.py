import sys
file_name = sys.argv[1]
out_dir = sys.argv[2]
with open(file_name) as f:
    tmp_string = ''
    header = ''
    out_name = ''
    keep_flag = False
    for line in f:
        if line[0] == '>':
            new_header = line.strip()
            new_header = new_header.replace('/', '-')
       	    new_header = new_header.replace('>', '')
            index = new_header.split('-')[1]
            new_header = new_header.split('-')[0]
            if index == '2':
                if new_header == header:
                    tmp_string += '>' + new_header + '-2\n'
                    keep_flag = True
                else:
                    keep_flag = False
                continue
            if tmp_string != '':
                fw = open(output_name, 'w+')
                fw.write(tmp_string)
                fw.close()
            header = new_header
            output_name = out_dir + '/' + header + '.fna'
            tmp_string = '>' + header + '-1\n'
            keep_flag = True
        else:
            if keep_flag:
               tmp_string += line
 
    if tmp_string != '':
        fw = open(output_name, 'w+')
        fw.write(tmp_string)
        fw.close()

