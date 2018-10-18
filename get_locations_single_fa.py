import re,sys,subprocess,pdb
from itertools import groupby

# Calculate the sequence similarity using Water from EMBOSS
# Group fasta files before pairwise alignment
# Compare the sequences from multiple sequences in single fasta files
# 2018.03.14

def get_command(program):
    if program=='needle':
        command="needle -asequence 01.fa -bsequence 02.fa -gapopen 10 -gapextend 0.5 -outfile temp.needle"
        output = "temp.needle"
    elif program=='water':
        command="water -asequence 01.fa -bsequence 02.fa -outfile temp.water"
        output = "temp.water"
    else:
        command="blastn -query 01.fa -subject 02.fa -out temp.blastn"
        output = "temp.needle"
    return command, output

def process_output(temp_out, program, f_out):
    with open(temp_out, 'r') as f:
        for line in f:
            f_out.write(line)

fa_file=sys.argv[1]
program=sys.argv[2]
out_file=sys.argv[3]
simi = 90

fa_dict1={}
title_list1=[]
contig_index=0
with open(fa_file,'r') as f:
    faiter=(x[1] for x in groupby(f,lambda line:line[0]=='>'))
    for header in faiter:
        contig_index+=1
        header=header.next().strip()
        seq="".join(s.strip() for s in faiter.next())

        if header.startswith('>'):
            title=header[1:].split()[0]
            #title='contig_'+str(contig_index)+'_'+str(len(seq))
            fa_dict1[title]=seq
            title_list1.append(title)
        else:
            print "Error!"
            pdb.set_trace()

f_out=open(out_file,'w')
##
for i in range(len(title_list1)-1):
    for j in range(i+1, len(title_list1)):
        f1=open('01.fa','w')
        f2=open('02.fa','w')
        seq1=fa_dict1[title_list1[i]]
        seq2=fa_dict1[title_list1[j]]
        f1.write('>'+title_list1[i]+'\n'+seq1)
        f2.write('>'+title_list1[j]+'\n'+seq2)
        f1.close()
        f2.close()
        command, temp_out = get_command(program)
        #command="blastn -query 01.fa -subject 02.fa -out temp.blastn"

        subprocess.call(command, shell=True)  #call will wait until the command completes
        #pdb.set_trace()

        process_output(temp_out, program, f_out)

f_out.close()
subprocess.call('rm 01.fa', shell=True)
subprocess.call('rm 02.fa', shell=True)
subprocess.call('rm '+temp_out, shell=True)

