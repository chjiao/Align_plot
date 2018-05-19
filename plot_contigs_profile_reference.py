import re,sys,pdb
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

contig_file=sys.argv[1]
vcf_file=sys.argv[2]
loc_file = sys.argv[3]

# get contigs aligned to references
lineno=0
loc_dict={}
simi_dict = {}
align_dict = {}
with open(loc_file,'r') as f:
    for line in f:
       lineno+=1
       if lineno%5==1:
           lmap1=line.strip().split()
           con, ref, similarity = lmap1[0], lmap1[1], float(lmap1[2])
       elif lineno%5==2:
           lmap2=line.strip().split()
       elif lineno%5==3:
           lmap3=line.strip().split()
           # contig_1_13812  89.6_9669:      99.9
           # contig_1_1381   13812   1       9650
           # 89.6    9669    2       9651
           con_len, align_start, align_end = int(lmap2[1]), int(lmap2[2]), int(lmap2[3])
           if align_start>align_end:
               pdb.set_trace()
           align_len = abs(align_end - align_start+1)
           if int(lmap2[1])<500 or float(align_len)/con_len<0.5:
               continue
           #if not ref in loc_dict:
           #    loc_dict[ref] = [con]
           #else:
           #    loc_dict[ref].append(con)

           ## only look at contigs that are fully aligned
           if not con in simi_dict and float(align_len)/con_len>0.7:
               simi_dict[con] = similarity
               align_dict[con] = ref
           elif align_len==con_len and simi_dict[con]<similarity:
               simi_dict[con] = similarity
               align_dict[con] = ref


for con in align_dict:
    ref = align_dict[con]
    if not ref in loc_dict:
        loc_dict[ref] = [con]
    else:
        loc_dict[ref].append(con)

count = 0
for ref in loc_dict:
    count+=len(loc_dict[ref])
print "The total number of aligned contigs is: ",count

pdb.set_trace()

contig_dict = {}
with open(contig_file, 'r') as f:
    for line in f:
        if line.startswith('>'):
            title=line[1:].strip().split()[0]
        else:
            seq=line.strip()
            assert not title in contig_dict, "Contig aready exist!"+'\t'+title
            contig_dict[title] = seq

for con in contig_dict:
    if not con in align_dict:
        print con

pdb.set_trace()

## 
con_profile_dict = {}
for con in contig_dict:
    con_profile = np.zeros(len(contig_dict[con]))
    con_profile_dict[con] = con_profile

with open(vcf_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        lmap=line.strip().split('\t')
        con,pos,info = lmap[0],lmap[1],lmap[7]
        m=re.search('DP=(\d+);',info)
        depth=int(m.group(1))
        con_profile_dict[con][int(pos)-1] = depth

## plot
num_contigs = len(contig_dict)
for ref in loc_dict:
    contigs = loc_dict[ref]
    num_figures = len(contigs)
    print num_figures
    plt.figure(figsize=(10, 10*num_figures))
    fig_idx = 1
    for con in contigs:
        assert con in con_profile_dict, "Contig is not in the profile!"
        plt.subplot(num_figures, 1, fig_idx)
        plt.plot(con_profile_dict[con])
        plt.title(con)
        plt.xlabel('Contig Position')
        plt.ylabel('Depth')
        fig_idx+=1
    figname = ref[:-1]+'.pdf'
    plt.savefig(figname,format='pdf',dpi=200)
    plt.close()
