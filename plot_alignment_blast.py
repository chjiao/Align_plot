import re,sys,pdb
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

# plot the alignment results for blastn
# adjust the annotation text for the blastn alignment results are more dispersed
# 2016.08.16

loc_file=sys.argv[1]

lineno=0
loc_dict={}
with open(loc_file,'r') as f:
    for line in f:
       lineno+=1
       if lineno%5==1:
           lmap1=line.strip().split()
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
           align_len = abs(align_end - align_start)
           if int(lmap2[1])<500 or float(align_len)/con_len<0.5:
               continue
           if not lmap1[1] in loc_dict:
               loc_dict[lmap1[1]]=[(lmap1,lmap2,lmap3)]
           else:
               loc_dict[lmap1[1]].append((lmap1,lmap2,lmap3))

## plot 
plt.figure(figsize=(14,35))
index=0
for genome in loc_dict:
    index+=1
    plt.subplot(len(loc_dict),1,index)
    contig_list=loc_dict[genome]
    con_idx=0
    
    wid_start_list=[]
    wid_end_list=[]
    gen_name = ''
    gen_len = 0
    for con in contig_list:
        con_name=con[0][0]
        gen_name=con[0][1]
        simi=con[0][2]
        con_len=int(con[1][1])
        con_al=int(con[1][2]) # left coordinate
        con_ar=int(con[1][3]) # right coordinate
        gen_len=int(con[2][1])
        gen_al=int(con[2][2]) # left
        gen_ar=int(con[2][3]) # right
        if gen_al>gen_ar:
            tmp = gen_al
            gen_al = gen_ar
            gen_ar = tmp

        con_start=1-(con_al-gen_al)
        con_end=con_start+con_len
        wid_start_list.append(con_start)
        wid_end_list.append(con_end)
    min_start = min(wid_start_list)
    if min_start>0:
        min_start = 0
    max_end = max(wid_end_list)


    #if con_len>=500:
    plt.plot((1,gen_len),(0.2,0.2),'k-',linewidth=3.0)
    plt.text(-1000+min_start,0.25,gen_name,horizontalalignment='center',verticalalignment='center',fontsize=8)
    plt.text(1-200,0.5,str(1),horizontalalignment='center',verticalalignment='center',fontsize=6,color='g') # left coordinate
    plt.text(gen_len+200,0.5,str(gen_len),horizontalalignment='center',verticalalignment='center',fontsize=6) # right coordinate
    xticks = range(-1000,10001,1000)
    plt.xticks(xticks, rotation='vertical')

    for con in contig_list:
        con_idx+=1
        con_name=con[0][0]
        gen_name=con[0][1]
        simi=con[0][2]
        con_len=int(con[1][1])
        con_al=int(con[1][2]) # left coordinate
        con_ar=int(con[1][3]) # right coordinate
        gen_len=int(con[2][1])
        gen_al=int(con[2][2]) # left
        gen_ar=int(con[2][3]) # right
        if gen_al>gen_ar:
            tmp = gen_al
            gen_al = gen_ar
            gen_ar = tmp

        con_start=1-(con_al-gen_al)
        con_end=con_start+con_len

        con_align_len=con_ar-con_al+1
        con_align_start=gen_al
        con_align_end=con_align_start+con_align_len-1
        con_align_mid = con_align_start+con_align_len/2

        plt.plot((con_start,con_end),(con_idx,con_idx),'y-',linewidth=3.0)
        plt.plot((con_align_start,con_align_end),(con_idx,con_idx),'b-',linewidth=4.0)
        plt.plot((con_align_start,con_align_end),(0.2,0.2),'b-',linewidth=4.0)  # plot on the genome
        plt.text(con_align_start-200,con_idx+0.2,str(con_al)+'('+str(con_align_start)+')',horizontalalignment='center',verticalalignment='center',fontsize=6,color='g') # left coordinate
        plt.text(con_align_end+200,con_idx+0.2,str(con_ar)+'('+str(con_align_end)+')',horizontalalignment='center',verticalalignment='center',fontsize=6) # right coordinate
        plt.text(con_align_mid-50,con_idx+0.2, str(simi)+'%',horizontalalignment='center',verticalalignment='center',fontsize=6,color='r', weight='bold')
        plt.text(-1000+min_start,con_idx,con_name,horizontalalignment='center',verticalalignment='center',fontsize=8)
        #if float(simi)>=99.0:
        #    plt.text(max_wid+1000,con_idx,simi,horizontalalignment='center',verticalalignment='center',fontsize=8,color='r',weight="bold")
        #else:
        #    plt.text(max_wid+1000,con_idx,simi,horizontalalignment='center',verticalalignment='center',fontsize=8)
    plt.grid()
figname='contig_alignments_test.png'
plt.savefig(figname,format='png',dpi=300)
plt.close()

