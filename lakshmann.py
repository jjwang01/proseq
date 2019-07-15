import os
os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"]=""

import pandas as pd
import numpy as np
import pyBigWig
k562_positive=pyBigWig.open('/oak/stanford/groups/akundaje/laks/proseq_chromputer/K562_unt.sort.bed.gz_plus.bw')
k562_negative=pyBigWig.open('/oak/stanford/groups/akundaje/laks/proseq_chromputer/K562_unt.sort.bed.gz_minus.bw')
h3k27ac_bigwig=pyBigWig.open('/oak/stanford/groups/akundaje/laks/proseq_chromputer/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig')
h3k27ac_peak_replicated_bed=pd.read_csv('/oak/stanford/groups/akundaje/laks/proseq_chromputer/wgEncodeBroadHistoneK562H3k27acStdPk.broadPeak',sep='\t',header=None)
total_sum=0
for i in range(0,22):
	if i==0:
		chrm='chrX'
	else:
		chrm='chr'+str(i)
	length=k562_positive.chroms(chrm)
	values=k562_positive.values(chrm,0,length)
	a=np.array(values)
	a[np.isnan(a)] = 0
	total_sum=total_sum+np.sum(a)

total_sum_negative=0
for i in range(0,22):
	if i==0:
		chrm='chrX'
	else:
		chrm='chr'+str(i)
	length=k562_negative.chroms(chrm)
	values=k562_negative.values(chrm,0,length)
	a=np.array(values)
	a[np.isnan(a)] = 0
	total_sum_negative=total_sum_negative+np.sum(a)

h3k27ac_peak_replicated_bed['diff']=h3k27ac_peak_replicated_bed[2]-h3k27ac_peak_replicated_bed[1]
h3k27ac_peak_replicated_bed['absolute_diff']=10000-h3k27ac_peak_replicated_bed['diff']
import math
h3k27ac_peak_replicated_bed['left']=h3k27ac_peak_replicated_bed['absolute_diff'].apply(lambda x : math.floor(float(x)/2))
h3k27ac_peak_replicated_bed['right']=h3k27ac_peak_replicated_bed['absolute_diff'].apply(lambda x : math.ceil(float(x)/2))
h3k27ac_peak_replicated_bed['left_new']=h3k27ac_peak_replicated_bed[1]-h3k27ac_peak_replicated_bed['left']
h3k27ac_peak_replicated_bed['right_new']=h3k27ac_peak_replicated_bed[2]+h3k27ac_peak_replicated_bed['right']
h3k27ac_peak_replicated_bed['new_diff']=h3k27ac_peak_replicated_bed['right_new']-h3k27ac_peak_replicated_bed['left_new']
h3k27ac_peak_replicated_bed['left_new']=h3k27ac_peak_replicated_bed['left_new'].astype(int)
h3k27ac_peak_replicated_bed['right_new']=h3k27ac_peak_replicated_bed['right_new'].astype(int)
proseq_positive_values=k562_positive.values("chr11",837429,847429)
proseq_negative_values=k562_negative.values("chr11",837429,847429)
h3k27ac_values=h3k27ac_bigwig.values("chr11",837429,847429)
import matplotlib.pyplot as plt
fig, ax1 = plt.subplots()
ax1.set_xlabel('genome')
ax1.set_ylabel('proseq read counts')
ax1.plot(proseq_positive_values/((np.abs(total_sum_negative)+total_sum ) / 1000000),color='blue',alpha=0.5,label='proseq +ve')
ax1.plot(proseq_negative_values/((np.abs(total_sum_negative)+total_sum ) / 1000000),color='orange',alpha=0.5,label='proseq -ve')
plt.show()
plt.close()

