import argparse
import fileinput
import numpy as np
import array 
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
parser = argparse.ArgumentParser(description='Process plotting histogram')
parser.add_argument('filenames', nargs='+')
args = parser.parse_args()
fig = plt.figure()
for i in range(0, 3):
  filename = args.filenames[i]
  mylabel = filename.rsplit('.')[0] + ' ($\AA$)'
  rawdata=np.loadtxt(fileinput.input(args.filenames[i])).transpose()
  ax1= [1, 2, 3] 
  ax2= [1, 2, 3]
  ax1[i] = fig.add_subplot(1, 3, i+1)
  ax2[i] = ax1[i].twinx()
  ax1[i].bar(rawdata[0] ,rawdata[1]/1000000, edgecolor='k', width=0.1)
  ax2[i].plot(rawdata[0] ,rawdata[1]/1000000, 'r-')
  ax1[i].set_xlabel(mylabel , fontname='Times New Roman', fontsize=14)
  plt.subplots_adjust(wspace=0, hspace=0)
  v1=max(rawdata[1])
  v2=max(rawdata[1])
  ax11=ax1[i]
  ax21=ax2[i]
  def align_yaxis(ax1, v1, ax2,v2):
     """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
     _, y1 = ax11.transData.transform((0, v1))
     _, y2 = ax21.transData.transform((0, v2))
     inv = ax21.transData.inverted()
     _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
     miny, maxy = ax21.get_ylim()
     ax2.set_ylim(miny+dy, maxy+dy)
  align_yaxis(ax11, 0, ax21, 0) 
ax1[i].yaxis.set_ticklabels([])
ax1[i].yaxis.set_visible(False)
ax2[i].yaxis.set_visible(False)
ax1= fig.add_subplot(1, 3, 1)
ax1.set_ylabel(ur'Count$\times$10$^{6}$', fontname="Times New Roman", fontsize=14)
ax1= fig.add_subplot(1, 3, 2)
ax1.yaxis.set_visible(False)

plt.show()
