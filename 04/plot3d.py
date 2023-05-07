import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import csv

# 设置显示中文字体

infile=open("dataout1.csv",'r')
import csv
table=[]
for row in csv.reader(infile):
    table.append(row)
infile.close()
u=np.eye(len(table[0])-2)
for r in range(len(table[0])-2):
    for c in range(len(table[0])-2):
        u[r][c]=float(table[r][c])



fig, ax = plt.subplots(figsize=(8, 6))

ax.set_title("")
ax.set_xlabel("$x$", fontsize=16)
ax.set_ylabel("$y$", fontsize=16)

x=np.linspace(0,1,len(table[0])-1)
y=np.linspace(0,1,len(table[0])-1)

X,Y=np.meshgrid(x,y)

norm = mpl.colors.Normalize(abs(u).min(), abs(u).max())
p = ax.pcolormesh(X, Y, u, norm=norm, cmap=mpl.cm.bwr)
ax.axis('tight')
ax.set_xlabel(r"$y$", fontsize=18)
ax.set_ylabel(r"$x$", fontsize=18)
ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(4))
ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(4))
cb = fig.colorbar(p, ax=ax)
cb.set_label(r"$u$", fontsize=18)

fig.savefig("pcolor.png")
