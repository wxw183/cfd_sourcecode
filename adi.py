import numpy as np
import matplotlib.pyplot as plt
import csv
from pylab import mpl
# 设置显示中文字体
mpl.rcParams["font.sans-serif"] = ["Microsoft YaHei"]

mpl.rcParams["axes.unicode_minus"] = False

#与时间步长关系
epsilon=np.zeros(4)

infile = open('04/epsilon_1', 'r')
for line in infile:
    epsilon[0] = float(line)
infile.close()

infile = open('05/epsilon_1', 'r')
for line in infile:
    epsilon[1] = float(line)
infile.close()

infile = open('06/epsilon_1', 'r')
for line in infile:
    epsilon[2] = float(line)
infile.close()

infile = open('07/epsilon_1', 'r')
for line in infile:
    epsilon[3] = float(line)
infile.close()

logdt=[np.log(0.005),np.log(0.01),np.log(0.02),np.log(0.04)]
fig=plt.figure()
ax = fig.add_axes((0.1, 0.1, 0.8, 0.8),facecolor="#e1e1e1")
ax.plot(logdt,epsilon)
ax.set_xlabel("$log\Delta t$")
ax.set_ylabel("$log\epsilon$")
ax.set_title("误差与时间步长关系")
formatter = mpl.ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-1,1))
ax.xaxis.set_major_formatter(formatter)
plt.savefig("误差与时间步长关系.pdf")

#与空间步长关系,时间步长Δt = 0.02， 网格间距为 x = 0.1，0.05，0.025的均匀网格
epsilon_1=np.zeros(3)

infile = open('03/epsilon_1', 'r')
for line in infile:
    epsilon_1[0] = float(line)
infile.close()

infile = open('02/epsilon_1', 'r')
for line in infile:
    epsilon_1[1] = float(line)
infile.close()

infile = open('01/epsilon_1', 'r')
for line in infile:
    epsilon_1[2] = float(line)
infile.close()

logdx=[np.log(0.025),np.log(0.05),np.log(0.1)]

fig=plt.figure()
ax = fig.add_axes((0.1, 0.1, 0.8, 0.8),facecolor="#e1e1e1")
ax.plot(logdx,epsilon_1)
ax.set_xlabel("$log\Delta x$")
ax.set_ylabel("$log\epsilon$")
ax.set_title("误差与空间步长关系")
formatter = mpl.ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-1,1))
ax.xaxis.set_major_formatter(formatter)
plt.savefig("误差与空间步长关系.pdf")




#与物理时间关系
fig=plt.figure()
ax = fig.add_axes((0.1, 0.1, 0.8, 0.8),facecolor="#e1e1e1")

infile=open("01/eh.csv",'r')
table=[]
for row in csv.reader(infile):
    table.append(row)
infile.close()
t=np.zeros(len(table))
e=np.zeros(len(table))
for r in range(0,len(table)):
    for c in range(0,len(table[0])):
        table[r][c]=float(table[r][c])
    t[r]=table[r][0]
    e[r]=table[r][1]
ax.plot(t,e,label="$\Delta x=0.1,\Delta y=0.02$",linewidth=0.4)

infile=open("02/eh.csv",'r')
table=[]
for row in csv.reader(infile):
    table.append(row)
infile.close()
t=np.zeros(len(table))
e=np.zeros(len(table))
for r in range(0,len(table)):
    for c in range(0,len(table[0])):
        table[r][c]=float(table[r][c])
    t[r]=table[r][0]
    e[r]=table[r][1]


ax.plot(t,e,label="$\Delta x=0.05,\Delta y=0.02$",linewidth=0.4)

infile=open("03/eh.csv",'r')
table=[]
for row in csv.reader(infile):
    table.append(row)
infile.close()
t=np.zeros(len(table))
e=np.zeros(len(table))
for r in range(0,len(table)):
    for c in range(0,len(table[0])):
        table[r][c]=float(table[r][c])
    t[r]=table[r][0]
    e[r]=table[r][1]
ax.plot(t,e,label="$\Delta x=0.025,\Delta y=0.02$",linewidth=0.4)

infile=open("04/eh.csv",'r')
table=[]
for row in csv.reader(infile):
    table.append(row)
infile.close()
t=np.zeros(len(table))
e=np.zeros(len(table))
for r in range(0,len(table)):
    for c in range(0,len(table[0])):
        table[r][c]=float(table[r][c])
    t[r]=table[r][0]
    e[r]=table[r][1]
ax.plot(t,e,label="$\Delta x=0.05,\Delta y=0.005$",linewidth=0.4)

infile=open("05/eh.csv",'r')
table=[]
for row in csv.reader(infile):
    table.append(row)
infile.close()
t=np.zeros(len(table))
e=np.zeros(len(table))
for r in range(0,len(table)):
    for c in range(0,len(table[0])):
        table[r][c]=float(table[r][c])
    t[r]=table[r][0]
    e[r]=table[r][1]
ax.plot(t,e,label="$\Delta x=0.05,\Delta y=0.01$",linewidth=0.4)

infile=open("06/eh.csv",'r')
table=[]
for row in csv.reader(infile):
    table.append(row)
infile.close()
t=np.zeros(len(table))
e=np.zeros(len(table))
for r in range(0,len(table)):
    for c in range(0,len(table[0])):
        table[r][c]=float(table[r][c])
    t[r]=table[r][0]
    e[r]=table[r][1]
ax.plot(t,e,label="$\Delta x=0.05,\Delta y=0.02$",linewidth=0.4)

infile=open("07/eh.csv",'r')
table=[]
for row in csv.reader(infile):
    table.append(row)
infile.close()
t=np.zeros(len(table))
e=np.zeros(len(table))
for r in range(0,len(table)):
    for c in range(0,len(table[0])):
        table[r][c]=float(table[r][c])
    t[r]=table[r][0]
    e[r]=table[r][1]
ax.plot(t,e,label="$\Delta x=0.05,\Delta y=0.04$",linewidth=0.4)


ax.set_xlabel("t(s)");ax.set_ylabel("$log_{10}E_h$")
ax.set_title("误差与时间关系曲线")
ax.legend()
plt.savefig("误差与时间关系曲线.pdf")







