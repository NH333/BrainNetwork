#!/usr/bin/python
# -*- coding: utf-8 -*-       
#coding=utf-8
# coding: unicode_escape
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')#这两句用来解决plot.show()后没有反应的问题
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif']=['SimHei'] #用来正常显示中文标签
plt.rcParams['axes.unicode_minus']=False #用来正常显示负号


"""
读取原文件txt的fmri数据，存储为array类型
"""
fp = open(r"C:\\Users\\84241\\Desktop\\vscode_Python\\课题\\rmsts116_0050952.txt")
line = fp.readline()
data_list = []
while line:
    num = list(map(float,line.split())) #用map 把str类型转变成float 类型
    data_list.append(num)
    line = fp.readline()
fp.close()
data_array = np.array(data_list)

"""
计算功能联通矩阵FC
计算Pearson相关系数
并且把<0的值置为0
"""
Matrix_preprocess = data_array[40:130,:] #预处理 剔除前40和后40个time point
FC_Matrix = np.zeros((116,116))
uX = np.mean(Matrix_preprocess,axis=0)
uY = uX
tmp0 = 0
tmp1 = 0
tmp2 = 0
Pearson = 0
for i in range(len(uX)):
    for j in range(len(uY)):
        tmp0 = (Matrix_preprocess[:,i] - uX[i]).dot(Matrix_preprocess[:,j] - uY[j])
        for k in range(Matrix_preprocess.shape[0]): 
            tmp1 = tmp1 + (Matrix_preprocess[:,i][k]-uX[i])**2
            tmp2 = tmp2 + (Matrix_preprocess[:,j][k]-uY[j])**2
        Pearson = tmp0/np.sqrt(tmp1*tmp2)
        FC_Matrix[i,j] = Pearson
        tmp0,tmp1,tmp2,Pearson = 0,0,0,0
FC_Matrix[FC_Matrix<0] = 0  

np.savetxt("C:\\Users\\84241\\Desktop\\vscode_Python\\课题\\FC_Matrix.txt",FC_Matrix)
"""
绘制FC联通矩阵的热力图
"""
f, ax1 = plt.subplots(figsize=(6,6),nrows=1)
sns.heatmap(FC_Matrix,annot=False,ax=ax1,vmax=1,vmin=0,cmap='rainbow')
plt.show()

"""
计算结点的度
"""
NodeDegree = []
Threshold = 0.6
tmp = 0
for i in range(FC_Matrix.shape[1]):
    for j in range(FC_Matrix.shape[0]):
        if FC_Matrix[i,j] > Threshold:
            tmp = tmp + 1
    NodeDegree.append(tmp)
    tmp = 0  

np.savetxt("C:\\Users\\84241\\Desktop\\vscode_Python\\课题\\NodeDegree.txt",NodeDegree)

xlocal = np.linspace(1,len(NodeDegree),len(NodeDegree))  
plt.bar(xlocal[40:60],NodeDegree[40:60],width=0.4,label='节点的度',fc='r')
#在柱状图上添加具体值
for a,b in zip(xlocal[40:60],NodeDegree[40:60]): 
    plt.text(a, b+0.05, '%.0f' % b, ha='center') 
# plt.xlabel('40~59个脑区')
# plt.ylabel('Node_Degree')
# plt.show()


"""
计算聚类系数
"""
tmp_Node = []
tmp_num = 0
Cluster_num = [] 
FC_Matrix_Cluster = FC_Matrix.copy()
FC_Matrix_Cluster[FC_Matrix_Cluster<=0.6] = 0
for j in range(FC_Matrix_Cluster.shape[1]):
    for i in range(FC_Matrix_Cluster.shape[0]):
        if FC_Matrix_Cluster[i,j] == 0:
            continue
        else:
            tmp_Node.append(i)
    if len(tmp_Node) <= 2 :
        Cluster_num.append(0)
        tmp_Node = []
        continue
    else:
        tmp_Node.remove(j)
        for k in range(len(tmp_Node)-1):
            for z in range(k+1,len(tmp_Node)):
                if FC_Matrix_Cluster[tmp_Node[k],tmp_Node[z]] > 0:
                    tmp_num = tmp_num + 1
        Cluster_num.append(2*tmp_num / ((len(tmp_Node)-1) * (len(tmp_Node)))) 
        tmp_Node = []
        tmp_num = 0

np.savetxt("C:\\Users\\84241\\Desktop\\vscode_Python\\课题\\ClusterCofficient.txt",Cluster_num)

xlocal = np.linspace(1,len(Cluster_num),len(Cluster_num)) 
Cluster_num_copy = Cluster_num.copy()
for i in range(40,60):
    xlocal[i] = xlocal[i] + 0.4
    Cluster_num_copy[i] = Cluster_num_copy[i] * 18 

plt.bar(xlocal[40:60],Cluster_num_copy[40:60],width=0.4,label='聚类系数',tick_label=range(40,60),fc='y')
for a,b in zip(xlocal[40:60],Cluster_num[40:60]): 
    plt.text(a, b+b*18-0.8, '%.4f' % b, ha='center') 
plt.xlabel('40~59个脑区',fontsize=20)
plt.ylabel('黄色：聚类系数*max(节点的度) \n 红色：节点的度 ',fontsize=20)
plt.legend()
plt.show() 




input('')
