import re
import numpy as np 
import csv
s = '2.65^{+0.96}_{-0.67}'
ex = re.findall('\{(.*?)\}',s)
ex = (float(ex[0])-float(ex[1]))/2
x = s[:s.index("^")]
print(x,ex)
clusters=[]
with open('tables/table-1.csv','r') as f:
    read = csv.reader(f)
    key = next(read)
    for row in read:
        clusters.append(row)

for cluster in clusters:
    name=cluster[0]

    x=float((cluster[2][:cluster[2].index("^")])) 
    err_x = (float(re.findall('\{(.*?)\}',cluster[2])[0])-float(re.findall('\{(.*?)\}',cluster[2])[1]))/2
    T=[T,err_T]
    
    x=float((cluster[3][:cluster[3].index("^")])) 
    err_x = (float(re.findall('\{(.*?)\}',cluster[3])[0])-float(re.findall('\{(.*?)\}',cluster[3])[1]))/2
    beta=[x,err_x]
    
    x=float((cluster[4][:cluster[4].index("^")]))
    err_x = (float(re.findall('\{(.*?)\}',cluster[4])[0])-float(re.findall('\{(.*?)\}',cluster[4])[1]))/2
    r_c=[x,err_x]    

    x=float((cluster[5][:cluster[5].index("^")]))*(1e4)
    err_x = (float(re.findall('\{(.*?)\}',cluster[5])[0])-float(re.findall('\{(.*?)\}',cluster[5])[1]))/2*(1e4)
    n_c=[x,err_x]    

    t_cool=float((cluster[6][:cluster[6].index("^")]))
    if t_cool<1.3:
        continue

    # n_c = float(cluster[5])*(10**4)
    # r_500 = float(cluster[-1])*1000