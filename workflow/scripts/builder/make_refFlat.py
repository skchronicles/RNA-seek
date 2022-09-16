
#while read a b c d e f g h i j k l;do
#a1=`grep -m1 $d hs38d1.annotate.isoforms.txt|awk '{print $2}'`;
#i1=`python add2all.py $l $b|sed "s/ //g"`
#j1=`python sum_lists.py $i1 $k|sed "s/ //g"`
#echo -e "$a1\t$d\t$a\t$f\t$b\t$c\t$g\t$h\t$i1\t$j1"
#done< genes.ref.bed|head
import sys
def add2all(a,b):
	lst=a.strip().split(",")
	lst.pop(-1)
	return ",".join([str(int(x)+int(b)) for x in lst])+","
def sum_lists(a,b):
	lst1=a.strip().split(",")
	lst1.pop(-1)
	lst2=b.strip().split(",")
	lst2.pop(-1)
	return ",".join([str(int(a)+int(b)) for a,b in zip(lst1,lst2)])+","

tid2genename={a[0]: a[1] for a in map(lambda x:x.strip().split("  "),open("annotate.isoforms.txt").readlines())}

for l in map(lambda x:x.strip().split("\t"),open("genes.ref.bed").readlines()):
	newl=[]
	newl.append(tid2genename["\""+l[3]+"\""].strip("\""))
	newl.append(l[3])
	newl.append(l[0])
	newl.append(l[5])
	newl.append(l[1])
	newl.append(l[2])
	newl.append(l[6])
	newl.append(l[7])
	newl.append(l[9])
	lst=add2all(l[11],l[1])
	newl.append(lst)
	newl.append(sum_lists(l[10],lst))
	print("\t".join(newl))
