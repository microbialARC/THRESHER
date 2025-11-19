# Adapted from CoreSimul (https://github.com/lbobay/CoreSimul):
# - extract_names.py
# - branch_length.py
# Note: Original variable names and code structure have been preserved from the source.
# This may result in naming conventions and coding style that differ from THRESHER standards.


import os
import sys


# Codes from extract_names.py
intermediate_dir = snakemake.params.intermediate_dir
simu_tree = snakemake.input.simu_tree

# Import input from snakemake
TREE = simu_tree

path = intermediate_dir

f=open(TREE,"r")
l=f.readline()

symbols=[",","(",")",":",";"]

link={}
bag=[]
tag=0
l=l.strip("\n")
i=1
while i < len(l):
	j=i-1
	if l[i] in symbols:
		if tag == 1:
			bag.append(taxa)
			link[taxa] = souvenir + taxa
			tag=0
	else:
		if l[j] == "(" or l[j] == ",":
			souvenir=l[j]
			taxa = l[i]
			tag=1
		else:
			taxa+=l[i]
	i+=1


print(bag)


memo=l


h=open(path + "names.txt","w")
nb=0
for stuff in bag:
	if 1==1:
		if 2==2:
			taxa = stuff
			#print(taxa)
			nb+=1
			if nb < 10:
				new = "BOZO000" + str(nb) + ":"
			elif nb < 100:
				new = "BOZO00" + str(nb)  + ":"
			elif nb < 1000:
				new = "BOZO0" + str(nb)  + ":"
			elif nb < 10000:
				new = "BOZO" + str(nb)  + ":"
			h.write(new.strip(":") + "\t" + taxa + "\n")
			taxa = taxa + ":"
			new = link[taxa.strip(":")][0] + new 
			thingy = link[taxa.strip(":")] + ":"
			memo = memo.replace(thingy,new)

h.close()

print("Found",nb,"taxa in the tree")


new=""
tmp=""
tag=0
for L in memo:
	if L == ")":
		tag=1
		tmp+=L
	if tag==0:
		new+=L
	elif tag==1:
		if L==":":
			tmp+=L
			new+=tmp
			tag=0
			tmp=""

new += tmp 
new+=";"

h=open(path + "renamed.tree","w")
h.write(new)
h.close()


# Codes from branch_length.py
# Import input from snakemake
TREE = os.path.join(intermediate_dir , "renamed.tree")
path = intermediate_dir

forest=[TREE]


for file in forest:
	unique=[]
	#print(file)
	new_file = file.split("/")[-1]
	new_file = new_file.split(".tree")[0] + ".txt"
	h=open(path + new_file,"w")
	f=open(path + file ,"r")
	l=f.readline()
	l = l.strip("\n")
	l = l.strip(";")
	#print l
	f.close()

	numbers=["0","1","2","3","4","5","6","7","8","9"]
	resu=""
	i=0
	while i < len(l):
		L = l[i]
		if l[i] in numbers:
			if i ==0:
				resu+=l[i]	
			elif l[i-1] == ")" or l[i-2] == ")" or l[i-3] == ")":
				pass
			else:
				resu+=l[i]
		else:
			resu+= l[i]
		i+=1


	resu = l


	I = 1
	
	dico={}
	detail={}
	arbre={}
	arbre[I]=resu
	first={}
	level={}
	length={}
	node=0
	while arbre[I].count("(") > 1:
		groupe=[]
		memo=[]
		resu=str(arbre[I])
		#print arbre[I]
		resu2 = resu
		#node=0
		i=0
		mark=0
		while i < len(resu):
			#print(resu)
			k = resu[i]
			if k == "B" or k == "n":
				name=k
				mark=0
			elif k != "," and k !="("  and k !=")" and k !=":":
				if mark == 0:
					name += k
				elif mark == 1:
					branch += k
			elif k == ':':
				branch = ''
				mark=1
			else:
				if k=="," or k ==')': 
					length[name] = branch
					if 'B' in name:
						if name not in unique:
							#print name
							h.write(name + "\t" + name + '\t' + str(abs(float(branch))) + '\ttip\n')
							unique.append(name)
					mark=0
				if k==",":
					#print name,' ',branch
					if resu[i+1] == "B" or  resu[i+1] =="n":
						if name not in memo:
							groupe=[name]
						else:
							groupe=[]
				if k == ")":
					if name not in groupe:
						if name not in memo:
							if len(groupe) == 1:
								groupe.append(name)
								memo.extend(groupe)
								node += 1
								NODE = "n" + str(I) + "_" + str(node)
								#print(NODE,groupe)
								mono = "(" + groupe[0] + ':' + length[groupe[0]] + "," + groupe[1] + ':' + length[groupe[1]] + ")" 
								if mono in resu2:
									level[NODE] = I
									#print NODE," ",resu2
									resu2 = resu2.replace(mono,NODE)
									I+=1
									dico[NODE] = [groupe[0] , groupe[1]]
									if "B" not in groupe[0]:
										detail[NODE] = list(detail[groupe[0]])
									else:
										detail[NODE] =[groupe[0]]
									if "B" not in groupe[1]:
										for stuff in detail[groupe[1]]:
											detail[NODE].append(stuff)
									else:
										detail[NODE].append(groupe[1])
									i=0
									break
								else:
									node = node - 1
			i+=1

		#I+=1
		arbre[I] = resu2
	



	MAX_I = int(I)

	#print l.strip("\n")
	#print "MAX_I ", arbre[MAX_I]
	
	ROOTS=[]
	nb=0
	if ":" in arbre[MAX_I]:
		a=arbre[MAX_I].strip("(").strip(")").split(",")
		for truc in a:
			nb+=1
			#if "S" in truc:
			sub_a = truc.split(":")
			#print "sub_a"," ",sub_a
			name=sub_a[0]
			#print file," root" + str(nb) + ": ",name
			if name not in unique:
				out = truc
				detail["root" + str(nb)] = [name]
				length[name] = sub_a[1].strip("(").strip(")")
				ROOTS.append(name)
				if name.startswith("B"):
					h.write(name + "\t" +  name + "\t" +  sub_a[1].strip("(").strip(")")  + "\t" + "tip" +  "\n")
			else:
				ROOTS.append(name)
	
	
	
	#h=open("../results/mono/core/mono_" + SP + ".txt","w")
	parent={}
	for NODE in detail:
		detail[NODE].sort()
		composit = NODE.split("_")[0] + "_" + "-".join(detail[NODE])
		try:
			LENGTH= length[NODE] 
		except KeyError:
			LENGTH='root'
		if "-".join(detail[NODE]) not in unique:
			#print "TAG2 ","-".join(detail[NODE])
			if LENGTH=="root":
				if ":" in "-".join(detail[NODE]):
					finish = "-".join(detail[NODE]).split(":")[0]
					#print("Finish= ",finish)
					if NODE in  ROOTS:
						tag="root"
					else:
						tag="branch"
					h.write(NODE + "\t" +  finish + "\t" +  LENGTH  + "\t" + tag +  "\n")
			else:
				if NODE in  ROOTS:
					tag="root"
				else:
					tag="branch"
				h.write(NODE.split("_")[0] + "\t" + "-".join(detail[NODE]) + "\t" +  LENGTH + "\t" + tag  + "\n")
			unique.append("-".join(detail[NODE]))
		for st in detail[NODE]:
			if st in parent:
				parent[st].append(NODE)
			else:
				parent[st] = [NODE]
	
	h.close()


	#print "ROOTS= ",ROOTS


	out_name = file.split(".tree")[0]
	h=open(path + "dichotomies.txt","w")
	for st in dico:
		#print st," ",dico[st]
		resu1,resu2,resu3=st,dico[st][0],dico[st][1]
		resu1,resu2,resu3=resu1.split("_")[0],resu2.split("_")[0],resu3.split("_")[0]
		h.write(resu1 + "\t" + resu2 + "\t" + resu3 + "\n")
	h.close()
	
	h=open(path + "/roots.txt","w")
	for root in ROOTS:
		h.write(root.split("_")[0] + "\n")
	h.close()
