import subprocess
import shlex
import glob
import re
import os


class seq_comp:
	################################################################ Function for Data Extraction fron Tab-separated File #################################################################
	def data_extract(self, infile, deletes):
		f = open(infile, 'r')
		lines = f.readlines()
		dict_array = {}
		if (len(lines) > 2):
			split_l = ((lines[3]).rstrip()).split('\t')
		else:
			split_l = ((lines[1]).rstrip()).split('\t')
		for i in range(len(split_l)):
			dict_array['arr_'+str(i)] = []
		for l in lines[deletes:]:
			split_l = (l.rstrip()).split('\t')
			for x in range(len(split_l)):
				dict_array['arr_'+str(x)].append(split_l[x])
		f.close(); #print(dict_array);
		if ('.ptt' in infile):
			genus = ((lines[0]).split())[0]; #print(genus);
			return (genus, dict_array)
		else:
			return dict_array

	def data_wrangle(self, inputs, deletes):
		dict_array = {}
		lines = inputs
		split_l = ((lines[3]).rstrip()).split('\t')
		for i in range(len(split_l)):
			dict_array['arr_'+str(i)] = []
		for l in lines[deletes:]:
			split_l = (l.rstrip()).split('\t')
			for x in range(len(split_l)):
				dict_array['arr_'+str(x)].append(split_l[x])
		#print(dict_array);
		return dict_array

	################################################################ Function to find Reverse Complement of any Genomic sequence ###########################################################
	def reverse_complement(self, seq):
		complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'Y':'R', 'R':'Y', 'W':'W', 'S':'S', 'K':'M', 'M':'K', 'D':'H', 'H':'D', 'V':'B', 'B':'V', 'X':'X', 'N':'N', '-':'-'};
		bases = list(seq) ; sequences = '';
		for base in bases:
			comp_base = complement[base] 
			sequences += "".join(comp_base); 
		return sequences[::-1];

	################################################################### Function to determine sequence composition bias ####################################################################
	def comp_bias(self, infile1, infile2, infile3, outfile1, outfile2, outfile3):
		x=0; sequence='';
		(genus, data_dict) = self.data_extract(infile1, 3); 
		coodrinates, strand = data_dict['arr_0'], data_dict['arr_1']
		genestart = [coord.split('..')[0] for i, coord in enumerate(coodrinates)]
		geneend = [coord.split('..')[-1] for i, coord in enumerate(coodrinates)]; #print(genus, genestart, geneend);

		file2 = open (infile2,"r");
		for i in range(2):
			line2 = file2.readline()
		while line2:
			line2 = line2.replace('\n','')
			sequence += line2
			line2 = file2.readline()
		#print(len(sequence))

		shadow_gene,gene,shadow_genst,shadow_genen,genest,geneen=[],[],[],[],[],[];
		dinuc = {'AT':1,'AG':2,'AC':3,'AA':4,'TA':5,'TG':6,'TC':7,'TT':8,'GA':9,'GT':10,'GC':11,'GG':12,'CA':13,'CT':14,'CG':15,'CC':16}; #print(dinuc);
		dinuc_keys = dinuc.keys(); dinuc_values = dinuc.values(); #print(dinuc_keys,dinuc_values);
		trinuc = {'AAA':1,'AGA':2,'TAA':3,'TGA':4,'GAA':5,'GGA':6,'CAA':7,'CGA':8,'AAC':9,'AAG':10,'AAT':11,'ATA':12,'ATC':13,'ATG':14,'ATT':15,'AGC':16,'AGG':17,
'AGT':18,'ACA':19,'ACC':20,'ACG':21,'ACT':22,'TAC':23,'TAG':24,'TAT':25,'TTA':26,'TTC':27,'TTG':28,'TTT':29,'TGC':30,'TGG':31,'TGT':32,'TCA':33,'TCC':34,
'TCG':35,'TCT':36,'GAC':37,'GAG':38,'GAT':39,'GTA':40,'GTC':41,'GTG':42,'GTT':43,'GGC':44,'GGG':45,'GGT':46,'GCA':47,'GCC':48,'GCG':49,'GCT':50,
'CAC':51,'CAG':52,'CAT':53,'CTA':54,'CTC':55,'CTG':56,'CTT':57,'CGC':58,'CGG':59,'CGT':60,'CCA':61,'CCC':62,'CCG':63,'CCT':64}
		trinuc_keys = trinuc.keys(); trinuc_values = trinuc.values(); #print(trinuc_keys,trinuc_values);

		for i in range(len(strand)):
			if strand[x] == '-':
				sgen = sequence[int(genestart[x])-1:int(geneend[x])]; 
				shadow_gene.append(self.reverse_complement(sgen)); #print(genestart[x],geneend[x],"\n",shadow_gene,"\n");
				shadow_genst.append(genestart[x]); shadow_genen.append(geneend[x]);
			elif strand[x] == '+':
				gene.append(sequence[int(genestart[x])-1:int(geneend[x])]); #print(genestart[x],geneend[x],"\n",gene,"\n");
				genest.append(genestart[x]); geneen.append(geneend[x]);
			x += 1

		################################################################# Determining Dinucleotide Bias ################################################################################
		command = 'perl mlcafe.pl --phylogenetic --genus ' + genus + ' --fasta Chimeric.fna Chimeric.ptt --verbose'
		#os.system(command) ## deprecated
		args = shlex.split(command) ## To break the above command into correct set of arguments
		p = subprocess.Popen(args)  ## To run the command
		if p.poll() != None:  ## To check if the child process is completed
			pass
		else:
			p.wait(timeout=None)  ## To wait and complete the child process if it is not completed

		file6 = open (outfile1,"w"); file7 = open (outfile2,"w");

		for g in range(len(shadow_gene)):
			length_shadow = len(shadow_gene[g]); #print("length = ",length_shadow,shadow_genst[g],shadow_genen[g],"\n",shadow_gene[g]); 
			if (length_shadow%2 == 0):
				count_dinucs=[0]*(len(dinuc_keys)+1); 
				for j in range(0,int(length_shadow),2):
					sgenst = shadow_genst[g]; sgenen = shadow_genen[g]; sgene = shadow_gene[g]; 
					dint = sgene[j:j+2]; #print(shadow_gene[g],"\n",j,j+1,dint); 
					if str(dint) in dinuc_keys:
						val = dinuc.get(str(dint)); #print(str(dint),val);
						count_dinucs[int(val)] += 1; #print(count_dinucs[int(val)]);
				count_dinucs.pop(0); #print(count_dinucs);
				file6.write(shadow_genst[g]+"\t"+shadow_genen[g]+"\t"+str(count_dinucs)+"\n")
		
			elif (length_shadow%2 != 0):
				count_dinucs=[0]*(len(dinuc_keys)+1);  
				for j in range(0,int(length_shadow-1),2):
					sgenst = shadow_genst[g]; sgenen = shadow_genen[g]; sgene = shadow_gene[g];
					dint = sgene[j:j+2]; #print(shadow_gene[g],"\n",j,j+1,dint);
					if str(dint) in dinuc_keys:
						val = dinuc.get(str(dint)); #print(str(dint),val);
						count_dinucs[int(val)] += 1; #print(count_dinucs[int(val)])
				count_dinucs.pop(0); #print(count_dinucs);
				file6.write(shadow_genst[g]+"\t"+shadow_genen[g]+"\t"+str(count_dinucs)+"\n")

		for h in range(len(gene)):
			length_gene = len(gene[h]); #print("length = ",length_gene,genest[h],geneen[h]);
			if (length_gene%2 == 0):
				count_dinucs=[0]*(len(dinuc_keys)+1);
				for j in range(0,int(length_gene),2):
					genst = genest[h]; genen = geneen[h]; genes = gene[h]; 
					dint = genes[j:j+2]; #print(gene[h],"\n",j,j+1,dint);
					if str(dint) in dinuc_keys:
						val = dinuc.get(str(dint)); #print(str(dint),val);
						count_dinucs[int(val)] += 1; #print(count_dinucs[int(val)])
				count_dinucs.pop(0); #print(count_dinucs);
				file6.write(genest[h]+"\t"+geneen[h]+"\t"+str(count_dinucs)+"\n")

			elif (length_gene%2 != 0):
				count_dinucs=[0]*(len(dinuc_keys)+1);
				for j in range(0,int(length_gene-1),2):
					genst = genest[h]; genen = geneen[h]; genes = gene[h];
					dint = genes[j:j+2]; #print(gene[h],"\n",j,j+1,dint);
					if str(dint) in dinuc_keys:
						val = dinuc.get(str(dint)); #print(str(dint),val);
						count_dinucs[int(val)] += 1; #print(count_dinucs[int(val)])
				count_dinucs.pop(0); #print(count_dinucs);
				file6.write(genest[h]+"\t"+geneen[h]+"\t"+str(count_dinucs)+"\n")

		################################################################ Determining Codon Usage Bias #################################################################################
		for g in range(len(shadow_gene)):
			length_shadow = len(shadow_gene[g]); #print("length = ",length_shadow,shadow_genst[g],shadow_genen[g]); 
			count_trinucs=[0]*(len(trinuc_keys)+1);
			for j in range(0,int(length_shadow),3):
				sgenst = shadow_genst[g]; sgenen = shadow_genen[g]; sgene = shadow_gene[g]; 
				trint = sgene[j:j+3]; #print(shadow_gene[g],"\n",j,j+2,trint); 
				if str(trint) in trinuc_keys:
					vals = trinuc.get(str(trint)); #print(trint,trinuc.get(str(trint)));
					count_trinucs[int(vals)] += 1; #print(count_trinucs[int(vals)])
			count_trinucs.pop(0); #print(count_trinucs);
			file7.write(shadow_genst[g]+"\t"+shadow_genen[g]+"\t"+str(count_trinucs)+"\n")

		for h in range(len(gene)):
			length_gene = len(gene[h]); #print("length = ",length_gene,genest[h],geneen[h]);
			count_trinucs=[0]*(len(trinuc_keys)+1);
			for j in range(0,int(length_gene),3):
				genst = genest[h]; genen = geneen[h]; genes = gene[h]; 
				trint = genes[j:j+3]; #print(gene[h],"\n",j,j+2,trint);
				if str(trint) in trinuc:
					vals = trinuc.get(str(trint)); #print(trint,trinuc.get(str(trint)));
					count_trinucs[int(vals)] += 1; #print(count_trinucs[int(vals)])
			count_trinucs.pop(0); #print(count_trinucs);
			file7.write(genest[h]+"\t"+geneen[h]+"\t"+str(count_trinucs)+"\n")

		file6.close(); file7.close();

		###############################################################################################################################################################################
		dicts1, dicts2 = self.data_extract(outfile1, 0), self.data_extract(outfile2, 0) 
		gene_start, gene_end, dinucleotides = dicts1['arr_0'], dicts1['arr_1'], dicts1['arr_2']
		gene_strt, gene_ed, trinucleotides = dicts2['arr_0'], dicts2['arr_1'], dicts2['arr_2']

		dicts3 = self.data_extract(infile3, 0)
		gene_st, gene_en, gene_func, phylo_val, marker_val = dicts3['arr_1'], dicts3['arr_2'], dicts3['arr_3'], dicts3['arr_4'], dicts3['arr_5']

		file8 = open (outfile3, "w");
		for i in range(len(gene_st)):
			if ((gene_st[i] in gene_start) and (gene_en[i] in gene_end)):
				idval = gene_start.index(gene_st[i]); #print("index = ", idval); 
				file8.write(gene_start[idval]+"\t"+gene_end[idval]+"\t"+gene_func[i]+"\t"+phylo_val[i]+"\t"+marker_val[i]+"\t"+dinucleotides[idval]+"\t"+trinucleotides[idval]+"\n")

if __name__ == "__main__":
	seq_comp().comp_bias("Chimeric.ptt", "Chimeric.fna", "Output_Files/Chimeric_AllGenes.txt", "Output_Files/Seq_Dinuc.txt", "Output_Files/Seq_Codon.txt", "Output_Files/GenomesFeatures.txt")



