import re
from sequence_composition import seq_comp

class Average_SeqComp:
	def Avg_Composition(self, infile1, infile2, infile3, outfile):
		seq_comp().comp_bias("Chimeric.ptt", "Chimeric.fna", "Output_Files/Chimeric_AllGenes.txt", "Output_Files/Seq_Dinuc.txt", "Output_Files/Seq_Codon.txt", "Output_Files/GenomesFeatures.txt")

########################################################## Average Dinucleotide Count ######################################################
		count = [0]*16; list1, list2 = [],[];
		dict1, dict2 = seq_comp().data_extract(infile1, 0), seq_comp().data_extract(infile2, 0)
		genest, genend, dinucs = dict1['arr_0'], dict1['arr_1'], dict1['arr_2']
		dinuc = [re.sub('[^0-9|,]','',dint) for i, dint in enumerate(dinucs)]; #print(dinuc);
		split_dinuc = [din.split(",") for i, din in enumerate(dinuc)]; #print(split_dinuc);
		count_gene = len(genest)
		for l in range(len(split_dinuc)):
			for k in range(16):
				count[k] += int((split_dinuc[l])[k])
		for i in range(16):
			list1.append(count[i]/count_gene); #print(count[i]); 

############################################################# Average Codon Count ##########################################################
		counts = [0]*64
		genesstart, genesend, trinucs = dict2['arr_0'], dict2['arr_1'], dict2['arr_2']
		trinuc = [re.sub('[^0-9|,]','',trint) for i, trint in enumerate(trinucs)]; #print(trinuc);
		split_trinuc = [trin.split(",") for i, trin in enumerate(trinuc)]
		counts_genes = len(genesstart)
		for l in range(len(split_trinuc)):
			for m in range(64):
				counts[m] += int((split_trinuc[l])[m])
		for j in range(64):
			list2.append(counts[j]/counts_genes); #print(counts[j])
		#print('\n', list2)

############################################################################################################################################
		file3 = open(infile3, 'r'); file4 = open(outfile, 'w');
		line3 = file3.readline(); line4 = '';
		while line3:
			line3 = line3.replace('\n','')
			line3 = line3+"\t"+str(list1)+"\t"+str(list2)+"\n"
			line4 += line3
			line3 = file3.readline()
		file3.close(); #print(line4);

		file4.write(line4)
		file4.close()

if __name__ == "__main__":
	Average_SeqComp().Avg_Composition("Output_Files/Seq_Dinuc.txt", "Output_Files/Seq_Codon.txt", "Output_Files/GenomesFeatures.txt", "Output_Files/Genomes_Training_Features.txt")



