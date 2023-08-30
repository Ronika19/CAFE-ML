import os
import re
from glob import glob
from sequence_composition import seq_comp
from Avg_dinuc_codon import Average_SeqComp

class GI_Data:
	def PosNeg_GI(self, infile1, infile2, outfile1, outfile2):
		################################################################## Extracting Positive and Negative GIs #######################################################################
		dict1, dict2 = seq_comp().data_extract(infile1, 1), seq_comp().data_extract(infile2, 0)
		pos_accession, pos_gistart, pos_giend = dict1['arr_0'], dict1['arr_1'], dict1['arr_2']
		neg_accession, neg_start, neg_end = dict2['arr_0'], dict2['arr_1'], dict2['arr_2']
		file2, file3 = open(outfile1, 'w'), open(outfile2, 'w')
		for i in range(len(pos_accession)):
			file2.write(pos_accession[i]+"\t"+pos_gistart[i]+"\t"+pos_giend[i]+"\n")

		for i in range(len(neg_accession)):
			file3.write(neg_accession[i]+"\t"+neg_start[i]+"\t"+neg_end[i]+"\n")

		file2.close(); file3.close();

	def Feature_Extraction(self, infile1, infile2, infile3, outfile1, outfile2, outfile3):
		Average_SeqComp().Avg_Composition("Output_Files/Seq_Dinuc.txt", "Output_Files/Seq_Codon.txt", "Output_Files/GenomesFeatures.txt", "Output_Files/Genomes_Training_Features.txt")
		self.PosNeg_GI("GI.txt", "NonGI.txt", "Output_Files/IPposGIs.txt", "Output_Files/IPnegGIs.txt")

		# Extracting training file from the bacterial chromosome folder(s)
		#path = glob('*.ptt'); split_path = (path[0]).split('/'); print(split_path[1]); # print(path); 
		split_path = 'Chimeric'+'1'
		c1, c2 =0,0; out1, out2=[],[];
		fl1 = open(outfile1, "w"); fl2 = open(outfile2, "w")

		dict3, dict4, dict5 = seq_comp().data_extract(infile1, 0), seq_comp().data_extract(infile2, 1), seq_comp().data_extract(infile3, 1)
		genest, genend, gene_func, phyloval, markerval, dinuc_value, trinuc_value, avg_dinuc, avg_trinuc = dict3['arr_0'], dict3['arr_1'], dict3['arr_2'], dict3['arr_3'], dict3['arr_4'], dict3['arr_5'], dict3['arr_6'], dict3['arr_7'], dict3['arr_8']

		###################################################################### Extracting GI Features #################################################################################
		accession_positive, gistart, giend = dict4['arr_0'], dict4['arr_1'], dict4['arr_2']
		accession = [re.sub('\.[0-9]', '', pos_accn) for i, pos_accn in enumerate(accession_positive)]; #print(split_path[-1]); ## Changes
		value, counter1, outcome1 = [],[],''
		for i in range(len(genest)):
			for j in range(len(accession)):
				if (accession[j] == re.sub('\.[0-9]','', split_path[:])): ## Changes
					if((int(genest[i])>=int(gistart[j])) and (int(genend[i])>int(gistart[j])) and (int(genest[i])<int(giend[j])) and (int(genend[i])<=int(giend[j]))):
						value.append("P"); #print(accession[j],gistart[j],giend[j],genest[i],genend[i],value[c1]);
						outs1=(str(accession[j])+"\t"+str(genest[i])+"\t"+str(genend[i])+"\t"+phyloval[i]+"\t"+markerval[i]+"\t"+dinuc_value[i]+"\t"+trinuc_value[i]+"\t"+avg_dinuc[i]+"\t"+avg_trinuc[i]+"\t"+"P"+"\n");
						out1.append(outs1);
						outcome1 += (accession[j]+"\t"+str(gistart[j])+"\t"+str(giend[j])+"\t"+str(genest[i])+"\t"+str(genend[i])+"\t"+value[c1]+"\n");
						c1 += 1; counter1.append(i); #print(outcome1)

		###################################################################### Extracting Non-GI Features ##############################################################################
		accession_negative, gi_start, gi_end = dict5['arr_0'], dict5['arr_1'], dict5['arr_2']
		accessions = [re.sub('\.[0-9]', '', neg_accn) for i, neg_accn in enumerate(accession_negative)]; ## Changes
		values, counter2, outcome2 = [],[],''
		for i in range(len(genest)):
			for j in range(len(accessions)):
				if accessions[j] == re.sub('\.[0-9]','', split_path[:]): ## Changes
					if((int(genest[i])>=int(gi_start[j])) and (int(genend[i])>int(gi_start[j])) and (int(genest[i])<int(gi_end[j])) and (int(genend[i])<=int(gi_end[j]))):
						values.append("N"); #print(accessions[j],gi_start[j],gi_end[j],genest[i],genend[i],values[c2]);
						outs2=(str(accessions[j])+"\t"+str(genest[i])+"\t"+str(genend[i])+"\t"+phyloval[i]+"\t"+markerval[i]+"\t"+dinuc_value[i]+"\t"+trinuc_value[i]+"\t"+avg_dinuc[i]+"\t"+avg_trinuc[i]+"\t"+"N"+"\n");
						out2.append(outs2); 
						outcome2 += (accessions[j]+"\t"+str(gi_start[j])+"\t"+str(gi_end[j])+"\t"+str(genest[i])+"\t"+str(genend[i])+"\t"+values[c2]+"\n");
						c2 += 1; counter2.append(i); #print(outcome2)

		outcome = outcome1+outcome2
		fl1.write(outcome)

		out = out1+out2
		finalout=[]; feature = '';
		for i in out:
			if i not in finalout:
				finalout.append(i)

		########################################################### Writing GI and Non-GI Features in Feature File #####################################################################
		q=0; iq=0; fl3 = open(outfile3, 'w')
		dict6 = seq_comp().data_wrangle(finalout, 0); 
		accession_no, genst, genn, phylo, markers, dint, trint, avg_dint, avg_trint, gival = dict6['arr_0'], dict6['arr_1'], dict6['arr_2'], dict6['arr_3'], dict6['arr_4'], dict6['arr_5'], dict6['arr_6'], dict6['arr_7'], dict6['arr_8'], dict6['arr_9']
		for i in range(len(finalout)):
			fl3.write(finalout[i])
			if (genst[q] in genst[:q] and genn[q] in genn[:q] and gival[q] == "M" and q != 0):
				genst.pop(q); genn.pop(q); accession_no.pop(q); phylo.pop(q); markers.pop(q); 
				dint.pop(q); trint.pop(q); gival.pop(q); avg_dint.pop(q); avg_trint.pop(q); #print(q-1); 
				q = q-1;
			if (genst[q] in genst[:q] and genn[q] in genn[:q] and ((gival[q] == "P") or (gival[q] == "N"))  and (q != 0)):
				st = genst.index(genst[q]); en = genn.index(genn[q]);  #print(st,q);
				genst.pop(st); genn.pop(st); accession_no.pop(st); phylo.pop(st); markers.pop(st); 
				dint.pop(q); trint.pop(q); gival.pop(st); avg_dint.pop(q); avg_trint.pop(q); 
				q = q-1;
			q += 1;	iq += 1;
		#print(q, iq);

		for m in range(q):
			#print(accession_no[m],genst[m],genn[m],phylo[m],markers[m],gival[m]);	
			feature += (str(accession_no[m])+"\t"+str(genst[m])+"\t"+str(genn[m])+"\t"+str(phylo[m])+"\t"+str(markers[m])+"\t"+str(dint[m])+"\t"+str(trint[m])+"\t"+str(avg_dint[m])+"\t"+str(avg_trint[m])+"\t"+str(gival[m])+"\n")
		fl2.write(feature);

		fl1.close(); fl2.close(); fl3.close();

if __name__ == "__main__":
	GI_Data().Feature_Extraction('Output_Files/Genomes_Training_Features.txt', 'Output_Files/IPposGIs.txt', 'Output_Files/IPnegGIs.txt', 'Output_Files/PosNeg.txt', 'Output_Files/GenesFeatures.txt', 'Output_Files/Features.txt')



