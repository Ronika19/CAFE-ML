import re
import os
from sequence_composition import seq_comp
from post_Training import GI_Data

class DataProcess:
	def data_preprocess(self, infile, outfile1, outfile2):
		GI_Data().Feature_Extraction('Output_Files/Genomes_Training_Features.txt', 'Output_Files/IPposGIs.txt', 'Output_Files/IPnegGIs.txt', 'Output_Files/PosNeg.txt', 'Output_Files/GenesFeatures.txt', 'Output_Files/Features.txt')

		fl1 = open(outfile1, 'w'); fl2 = open(outfile2, 'w');

		dicts1 = seq_comp().data_extract(infile, 0)
		accession, genestrt, gened, phylo, marker, dint, trint, avg_dint, avg_trint, label = dicts1['arr_0'], dicts1['arr_1'], dicts1['arr_2'], dicts1['arr_3'], dicts1['arr_4'], dicts1['arr_5'], dicts1['arr_6'], dicts1['arr_7'], dicts1['arr_8'], dicts1['arr_9']
		genest, genend = [int(val) for i, val in enumerate (genestrt)], [int(value) for i, value in enumerate (gened)]	
		
		genestart = [];
		for items in genest:
			genestart.append(int(items))
		genestart.sort();
		for i in range(len(genest)):
			if (genestart[i] in genest):
				j = genest.index(genestart[i])
				dinuc, trinuc = re.sub(',', '\t', re.sub('\[|\]', '', str(dint[j]))), re.sub(',', '\t', re.sub('\[|\]', '', str(trint[j])))
				avg_dinuc, avg_trinuc = re.sub(',', '\t', re.sub('\[|\]', '', str(avg_dint[j]))), re.sub(',', '\t', re.sub('\[|\]', '', str(avg_trint[j])))
				# If you want to add avg_dint[j] & avg_trint[j], use this file.
				fl1.write(str(accession[j])+"\t"+str(genest[j])+"\t"+str(genend[j])+"\t"+str(phylo[j])+"\t"+str(marker[j])+"\t"+str(dinuc)+"\t"+str(trinuc)+"\t"+str(avg_dinuc)+"\t"+str(avg_trinuc)+"\t"+str(label[j])+"\n")

				# If you do not want to add avg_dint[j] & avg_trint[j], use this file.
				fl2.write(str(accession[j])+"\t"+str(genest[j])+"\t"+str(genend[j])+"\t"+str(phylo[j])+"\t"+str(marker[j])+"\t"+str(dinuc)+"\t"+str(trinuc)+"\t"+str(label[j])+"\n")
		fl1.close(); fl2.close();
		
if __name__ == "__main__":
	DataProcess().data_preprocess('Output_Files/GenesFeatures.txt', 'Output_Files/FeaturesProcess.txt', 'Output_Files/GeneFeat.tsv')




