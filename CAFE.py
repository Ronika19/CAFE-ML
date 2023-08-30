from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB, MultinomialNB
from sklearn.svm import SVC
from sklearn.decomposition import PCA, KernelPCA
from sklearn.inspection import permutation_importance
from sklearn.metrics import confusion_matrix, recall_score, precision_score, f1_score
from sklearn.metrics import roc_curve, precision_recall_curve, roc_auc_score
from imblearn.over_sampling import RandomOverSampler, SMOTE ### Need to install imbalanced-learn
from imblearn.under_sampling import RandomUnderSampler
from collections import Counter
from matplotlib import pyplot
import seaborn as sns
import pandas as pd
import numpy as np
import argparse
import re
import os
import sys
import glob

from sequence_composition import seq_comp
from data_process import DataProcess

class GI_Genome_Prediction:
	def FeatureSelection(self, feature_model, clf, clf_key, X_train, y_train):
		#if ((clf_key == 'RF-EST10') or (clf_key == 'ABC') or (clf_key == 'DT') or (clf_key == 'ETC')):
		#feature_importance = (clf.feature_importances_).tolist()
		feature_importance = (feature_model.feature_importances_).tolist()
		feature_importance_sort = sorted(feature_importance, reverse=True); #print(feature_importance);
		top_column_indices = []
		for i in range(len(feature_importance_sort)):
			if feature_importance_sort[i] in feature_importance:
				indexes = int(feature_importance.index(feature_importance_sort[i]))
				top_column_indices.append(indexes)
		return top_column_indices

	def gi_predict(self):
		DataProcess().data_preprocess('Output_Files/GenesFeatures.txt', 'Output_Files/FeaturesProcess.txt', 'Output_Files/GeneFeat.tsv')

		# Initialize parser
		parser = argparse.ArgumentParser()
 		# Adding optional argument
		parser.add_argument('-i', '--Input', required = True, type=str)
		parser.add_argument('-ptt', '--Annot', required = True, type=str)
		parser.add_argument("-s", "--Parameter", type=str) ### PCA or kPCA
		parser.add_argument("-o", "--Output", required = True, type=str)
		
		# Read arguments from command line
		args = parser.parse_args(); #print(args); print(args.Input, args.Output, args.Annot)

		inputf = pd.read_csv(args.Input, sep='\t')
		dataset = inputf.iloc[:, 4:-1].values; #print(eval_dataset[:5], '\n\n\n'); #eval_dataset = X_imb_data
				
		data = pd.read_csv('concatenated_Features.tsv', sep='\t', header=None); #print(len(data.columns), data.head());
		dinuc = {'AT':1,'AG':2,'AC':3,'AA':4,'TA':5,'TG':6,'TC':7,'TT':8,'GA':9,'GT':10,'GC':11,'GG':12,'CA':13,'CT':14,'CG':15,'CC':16}; #print(dinuc);
		dinuc_keys = list(dinuc.keys()); dinuc_values = list(dinuc.values()); #print(dinuc_keys,dinuc_values);
		trinuc = {'AAA':1,'AGA':2,'TAA':3,'TGA':4,'GAA':5,'GGA':6,'CAA':7,'CGA':8,'AAC':9,'AAG':10,'AAT':11,'ATA':12,'ATC':13,'ATG':14,'ATT':15,'AGC':16,'AGG':17,
'AGT':18,'ACA':19,'ACC':20,'ACG':21,'ACT':22,'TAC':23,'TAG':24,'TAT':25,'TTA':26,'TTC':27,'TTG':28,'TTT':29,'TGC':30,'TGG':31,'TGT':32,'TCA':33,'TCC':34,
'TCG':35,'TCT':36,'GAC':37,'GAG':38,'GAT':39,'GTA':40,'GTC':41,'GTG':42,'GTT':43,'GGC':44,'GGG':45,'GGT':46,'GCA':47,'GCC':48,'GCG':49,'GCT':50,
'CAC':51,'CAG':52,'CAT':53,'CTA':54,'CTC':55,'CTG':56,'CTT':57,'CGC':58,'CGG':59,'CGT':60,'CCA':61,'CCC':62,'CCG':63,'CCT':64}
		trinuc_keys = list(trinuc.keys()); trinuc_values = list(trinuc.values()); #print(trinuc_keys,trinuc_values);
		columnn_names = ['Accession', 'GeneStart', 'GeneEnd', 'AberrantPhyleticPattern', 'Marker', 'Label']

		df_cols = []
		for i in range(5):
			df_cols.append(columnn_names[i])

		for i in range(len(dinuc_keys)):
			df_cols.append(dinuc_keys[i])
			
		for i in range(len(trinuc_keys)):
			df_cols.append(trinuc_keys[i])

		df_cols.append(columnn_names[5])

		data.columns = df_cols; #print(len(data.columns)); print(data.head())		

		X = data.iloc[:, 4:-1].values;
		y = data.iloc[:, -1].values; 
		for i, value in enumerate(y):
				if y[i] == 'P': 
					y[i] = 1
				elif y[i] == 'N':
					y[i] = 0
		
		y=y.astype('int'); #print(y); print(X);

		acc_datafm = data.iloc[:,0]; genest_datafm = data.iloc[:,1]; genend_datafm = data.iloc[:,2];

		X_new = [value for i, value in enumerate (X)]; #print(type(X_new), type(acc_datafm));
		X_new_data, acc_data, genest_data, genend_data = [],[],[],[]
		for i in range(len(X_new)):
			X_dataset = []
			for j in range(len(X_new[i])):
				X_dataset.append((X_new[i])[j])
			X_array = (np.array(X_dataset)).tolist(); 
			X_new_data.append(X_array)
		
		print(f"The training dataset has {sorted(Counter(y).items())[0][1]} records for the majority class and {sorted(Counter(y).items())[1][1]} records for the minority class.")
		
		################################################################### Solving Dataset Imbalance ##################################################################################
		# Randomly up/over sample the minority class
		ros = RandomOverSampler(random_state=50)
		X_ros, y_ros= ros.fit_resample(X, y); print(sorted(Counter(y_ros).items())) # Check the number of records after over sampling		

		# SMOTE oversampling of minority class
		X_smote, y_smote = SMOTE().fit_resample(X, y); print(sorted(Counter(y_smote).items()))
		
		imb_dataset = {'ROS':X_ros}; 
		imb_class = {'ROS':y_ros}
		
		for imb, imb_data in imb_dataset.items():
			for t in range(len(imb_data)):
				if ((imb_data[t]).tolist() in X_new_data):
					data_indices = int(X_new_data.index((imb_data[t]).tolist())); #print(data_indices)
					acc_data.append(acc_datafm[data_indices]); genest_data.append(genest_datafm[data_indices]); genend_data.append(genend_datafm[data_indices]);

		######################################################## Training ML models for GI gene prediction & K-fold Cross Validation #################################################
		# RandomForestClassifier(n_estimators=100)
		classifiers_dict = {'RF-EST10':RandomForestClassifier(n_estimators=10), 'LDA':LinearDiscriminantAnalysis(), 'LR-LFGS':LogisticRegression(solver='lbfgs'), 'LR-LIBLINEAR':LogisticRegression(solver='liblinear'), 'ETC':ExtraTreesClassifier(n_estimators=10), 'ABC':AdaBoostClassifier(n_estimators=10), 'DT':DecisionTreeClassifier(max_depth=None, min_samples_split=2), 'GB':GaussianNB(), 'SVC-RBF':SVC(kernel='rbf', probability=True)}
		feature_imp_dict, feature_imp_sorted_dict = {},{}
		for X_imb, X_imb_data in imb_dataset.items():
			for key, classifier in classifiers_dict.items():
				# KFold Cross Validation approach
				kf = KFold(n_splits=10, shuffle=True)
				kf_split = kf.split(X_imb_data) 

				total_recall, total_precision, total_f1_score, counter = 0,0,0,0; X_impfeat_datafm = pd.DataFrame()
				test_indices, test_labels_y, test_labels_ypred = [],[],[]
				recall_list, precision_list, f1score_list = [],[],[]

				# Iterate over each train-test split
				for train_index, test_index in kf_split:
					# Split train-test
					X_train, X_test = X_imb_data[train_index], X_imb_data[test_index]
					y_train, y_test = (imb_class[X_imb])[train_index], (imb_class[X_imb])[test_index]

					# Feature Selection
					# Fit the model
					feature_model = RandomForestClassifier(n_estimators=10)
					model = feature_model.fit(X_train, y_train)
					feature_imp_column_indices = self.FeatureSelection(feature_model, classifier, key, X_train, y_train); print('Feature Column Indices : ', feature_imp_column_indices);
					df = pd.DataFrame(X_imb_data); 
					if (feature_imp_column_indices is None):
						X_impfeat_datafm = df.iloc[:, :]
					else:
						X_impfeat_datafm = df.iloc[:, feature_imp_column_indices[:25]]; print('Feature Column Indices : ', feature_imp_column_indices);
					X_impfeat_array = X_impfeat_datafm.to_numpy(); #print(X_impfeat_array);
					X_impfeat_train, X_impfeat_test = X_impfeat_array[train_index], X_impfeat_array[test_index]
					dataset_df = pd.DataFrame(dataset) 
					sample_impfeat = (dataset_df.iloc[:, feature_imp_column_indices[:25]]).to_numpy() 

					# Train the model
					model = classifier.fit(X_impfeat_train, y_train)

					# Testing the model on the selected important Features	
					y_pred_test = classifier.predict(X_impfeat_test)
					
					# Store the test indices
					test_indices.append(test_index)
					test_labels_y.append((imb_class[X_imb])[test_index]); #print(test_labels_y);
					test_labels_ypred.append(y_pred_test)

					# Accuracy determination
					tn, fp, fn, tp = confusion_matrix(y_test, y_pred_test, labels=[0, 1]).ravel(); 
					#test_positives = np.sum(y_test); #print('TN = ', tn, 'FP = ', fp, 'FN = ', fn, 'TP = ', tp, 'True Positives = ', test_positives);
					rf_testprecision, rf_testrecall = precision_score(y_test, y_pred_test), recall_score(y_test, y_pred_test)
					rf_testf1 = f1_score(y_test, y_pred_test); 
					print("Recall_"+str(key)+'_'+str(X_imb)+" = ", rf_testrecall, "Precision_"+str(key)+'_'+str(X_imb)+" = ", rf_testprecision, "F1 Score_"+str(key)+'_'+str(X_imb)+" = ", rf_testf1)
					counter += 1
					total_recall += rf_testrecall
					total_precision += rf_testprecision
					total_f1_score += rf_testf1
	
					recall_list.append(float(rf_testrecall))
					precision_list.append(float(rf_testprecision))
					f1score_list.append(float(rf_testf1))

					rf_kfprecision, rf_kfrecall, _  = precision_recall_curve(y_test,y_pred_test)
					# Precision-Recall curve
					no_skill = len(y_train[y_train==1]) / len(y_train)
					pyplot.plot([0, 1], [no_skill, no_skill], linestyle='--', label='No Skill')
					pyplot.plot(rf_kfrecall, rf_kfprecision, marker='.', label=str(key))
					# axis labels
					pyplot.ylabel('Recall')
					pyplot.xlabel('Precision')
					# show the legend
					pyplot.legend()
					# show the plot
					#pyplot.show()
		
				########################################################## Predicting genes harbored by GIs (GI genes) ########################################################
				eval_dataset = sample_impfeat
				acc_eval = inputf.iloc[:, 0]; genest_eval = inputf.iloc[:, 1]; genend_eval = inputf.iloc[:, 2]; 
				prediction = classifier.predict(eval_dataset); #print(prediction)
				pos_gi_genestart = []; outfl = open('ML_Output_Files/'+str(key)+'_'+str(X_imb)+'_Prediction.txt', 'w')
				for p in range(len(prediction)):
					if (int(prediction[p]) == 1):
						pos_gi_genestart.append(str(genest_eval[p])); #print(prediction[p], genest_eval[p]);
						outfl.write(str(acc_eval[p])+'\t'+str(genest_eval[p])+'\t'+str(genend_eval[p])+'\t'+str(prediction[p])+'\n')
				#print(pos_gi_genestart)

				ptt_file = open(args.Annot, 'r')
				lines = ptt_file.readlines()
				genome_size = (((((lines[0]).rstrip()).split())[-1]).split('..'))[-1]; #print(genome_size)
				genestart, geneend = [],[]
				for ln in lines[3:]:
					split_ln = ((ln.split('\t'))[0]).split('..')
					genestart.append(split_ln[0])
					geneend.append(split_ln[1])

				################################################################### Determining putative GIs ##################################################################
				outfile = open('ML_Output_Files/'+str(key)+'_'+str(X_imb)+'_All_GIs.txt','w')
				outfile.write('Accession'+'\t'+'GI_Start'+'\t'+'GI_End'+'\t'+'GI_Gene_Count'+'\t'+'Total_Gene_Count'+'\n')
				gstart, gend = 0,0;
				for x in range(len(prediction)):
					count_genes = 0; count_gi_genes = 0; #print(acc_eval[x], genest_eval[x], genend_eval[x], prediction[x])
		
					# If the gene is a GI gene, as predicted by Random Forest model
					if (int(prediction[x]) == 1): 
						for j in range(len(genestart)):
				
							# Check for a particular genome
							# Check to see if the gene has 4000 kb upstream and downstream
							if ((int(genest_eval[x]) >= 4000) and ((int(genend_eval[x]) >= 8000) or (int(genend_eval[x]) <= 8000))): 
								if (int(genome_size)-int(genend_eval[x]) >= 4000):
									gstart = int(genest_eval[x])-4000
									gend = int(genend_eval[x])+4000
									if ((gstart <= int(genestart[j])) and (gstart < int(geneend[j])) and (gend > int(genestart[j])) and (gend >= int(geneend[j]))):
										count_genes += 1; #print(genestart[j]) # Total genes in the 8000 kb region of predicted GI gene
										if (str(genestart[j]) in pos_gi_genestart): # If the gene is a GI gene
											count_gi_genes += 1; #print(genestart[j])

								if (int(genome_size)-int(genend_eval[x]) < 4000):
									gend = int(genome_size)
									gstart = int(genest_eval[x])-(8000-(gend-int(genend_eval[x])))
									if ((gstart <= int(genestart[j])) and (gstart < int(geneend[j])) and (gend > int(genestart[j])) and (gend >= int(geneend[j]))):
										count_genes += 1; # Total genes in the 8000 kb region of predicted GI gene
										if (str(genestart[j]) in pos_gi_genestart): # If the gene is a GI gene
											count_gi_genes += 1; 

							# If the gene has less than 4000 kb upstream and has 4000 kb downstream
							if ((int(genest_eval[x]) <= 4000) and ((int(genend_eval[x]) <= 4000) or (int(genend_eval[x]) >= 4000))):
								gstart = 1
								gend = int(genend_eval[x])+(8000-(int(genest_eval[x])-gstart))
								if ((gstart <= int(genestart[j])) and (gstart < int(geneend[j])) and (gend > int(genestart[j])) and (gend >= int(geneend[j]))):
									count_genes += 1; # Total genes in the 8000 kb region of predicted GI gene
									if (str(genestart[j]) in pos_gi_genestart): # If the gene is a GI gene
										count_gi_genes += 1; 

						outfile.write(str(acc_eval[x])+'\t'+str(gstart)+'\t'+str(gend)+'\t'+str(count_gi_genes)+'\t'+str(count_genes)+'\n')
						gstart, gend = 0,0;

				outfile.close(); outfl.close()
			
				################################################################ Merging contiguous GI segments ###############################################################
				dict1 = seq_comp().data_extract('ML_Output_Files/'+str(key)+'_'+str(X_imb)+'_All_GIs.txt', 1)
				accessn_list, gi_start_list, gi_end_list, gi_gene_count_list, total_gene_count_list = dict1['arr_0'], dict1['arr_1'], dict1['arr_2'], dict1['arr_3'], dict1['arr_4'] 
				accessn, gi_start, gi_end, gi_gene_count, total_gene_count = [],[],[],[],[];
				for i in range(len(accessn_list)):
					if ((float(gi_gene_count_list[i]) > 0) and (float(total_gene_count_list[i]) > 0)):
						gi_gene_percent = float(gi_gene_count_list[i])/float(total_gene_count_list[i])
					if (float(gi_gene_percent) > 0.5): # try 0.9 also, it gives lower recall but very high f-measure and precision
						accessn.append(accessn_list[i])
						gi_start.append(gi_start_list[i])
						gi_end.append(gi_end_list[i])
						gi_gene_count.append(gi_gene_count_list[i])
						total_gene_count.append(total_gene_count_list[i])

				outfile1 = open('ML_Output_Files/'+str(key)+'_'+str(X_imb)+'_Merged_GIs.txt', 'w')
				accession1, gi_start1, gi_end1 = [],[],[];
				for x in range(len(gi_start)):
					try:
						if (accessn[x+1] == accessn[x]):
							if (int(gi_start[x+1]) <= int(gi_end[x])):
								gi_start[x] = gi_start[x]
								gi_end[x] = gi_end[x+1]; #print(accessn[x], gi_start[x], gi_end[x], gi_start[x+1], gi_end[x+1])
								gi_start[x+1] = gi_start[x] 
								gi_end[x+1] = gi_end[x+1]; #print(accessn[x], gi_start[x], gi_end[x+1])
								accession1.append(accessn[x])
								gi_start1.append(gi_start[x])
								gi_end1.append(gi_end[x+1])
								outfile.write(accessn[x]+'\t'+gi_start[x]+'\t'+gi_end[x+1]+'\n')
							if (int(gi_start[x+1]) > int(gi_end[x])):
								accession1.append(accessn[x])
								gi_start1.append(gi_start[x])
								gi_end1.append(gi_end[x])
								outfile1.write(accessn[x]+'\t'+gi_start[x]+'\t'+gi_end[x]+'\n'); #print(accessn[x], gi_start[x], gi_end[x])
					except:
						continue

				if (len(accessn) > 0):
					accession1.append(accessn[-1])
					gi_start1.append(gi_start[-1])
					gi_end1.append(gi_end[-1])
					outfile1.write(accessn[-1]+'\t'+gi_start[-1]+'\t'+gi_end[-1]+'\n')

					gi_start2, accession2, gi_end2 = [],[],[]
					for i in range(len(gi_start1)):
						if (gi_start1[i] not in gi_start2):
							gi_start2.append(gi_start1[i])
							gi_end2.append(gi_end1[i])
							accession2.append(accession1[i])

					accession3, gis_start3, gi_end3, gis_end = [],[],[],[]
					outfile2 = open('ML_Output_Files/'+str(key)+'_'+str(X_imb)+'_GIs.txt','w')
					for i in range(len(gi_start2)):
						for j in range(len(gi_start1)):
							if (accession1[j] == accession2[i]):
								if (int(gi_start1[j]) == int(gi_start2[i])):
									gis_end.append(int(gi_end1[j]))
						accession3.append(accession2[i])
						gis_start3.append(int(gi_start2[i]))
						gi_end3.append(int(max(gis_end))); print(accession2[i],gi_start2[i],max(gis_end));
						outfile2.write(accession2[i]+'\t'+str(gi_start2[i])+'\t'+str(max(gis_end))+'\n')
						gis_end = []
		
GI_Genome_Prediction().gi_predict()



