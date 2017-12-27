#!/usr/bin/env python

'''
Author: Anna Fei
Date: 21/March/2017
Usage: python gbk_to_fasta.py input.gk
Function: convert genbank file to fasta file. with header ">D1|Jamaica|Jamaica|DENV1-17388|2012|KU509249|10582", output is output_filter.fa
Version: python 2.7

'''


from Bio import SeqIO
import re
import datetime
import string
import sys
import os.path


def ConvertGB(input_handle,output_handle):
	#Short version:
	#SeqIO.write(SeqIO.parse(input_handle, "genbank"), output_handle, "fasta")

	#Long version, allows full control of fasta output
	for seq_record in SeqIO.parse(input_handle, "genbank") :
		if seq_record.annotations['data_file_division'] == 'VRL':
			for seq_feature in seq_record.features :
				if seq_feature.type=="source" :					
					try:
						country = seq_feature.qualifiers['country'][0]
						country=country.replace(" ", "_")
					except KeyError:
						country =''
						
				
					try:
						organism = seq_feature.qualifiers['organism'][0]
						if organism =="Dengue virus":
							Organism = ""
						elif organism =="Dengue virus 1":
							Organism = "DENV1"
						elif organism =="Dengue virus 2":
							Organism = "DENV2"
						elif organism =="Dengue virus 3":
							Organism = "DENV3"
						elif organism =="Dengue virus 4":
							Organism = "DENV4"
						else:
							Organism ="" 
					except KeyError:
						organism = ""
						Organism = ""
								

					try:
						collection_date = seq_feature.qualifiers['collection_date'][0]
						if "-" in collection_date:
							if re.search ('[a-zA-Z]', collection_date):
								try:
									collection_date= datetime.datetime.strptime(collection_date, '%d-%b-%Y').strftime('%Y%m%d')
								except:
									collection_date= datetime.datetime.strptime(collection_date, '%b-%Y').strftime('%Y%m')
							else:
								try:
									collection_date= datetime.datetime.strptime(collection_date, '%Y-%m-%d').strftime('%Y%m%d')
								except:
									collection_date= datetime.datetime.strptime(collection_date, '%Y-%m').strftime('%Y%m%d')
					except KeyError:
						collection_date = ""
						
					
					try:
						mol_type = seq_feature.qualifiers['mol_type'][0]
					except KeyError:
						mol_type =""
					

					try:
						host = seq_feature.qualifiers['host'][0]
						
						human_array = ["Homo sapiens","Homo Sapiens","Homo sapines","homo","homo sapiens","Homo sapien","Homo sapience","patient"]
						mosquito_array = ["Aedes","Anopheles","Armigeres","Culex","mosquito","mosquitoes"]
						if any (substring in host for substring in human_array):
							host ="homo"
						elif any (substring in host for substring in mosquito_array):
							host = "mosquito"
						else:
							host = host
					except KeyError:
						host = ""


					try:
						serotype = seq_feature.qualifiers['serotype']
						serotype1 = ['1','I','D1','D 1','Den-1','DENV-1','type 1','DEN-1','DENV1']
						serotype2 = ['2','II','D2','Den-2','DEN-2','DENV-2','DENV2',"Den'2'"]
						serotype3 = ['3','III','D3','dengue3','dengue 3','DENV3','DENV-3','Dengue 3']
						serotype4 = ['4','IV','D4','dengue4','DENV4','DENV-4','Dengue 4']
						if any (x in serotype1 for x in serotype):
							Serotype = 'DENV1'
						if any (x in serotype2 for x in serotype):
							Serotype = 'DENV2'
						if any (x in serotype3 for x in serotype):
							Serotype = 'DENV3'
						if any (x in serotype4 for x in serotype):
							Serotype = 'DENV4'					
					except KeyError:
						serotype = ""
						Serotype = ""


					try:
						strain = seq_feature.qualifiers['strain'][0]
						strain = strain.replace(" ", "_")
					except KeyError:
						strain =""					

							
					if Organism=="":
						if Serotype == "":
							break
						else:
							Virus_type=Serotype
					else:	
						Virus_type=Organism
									

				if seq_feature.type=="CDS" :			
					try:			
						product = seq_feature.qualifiers['product']
						product = [item.replace(" ", "_") for item in product]
					except KeyError:
						product = ""
					
					
		length = len(seq_record.seq)
		
		#if host=="":
		#	continue
		if collection_date=="":
			continue
		if country=="":
			continue
		
		output_handle.write(">%s|%s|%s|%s|%s|%s|%s|%s\n%s\n" % (
			Virus_type,
			country.split(":",1)[0],
			country,
			strain,
			collection_date,			
			seq_record.name,
			length,
			host,			
							
			seq_record.seq))
				

if __name__=="__main__":
	in_filePath = sys.argv[1]
	out_folder = sys.argv[2]

	path,in_fileName = os.path.split(in_filePath)
	FileBase = in_fileName.split(".",1)[0]
	out_filePath = out_folder + "/" + FileBase + "_filter.fa"

	input_handle  = open(in_filePath, "r")
	output_handle = open(out_filePath, "w")
	ConvertGB(input_handle,output_handle)
			
	output_handle.close()
	input_handle.close()
