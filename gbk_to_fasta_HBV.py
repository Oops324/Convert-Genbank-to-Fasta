#!/usr/bin/env python

'''
Author: Anna Fei
Date: 21/March/2017
Usage: python gbk_to_fasta_HBV.py input.gk
Function: convert genbank file to fasta file. with header ">D1|Jamaica|Jamaica|DENV1-17388|2012|KU509249|10582", output_filter.fa
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

					hbv_array=["HBV","Hepatitis B virus"]
					try:
						organism = seq_feature.qualifiers['organism'][0]
						if any (substring in organism for substring in hbv_array):
							Organism = "HBV"
						else:
							Organism =""
					except KeyError:
						organism = ""
						Organism = ""
				
				
					try:
						collection_date = seq_feature.qualifiers['collection_date'][0]
						if "/" in collection_date:
							collection_date = collection_date.split('/',1)[0]			#2015-06/2015-08
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
						if any (substring in host for substring in human_array):
							host ="homo"
						else:
							host =""
					except KeyError:
						host = ""


					try:
						serotype = seq_feature.qualifiers['serotype'][0]
						Serotype = serotype[:1]					
					except KeyError:
						serotype = ""
						Serotype = ""


					try:
						note = seq_feature.qualifiers['note'][0]
						type_array = ["genotype","subgenotype","subtype","subgenotpe"]
					
						if any (substring in note for substring in type_array):
							segment = note.count(";")
							subtype_array = ["subgenotype", "subtype","subgenotpe"]
							if "not genotype" in note:
								genotype = ""
							elif "recombinant" in note:
								genotype = "recom"	
							elif segment ==0:	
								genotype = note.rsplit(None,1)[-1][:1].strip()
								if 6<=len(note.split()):
									genotype = ""
							elif segment==1:
								element1 = note.split(";",1)[0].strip()
								element2 = note.split(";",1)[1].strip()

								if element1.startswith("genotype"):
									genotype = element1.rsplit(None,1)[-1][:1]
								elif element2.startswith("genotype"):
									genotype = element2.rsplit(None,1)[-1][:1]	
								else:
									genotype =""
							elif segment ==2:
								element1 = note.split(";",2)[0].strip() 
								element2 = note.split(";",2)[1].strip()
								element3 = note.split(";",2)[2].strip()

								if element1.startswith("genotype"):
									genotype=element1.rsplit(None,1)[-1][:1]	
								elif element2.startswith("genotype"):
									genotype=element2.rsplit(None,1)[-1][:1]
								elif element3.startswith("genotype"):
									genotype=element3.rsplit(None,1)[-1][:1]
								else:
									genotype =""	
					except KeyError:
						genotype = ""
				


					try:
						strain = seq_feature.qualifiers['strain'][0]
						strain = strain.replace(" ", "_")
					except KeyError:
						strain =""
					
	
					if genotype=="":
						if Serotype =="":
							if Organism == "":
								break
							else:
								Virus_type=Organism
								Virus_type="UNDEF"
						else:
							Virus_type=Serotype.upper()
					else:	
						Virus_type=genotype.upper()
									

			
				if seq_feature.type=="CDS" :			
					try:	
						product = seq_feature.qualifiers['product']
						product = [item.replace(" ", "_") for item in product]
					except KeyError:
						product = ""
											
			length = len(seq_record.seq)
			
			if host=="":
				continue
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

