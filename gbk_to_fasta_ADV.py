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
						
				
					type_sub=["Human adenovirus 3p","Human adenovirus 4a","Human adenovirus 7a","Human adenovirus 7d","Human adenovirus 7d2","Human adenovirus 7h","Human adenovirus 7l","Human adenovirus 7m","Human adenovirus 7p","Human adenovirus 8E","Human adenovirus 11a","Human adenovirus 14p","Human adenovirus 14p1","Human adenovirus 16p","Human adenovirus 19a","Human adenovirus 19C","Human adenovirus 19p","Human adenovirus 21a","Human adenovirus 60a"]
					
					try:
						organism = seq_feature.qualifiers['organism'][0]
						if "Human" in organism or "human" in organism:
							if "adenovirus" in organism:
								Organism="Others"
								for i in range(1,85):
									if organism =="Human adenovirus "+ str(i):
										Organism = "ADV" + str(i)
								for ii in type_sub:
									if organism == ii:
										type_sub = organism.split(" ",2)[2]
										Organism = "ADV" + re.search(r'\d+',type_sub).group()
							else:
								Organism ="Others"
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
						if any (substring in host for substring in human_array):
							host ="homo"
						else:
							host = ""
					except KeyError:
						host = ""

				
					try:
						serotype = seq_feature.qualifiers['serotype']
						serotype3 = ['3','Human adenovirus 3','3a variant','3a','3a2','Type 3','HAdV-3','HAdV 3']
						serotype4 = ['4','human adenovirus 4','type 4','HAdV4','HAdV-4','HAdV 4']
						serotype5 = ['5','Human adenovirus type 5','H5N1']
						serotype6 = ['6']
						serotype7 = ['7','HAdV-7']
						if any (x in serotype3 for x in serotype):
							Serotype = 'ADV3'
						elif any (x in serotype4 for x in serotype):
							Serotype = 'ADV4'
						elif any (x in serotype5 for x in serotype):
							Serotype = 'ADV5'
						elif any (x in serotype6 for x in serotype):
							Serotype = 'ADV6'
						elif any (x in serotype7 for x in serotype):
							Serotype = 'ADV7'
						else:
							Serotype ="Others"				
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
		
		if host=="":
			continue
		if collection_date=="":
			continue
		if country=="":
			continue
		if Virus_type=="Others":
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
