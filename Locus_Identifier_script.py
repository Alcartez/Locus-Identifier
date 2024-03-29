###################### LOCUS IDENTIFIER ######################

### PACKAGES ###

import os
import sys
import pandas as pd 
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
import csv
import subprocess

###### MAKE LOCAL DATABASE ######

db_fasta_raw_filename = sys.argv[2]
local_dbname = "local_db\local_alu_db"
print("\n Creating Local Blast Database ...")

makeblastdb_cline = NcbimakeblastdbCommandline(dbtype="nucl", input_file= db_fasta_raw_filename ,out = local_dbname )
print("\n CMD Command : " , makeblastdb_cline)
stdout,stderr = makeblastdb_cline()

###### RUN BLASTN ######

input_path = sys.argv[1]
input_file = os.path.basename(input_path)
seq_filename = input_file.split(".")[0]
seq_ext = input_file.split(".")[1]

# command_line = ['blastn','-query',
#                seq + '.fasta','-out',
#                seq + '_blout', '-outfmt',
#                '6','-db','local_alu_db']

command_line = "blastn -db local_db\local_alu_db -query " + seq_filename + "." + seq_ext + " " +"-out " + seq_filename + ".xml -outfmt 16"

# Sometimes it displays only 1 HSP even though there might be multiple instances. Will need to debug later.

print("\n CMD Command : " , command_line)
print("\n Result :")
subprocess.call(command_line , shell=True)

###### DEFINE BLAST PARSER ######

def parse_blast(resultfile): #takes in the BLAST result, outputs list that can be made into csv rows
    print("\n Running BLAST Parser ... \n")
    result_handle = open(resultfile)
    print("\n Opening " + resultfile + " ... \n")
    blast_records = NCBIXML.parse(result_handle)
    csv_list = []
    
    header = [  'Query',
                'Name', 'Length', 'Score', 'Expect',
                'QueryStart', 'QueryEnd',
                'SubjectStart', 'SubjectEnd',
                'QueryStrand' , 'HitStrand'
            ]
    
    csv_list.append(header)
    count = 0
    for blast_record in blast_records:
        '''help(blast_record.alignments[0].hsps[0])''' # these give help info for the parts 
        '''help(blast_record.alignments[0])        '''
        count +=1
        
        query = blast_record.query
        for alignment in blast_record.alignments:

            name = alignment.title
            length = alignment.length
            ### INCREASE HSP IF REQUIRED ###
            hsp = alignment.hsps[0] # I don't know if we will ever have more than one, so might as well take the first one.
            score = hsp.score
            expect = hsp.expect
            querystart = hsp.query_start
            queryend = hsp.query_end
            subjectstart = hsp.sbjct_start
            subjectend = hsp.sbjct_end
            querystrand = hsp.strand[0]
            hitstrand = hsp.strand[1]
            row = [query,name,length,score,expect,querystart,queryend,subjectstart,subjectend,querystrand,hitstrand]
            csv_list.append(row)
            
    result_handle.close()
    return csv_list

###### PARSED OUTPUT TO A CSV ######

blast_ot =  seq_filename + ".xml"
output_csv = parse_blast(blast_ot)
filename = "output.csv"
print("\n Creating Standard Blast Output in CSV Format .... \n")
with open(filename, 'w') as csvfile: 
    csvwriter = csv.writer(csvfile)
    for line in output_csv:
        csvwriter.writerow(line)

###### ALU LENGTH DATABASE ######

record_list = list(SeqIO.parse(db_fasta_raw_filename , "fasta"))
df = pd.DataFrame()
Seq_id_list = []
Seq_Length_list = []

for record in record_list :
    Seq_Length_list.append(len(record.seq))
    Seq_id_list.append(record.id)

df["Sequence Name"] = Seq_id_list
df["Sequence Length"] = Seq_Length_list

df.to_csv("Alu_Length_Database.csv")

#### Intact or Truncated ####
print("\n Scanning strands to see if they are intact or truncated ... \n ")
output_csv = pd.read_csv("output.csv")
ele_df = pd.read_csv("Alu_Length_Database.csv")

SS_list = []
SE_list = []


for i in range(0,len(output_csv['QueryStart'])) :

    ele_name = output_csv.at[i , "Name"].split(" ")[-1]

    if output_csv.at[i, "HitStrand"] == 'Plus' :

        if output_csv.at[i, "SubjectStart"] > 9 :
            SS_list.append('tr') # truncated
        else :
            SS_list.append('i') # intact

        if output_csv.at[i, "SubjectEnd"] < int(ele_df[ele_df["Sequence Name"]==ele_name]["Sequence Length"]) - 10 : # Original Length of element minus 10
            SE_list.append("tr")
        else:
            SE_list.append("i")

    elif output_csv.at[i, "HitStrand"] == 'Minus' :

        if output_csv.at[i, "SubjectEnd"] > 9 :
            SE_list.append('tr') # truncated
        else :
            SE_list.append('i') # intact

        if output_csv.at[i, "SubjectStart"] < int(ele_df[ele_df["Sequence Name"]==ele_name]["Sequence Length"]) - 10 :
            SS_list.append("tr")
        else:
            SS_list.append("i")

output_csv["Subject Quality 3` End"] = SS_list
output_csv["Subject Quality 5` End"] = SE_list

print("\n Creating output_2.csv ..... \n This file contain more scanned data ... \n")
output_csv.to_csv("output_2.csv")   

###### FINDING THE LOCUS ######

output_csv = pd.read_csv("output.csv")
Alu_start_position_list = output_csv.QueryStart
Alu_end_position_list = output_csv.QueryEnd
alu_locations_df = pd.DataFrame(list(zip(output_csv.Name, output_csv.QueryStart, output_csv.QueryEnd , output_csv.Score , output_csv.Length , output_csv.Expect , output_csv.QueryStrand , output_csv.HitStrand)),columns =['Name' ,'Start Position', 'End Position', 'Score' , 'Length' , 'E-Value','QueryStrand','HitStrand'])
seq_record = SeqIO.parse(open(str(seq_filename+'.fasta')),'fasta')


# Locus Finder
columns = ['Genome ID','Strand' ,'Score' ,'Start','End', 'Length' , 'E-Value','Upstream (<100)', 'Insertion Sequence (Query)','Downstream (<100)','Locus', 'Blast Query Match' ]
Locus_df = pd.DataFrame(columns)

for seq_object in seq_record :
    for i in range(0,alu_locations_df.shape[0]):
        seq = str(seq_object.seq)
        Query_Start = alu_locations_df.at[i, "Start Position"] - 1
        Query_End = alu_locations_df.at[i, "End Position"] - 1
        Query_Match = alu_locations_df.at[i, "Name"]
        Locus_Downstream = seq[max(0,Query_Start-100):Query_Start]
        Locus_Upstream = seq[Query_End + 1:min(Query_End + 101 , len(seq))]
        #print("\n Start Position : " , Query_Start , "\n End Position : " , Query_End ,"\n Match : " , Query_Match ,"\n Transposon : " , seq[Query_Start:Query_End] , "\n")
        Locus = Locus_Downstream + Locus_Upstream
        #print("\n Locus : " , Locus)
        if alu_locations_df.at[i,'HitStrand'] == 'Plus' :
            Locus_df = Locus_df.append({'Start': Query_Start,'End' : Query_End , 'Locus' : Locus, 'Blast Query Match'  : Query_Match , 'Upstream (<100)' : Locus_Upstream, 'Insertion Sequence (Query)': seq[Query_Start:Query_End],'Downstream (<100)' : Locus_Downstream , 'Genome ID' : seq_object.id , 'Strand' : "Positive (+ve)" ,'Score' : alu_locations_df.Score[i] , 'Length' : alu_locations_df.Length[i] , 'E-Value' : alu_locations_df['E-Value'][i] }, ignore_index = True)
        elif alu_locations_df.at[i,'HitStrand'] == 'Minus' : # If hit strand is minus , then the Upstream and Downstream swaps positions.
            Locus_df = Locus_df.append({'Start': Query_Start,'End' : Query_End , 'Locus' : Locus, 'Blast Query Match'  : Query_Match , 'Upstream (<100)' : Locus_Downstream, 'Insertion Sequence (Query)': seq[Query_Start:Query_End],'Downstream (<100)' : Locus_Upstream , 'Genome ID' : seq_object.id , 'Strand' : "Positive (+ve)" ,'Score' : alu_locations_df.Score[i] , 'Length' : alu_locations_df.Length[i] , 'E-Value' : alu_locations_df['E-Value'][i] }, ignore_index = True)
            


Locus_df.to_csv('LocusIdentifier_Output.csv', index = True)

# Upstream and Downstream

lf_data = pd.read_csv("LocusIdentifier_Output.csv")
Blast_match_idlist = lf_data['Blast Query Match']
Upstream_seqlist = lf_data['Upstream (<100)'] 
Downstream_seqlist = lf_data['Downstream (<100)']

#### UPSTREAM ####
print("\n Creating text file with Upstream sequences with custom sample IDs ... \n")
with open(seq_filename + '_upstream.txt', 'w') as f:
    for i in range(0,len(Upstream_seqlist)) :
        f.write('>SMPL' + str(i) + "_UP " + str(Blast_match_idlist[i]))
        f.write("\n" + str(Upstream_seqlist[i]) + "\n")
print("\n Done ...")

#### DOWNSTREAM ####
print("\n Creating text file with Downstream sequences with custom sample IDs ... \n")
with open(seq_filename + '_downstream.txt', 'w') as f:
    for i in range(0,len(Downstream_seqlist)) :
        f.write('>SMPL' + str(i) + "_DWN " + str(Blast_match_idlist[i]))
        f.write("\n" + str(Downstream_seqlist[i]) + "\n")
print("\n Done ...")

        