###################### LOCUS IDENTIFIER ######################

###### MAKE LOCAL DATABASE ######
from Bio.Blast.Applications import NcbimakeblastdbCommandline
db_fasta_raw_filename = "13396_Alu_Elements.fsa"
local_dbname = "local_alu_db"
print("\n Creating Local Database ...")

makeblastdb_cline = NcbimakeblastdbCommandline(dbtype="nucl", input_file= db_fasta_raw_filename ,out = local_dbname )
print("\n CMD Command : " , makeblastdb_cline)
stdout,stderr = makeblastdb_cline()

###### RUN BLASTN ######

from Bio.Blast.Applications import NcbiblastnCommandline
import subprocess

seq = "blast_input" 

print("0 = pairwise")
print("1 = query-anchored showing identities")
print("2 = query-anchored no identities")
print("3 = flat query-anchored, show identities")
print("4 = flat query-anchored, no identities")
print("5 = XML Blast output")
print("6 = tabular")
print("7 = tabular with comment lines")
print("8 = Text ASN.1")
print("9 = Binary ASN.1")
print("10 = Comma-separated values")

# command_line = ['blastn','-query',
#                seq + '.fasta','-out',
#                seq + '_blout', '-outfmt',
#                '6','-db','local_alu_db']

command_line = "blastn -db local_alu_db -query " + seq + ".fasta " +"-out " + seq + ".xml -outfmt 5 "

# Sometimes it displays only 1 HSP even though there might be multiple instances. Will need to debug later.

print("\n CMD Command : " , command_line)
print("\n Result :")
subprocess.call(command_line)

###### DEFINE BLAST PARSER ######

def parse_blast(resultfile): #takes in the BLAST result, outputs list that can be made into csv rows
    from Bio.Blast import NCBIXML
    result_handle = open(resultfile)
    blast_records = NCBIXML.parse(result_handle)
    csv_list = []
    
    header = [  'Query',
                'Name', 'Length', 'Score', 'Expect',
                'QueryStart', 'QueryEnd',
                'SubjectStart', 'SubjectEnd'
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
    
            hsp = alignment.hsps[0] # I don't know if we will ever have more than one, so might as well take the first one.
            score = hsp.score
            expect = hsp.expect
            querystart = hsp.query_start
            queryend = hsp.query_end
            subjectstart = hsp.sbjct_start
            subjectend = hsp.sbjct_end
            row = [query,name,length,score,expect,querystart,queryend,subjectstart,subjectend]
            csv_list.append(row)
            
    result_handle.close()
    return csv_list

###### PARSED OUTPUT TO A CSV ######

import pandas as pd
import csv
blast_ot =  "blast_input.xml"
output_csv = parse_blast(blast_ot)
filename = "output.csv"
print("\n Creating Blast Output in CSV Format ....")
with open(filename, 'w') as csvfile: 
    csvwriter = csv.writer(csvfile)
    for line in output_csv:
        print(line)
        csvwriter.writerow(line)

###### FINDING THE LOCUS ######

import pandas as pd 
from Bio import SeqIO

output_csv = pd.read_csv("output.csv")
Alu_start_position_list = output_csv.QueryStart
Alu_end_position_list = output_csv.QueryEnd
alu_locations_df = pd.DataFrame(list(zip(output_csv.Name, output_csv.QueryStart, output_csv.QueryEnd)),columns =['Name' ,'Start Position', 'End Position'])
seq_record = SeqIO.parse(open('blast_input.fasta'),'fasta')


# Locus Finder

Locus_df = pd.DataFrame(columns = ['Locus Start Position', 'Locus End Position' , 'Locus', 'Query Match'])

for seq_object in seq_record :
    for i in range(0,alu_locations_df.shape[0]):
        seq = str(seq_object.seq)
        Query_Start = alu_locations_df.at[i, "Start Position"] - 1
        Query_End = alu_locations_df.at[i, "End Position"] - 1
        Query_Match = alu_locations_df.at[i, "Name"]
        Locus_Downsteam = seq[max(0,Query_Start-100):Query_Start]
        Locus_Upstream = seq[Query_End + 1:min(Query_End + 100 , len(seq))]
        #print("\n Start Position : " , Query_Start , "\n End Position : " , Query_End ,"\n Match : " , Query_Match ,"\n Transposon : " , seq[Query_Start:Query_End] , "\n")
        Locus = Locus_Downsteam + Locus_Upstream
        #print("\n Locus : " , Locus)
        Locus_df = Locus_df.append({'Locus Start Position' : max(0,Query_Start-100), 'Locus End Position' : min(Query_End + 100 , len(seq)), 'Locus' : Locus, 'Query Match' : Query_Match}, ignore_index = True)

Locus_df.to_csv('LocusIdentifier_Output.csv', index = True)