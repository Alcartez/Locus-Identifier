# Locus Finder Documentation

## Pre-Requisites :
* Blast+
* Biopython
* Python packages : pandas , numpy , csv 

## Methodology :

### Step 1 : Import the Packages

```python
import os
import sys
import pandas as pd 
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
import csv
import subprocess
```

Sys package provides various functions and variables that are used to manipulate different parts of the Python runtime environment. 

Pandas is used for working with Dataframes.

BioPython is used to work with  Biological data like sequences , protein coordinates etc. 

{INSERT_HERE}

### Step 2 : Create Blast Database
```python
db_fasta_raw_filename = sys.argv[2]
local_dbname = "local_alu_db"

makeblastdb_cline = NcbimakeblastdbCommandline(dbtype="nucl", input_file= db_fasta_raw_filename ,out = local_dbname )
stdout,stderr = makeblastdb_cline()
```

The algorithm parses the input filenames as arguments for the python command. The third argument i.e. sys.argv[2] imports the file that is responsible for creating the local blast database. 

We name the database as local_alu_db.
Then the algorithm creates a CLI command to make a local blast database out of a fasta file and executes it through the shell.

### Step 3 : Running Blastn

```python
input_path = sys.argv[1]
input_file = os.path.basename(input_path)
seq_filename = input_file.split(".")[0]
seq_ext = input_file.split(".")[1]
```
The input file is imported through an argument in the python command in the terminal by the user. The input file path string is then broken down into the file name and the extension. 

```python
command_line = "blastn -db local_alu_db -query " + seq_filename + "." + seq_ext + " " +"-out " + seq_filename + ".xml -outfmt 16"
subprocess.call(command_line , shell=True)
```

A string variable is then created using the file name and extension recieved from the user , this string is the CLI input for blastn which runs blastn on the input file against the local database we created for it. Then the string command is run through a subprocess in the system shell which calls the standalone blast to perform the alignment. 

### Step 4 : Parsing the Blast output

```python
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
```

This is a function that parses the blast XML output file and creates a list of rows that can be written into a CSV file. 

```python
blast_ot =  seq_filename + ".xml"
output_csv = parse_blast(blast_ot)
filename = "output.csv"
with open(filename, 'w') as csvfile: 
    csvwriter = csv.writer(csvfile)
    for line in output_csv:
        csvwriter.writerow(line)
```
This code snippet imports the blast XML file , parses it using the parse_blast() function and writes the list of rows into a CSV file. 

### Step 5 : Analyse Blast Output

```python
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
```

This snippet creates a CSV file containing the names and lengths of every Alu element in the local blast database. 

```python
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

output_csv.to_csv("output_2.csv")   
```

This code snippet checks if the length of the blast subject matches with the Alu Elements in the database and checks if the ends of the matches are intact or truncated based of missing nucleotides on the ends. If any of the ends is missing more than 10 nucleotides , it is declared as truncated. 

The snippet then proceeds to create another CSV file named output_2.csv that contains the values for the nucleotide ends. 

### Step 6 : Identifying the locus 

```python
output_csv = pd.read_csv("output.csv")
Alu_start_position_list = output_csv.QueryStart
Alu_end_position_list = output_csv.QueryEnd
alu_locations_df = pd.DataFrame(list(zip(output_csv.Name, output_csv.QueryStart, output_csv.QueryEnd , output_csv.Score , output_csv.Length , output_csv.Expect , output_csv.QueryStrand , output_csv.HitStrand)),columns =['Name' ,'Start Position', 'End Position', 'Score' , 'Length' , 'E-Value','QueryStrand','HitStrand'])
seq_record = SeqIO.parse(open(str(seq_filename+'.fasta')),'fasta')
```

The snippet reads the output.csv file and makes a list of the start position and the end position of the insertions/blast hits.

It then creates a dataframe containing various information about the insertions like Name ,Start Position, End Position, Score , Length , E-Value, QueryStrand direction and HitStrand direction.

It also parses the input fasta file using SeqIO.

```python
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
```

This snippet extracts 100 basepairs upstream i.e. 3' end of the Alu insertion and 100 basepairs downstream i.e. 5' end of the Alu insertion. The former is called the 'upstream' element and the later is called the 'downstream' element. The 'Locus' which is the preinsertion locus is simply the upstream and downstream element joined together. from 3' to 5'. 

These values are then recorded in a CSV file called the LocusIdentifier_Output.csv. 

### Step 7 : Creating Upstream Downstream Files in fasta format. 

```python
lf_data = pd.read_csv("LocusIdentifier_Output.csv")

Upstream_seqlist = lf_data['Upstream (<100)'] 
Downstream_seqlist = lf_data['Downstream (<100)']
```

This part of the script imports the LocusIdentifier_Output.csv as input and extracts the columns containing the upstream and downstream sequences. 

```python
print("\n Creating text file with Upstream sequences with custom sample IDs ... \n")
with open(seq_filename + '_upstream.txt', 'w') as f:
    for i in range(0,len(Upstream_seqlist)) :
        f.write('>SMPL' + str(i) + "_UP " + str(Blast_match_idlist[i]))
        f.write("\n" + str(Upstream_seqlist[i]) + "\n")
print("\n Done ...")
```

The script then extracts each element from the upstream column and writes them as fasta in a text file assigning a SMPL id on the output along with it's original blast match ID. 

```python
print("\n Creating text file with Downstream sequences with custom sample IDs ... \n")
with open(seq_filename + '_downstream.txt', 'w') as f:
    for i in range(0,len(Downstream_seqlist)) :
        f.write('>SMPL' + str(i) + "_DWN " + str(Blast_match_idlist[i]))
        f.write("\n" + str(Downstream_seqlist[i]) + "\n")
print("\n Done ...")
```

The same is repeated for downstream elements i.e. the script then extracts each element from the downstream column and writes them as fasta in a text file assigning a SMPL id on the output along with it's original blast match ID. 