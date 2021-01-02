#!/usr/bin/python3
import os,subprocess,shutil,re,sys
import matplotlib.pylab as plt
from collections import Counter

#########################################################################
#This program uses 4 functions. These are listed below:

#FUNCTION1: Function for yes and no questions, with error proofs
def yesno(question): #The input of this function is question 
  while True:
    ans = input(question).lower().rstrip() #asks user for the input. The input string is lowercased and spaces at the end are removed 
    if ans == 'y': #if yes, the function returns True
      return True
    if ans == 'n': #if no, the function returns False
      return False
    else: #If the input is other than 'y' or 'n', error message is printed and asks for the input again
      print('Incorrect input. Please enter \'y\' for yes and \'n\' for no.'+ '\n') 

#########################################################################
#FUNCTION2: Function for fishing out user-specified number of most similar sequences  
def topseq(question):
  topseq=[] #Contains user-specified number of most similar sequences 
  while True:
    num = input(question).rstrip() #Asks user to enter the number
    if re.search(r'\d',num): #if the input is a digit:
      for n in range(0,int(num)): #For range from 0 to user-specified number:
        accession = accession250_list[n] #moves from first accession number to user-specified position of accession number
        for sequence in sequence_split: #For each sequence in the dataset
          if accession in sequence: #If accession number is found, the sequence is stored in the list 'topseq'
            topseq.append('>'+sequence) #'>' is added to the beginning of every sequence
            ''.join(topseq) #list is joint into a string
      break #loop is broken
    print('Incorrect input. Please enter an integer.' +'\n') #Error message 
  return topseq #the function returns user-specified number of most similar sequences
 
####################################################################### 
#FUNCTION3: Function for counting number of sequences per species, removing user-specified species from the dataset.   
def species(dataset):

  species_list= [] #contains the name of all species with repeats 
  name = re.findall(r'\[.+ .+\]',str(dataset)) #Finds all species names within the dataset
  for n in name:
    species_list.append(n.replace('[','').replace(']','')) #only the species names without square brackets are stored inthe list 'species'
  norepeat_species_list = set(species_list) #contains non-rebundant list of species 
  print('\n'+'*****There are ' +str(len(norepeat_species_list))+ ' different species present in the dataset.') #shows the user the number of different species by counting the length of non-redundant list 
    
  #Asks user if s/he would like to see the number of sequences per species
  choice = yesno("Would you like to see all " +str(len(norepeat_species_list))+ " different species and the number of sequences per species? (y/n) ")
  if choice == True: #If yes:
    counted_species =Counter(species_list) #Counter function counts the number of items with the same name in the species list and generates a dictionary
    for n in counted_species:
      print('\t'+ n +'\t'+ str(counted_species[n])) #prints the species and the number of sequences per species
  else: 
    pass
    
  #Asks user if s/he would like to exclude any species from the current dataset 
  choice = yesno('Would you like to remove any species from the dataset? (y/n) ') 
  if choice == False:
    print('\n'+ '*****There are currently ' +str(dataset.count('>'))+ ' sequences in the dataset.')
    seq_modified = dataset #The original dataset without any modifications is saved into the variable seq_modified
    pass
  if choice == True:#If the answer is yes, the program asks user to enter the name of the species s/he would like to remove from the dataset

    seq_modified_split= re.split(r'>',dataset) #the fasta file is split into list by '>'
    seq_modified_split.pop(0) #the first element in the list is blank. This is removed
   
    while True:
      name = input('Please enter the name of the species you would like to remove. If you are done, please enter \'n\' for no: ').rstrip()
      count=0 #counting
      if name == 'n': #If the user is done removing species from dataset, the inner loop is broken 
        break
      if name in norepeat_species_list: #if the specified species name is in the species list:
        norepeat_species_list.remove(name) #The species is removed from the non-rebundant list of species 
        for sequence in seq_modified_split[:]: 
          if re.search(name,str(sequence)): #For every seqence, if the specified species name is present:
            count+=1 #1 is added to the variable count
            seq_modified_split.remove(sequence)
            
        print('Species ' +name+ ' has been removed from the dataset. '+str(count)+ ' sequences have been removed.')
        print('\n'+ '*****There are currently ' +str(len(seq_modified_split))+ ' sequences in the dataset.' +'\n') #number of sequences after modification
      else:
        print('The species ' +name+ ' is incorrect. Please enter the name of species again. Please make sure that the first letter is in upper case. If you are done, enter \'n\' for no.' +'\n') #Error message
            
    seq_modified_split = ['>'+seq for seq in seq_modified_split] #Adding '>' to beginning of every sequence in the new modified dataset
    seq_modified = ''.join(seq_modified_split) #seq_modified contains user-specified number of species 

  return seq_modified  
  
#####################################################################
#FUNCTION4: Function for counting number of sequences per genus, removing user-specified genus from the dataset.
def genus(dataset):
  genus_list = [] #rebundant list of genus
  name = re.findall(r'\[.+ .+\]',dataset) #Finds all species names within the dataset
  for n in name:
    genus = n.split(' ')[0].replace('[','') #genus is sliced from the species name
    genus_list.append(genus) #genus is added to the list
  norepeat_genus_list = set(genus_list) #contains non-rebundant list of genus
  print('\n'+'*****There are ' +str(len(norepeat_genus_list))+ ' different genus present in the dataset.')
      
  #Asks the user if s/he would like to see the number of sequences per genus
  choice = yesno("Would you like to see all " +str(len(norepeat_genus_list))+ " different genus and the number of sequences per genus? (y/n) ")
  if choice == True:
    counted_genus =Counter(genus_list) #Counter function counts the number of items with the same name in the genus list and generates a dictionary 
    for n in counted_genus:
      print('\t'+ n +'\t'+ str(counted_genus[n])) #prints the genus and the number of sequences per genus 
  else:
    pass

  #Asks the user if s/he would like to remove any genus from the dataset
  choice = yesno('Would you like to remove any genus from the dataset? (y/n) ') 
  if choice == False:
    print('\n'+ '*****There are currently ' +str(dataset.count('>'))+ ' sequences in the dataset.')
    seq_modified = dataset #The original dataset without any modifications is saved onto the variable seq_modified
    pass
  if choice == True:
    seq_modified_split= re.split(r'>',dataset) #the fasta file is split into list by '>'
    seq_modified_split.pop(0) #the first element in the list is blank. This is removed

  #Asks user to specify the name of genus s/he would like to remove from the dataset
    while True:
      name = input('Please enter the name of the genus you would like to remove. If you are done, please enter \'n\' for no: ').rstrip()
      count=0 #counting
      if name == 'n': #If the user is done removing genus from dataset, the inner loop is broken 
        break
      if name in norepeat_genus_list: #if the specified genus is in the genus list:
        norepeat_genus_list.remove(name) #The genus is removed from the non-rebundant list of genus 
        for sequence in seq_modified_split[:]: 
          if re.search(name,str(sequence)): #For every seqence, if the specified species name is present:
            count+=1 #1 is added to the variable count
            seq_modified_split.remove(sequence)
        print('Genus ' +name+ ' has been removed from the dataset. '+str(count)+ ' sequences have been removed.')
        print('\n'+ '*****There are currently ' +str(len(seq_modified_split))+ ' sequences in the dataset.' +'\n') #number of sequences after modification
      else:
        print('The species ' +name+ ' is incorrect. Please enter the name of species again. Please make sure that the first letter is in upper case. If you are done, enter \'n\' for no.' +'\n') #Error message
            
    seq_modified_split = ['>'+seq for seq in seq_modified_split] #Adding '>' to beginning of every sequence in the new modified dataset
    seq_modified = ''.join(seq_modified_split) #seq_modified contains user-specified number of genus
    
  return seq_modified

####################################################
################# STEP 1 ###########################
####################################################

#clearing the screen 
os.system('clear')

#Program starts here!
print('*****Welcome! This program analyses protein sequences across different genus or species!')

#Asks user for the full path to the space s/he would like to save the output files 
#While loop checks that the path exists and if not, it asks user for the path again until correct one is given 
while True:
    outpath = input('Enter the absolute path to the file you would like to save all outputs: ').rstrip()
    if os.path.exists(outpath): #If the path provided exists, the user specified path is saved in the variable 'outpath'
      print('*****All output files will be saved in ' +outpath + '\n' )
      break #while loop is broken
    print('*****The path does not exist. Please try again. Don\'t forget \'/\' at the beginning and end of the path. Example: /localdisk/home/UNN/file1/'+ '\n' ) #If not, error message is given and asks for the correct path

#checks the the file 'analysis' doesn't exist. If it does, removes the file 
if os.path.isdir(outpath+'/analysis'):
  shutil.rmtree(outpath+'/analysis') 
else:
  pass
  
#Creates 'analysis' file in user-specified location 
os.mkdir(outpath+'/analysis') 

###################################################

#User specifies the taxonomy, protein and maximum number of sequences
#yesno function double checks these with the user
organism = input("Please enter the taxonomic group: ") # Asks to specify the taxonomic group
choice = yesno("Is " +organism+ ' correct? (y/n)') #double checks with the user if the specified taxonomic group is correct
if choice == True: #If 'y', the first input is stored in the variable 'organism'
  pass
else:
  organism = input("Please enter the correct taxonomic group: ") #If 'n', program asks user for the taxonomic group again
     
protein= input("Please enter the protein family: ") #Asks to specify the protein family
choice = yesno('Is ' +protein+ ' correct? (y/n)')#Double checks that the specified protein family is correct 
if choice == True: #If 'y', the first input is stored in the variable 'protein'
  pass
else:
  protein = input("Please enter the correct protein family: ") #If 'n', program asks user for the protein again
    
max_seq = input("Please enter maximum number of sequences: ") #Asks user to specify the maximum number of sequences 
choice = yesno('Is ' +max_seq+ ' the correct maximum number of sequencies? (y/n)') #Double checks that the specified max number is correct
if choice == True: #If 'y', the first input is stored in the variable 'max_seq'
  pass
else:
  max_seq = input("Please enter the correct maximum number of sequences: ") #If 'n', program asks user for the max number again

print('\n'+ '*****Great, let\'s create the dataset')

#################################################

#Refining the dataset further: Gives the user an option to include or exclude low quality, partial, hypothetical, predicted and isoform of sequences.
#The yesno function asks the user until corret input is given 
#partial sequence
choice = yesno("Would you like to include partial sequences? (y/n) ")
if choice == True: #if yes, nothing is assigned to the variable 
  partial = ' '
else:
  partial = 'NOT partial' #If no, NOT partial is assigned to the variable. This variable will be used later in esearch.

#low-quality sequences
choice = yesno("Would you like to include low-quality sequences? (y/n) ")
if choice == True:
  low = ' ' #if yes, nothing is assigned to the variable 
else:
  low = 'NOT low quality' #If no, NOT low quality is assigned to the variable. This variable will be used later in esearch.

#Hypothetical sequences
choice = yesno("Would you like to include hypothetical sequences? (y/n) ")
if choice == True:
  hypoth = ' ' #if yes, nothing is assigned to the variable 
else:
  hypoth = 'NOT hypothetical' #If no, NOT hypothetical is assigned to the variable. This variable will be used later in esearch.

#Predicted sequences
choice = yesno("Would you like to include predicted sequences? (y/n) ")
if choice == True:
  predict = ' ' #if yes, nothing is assigned to the variable 
else:
  predict = 'NOT predicted' #If no, NOT predicted is assigned to the variable. This variable will be used later in esearch.

#Isoform of sequences
choice = yesno("Would you like to include isoform of sequences? (y/n) ")
if choice == True:
  isoform = ' ' #if yes, nothing is assigned to the variable 
else:
  isoform = 'NOT isoform' #If no, NOT isoform is assigned to the variable. This variable will be used later in esearch.

##################################################

#Following user's choices, program downloads sequences from NCBI
print('\n'+ "*****Currently downloading sequences from NCBI. Please wait a few moments.")
subprocess.call('esearch -db protein -query " {} [ORGN] {}*[PROT] {} {} {} {} {} " |efetch -format fasta > {}/analysis/sequences.fasta'.format(organism, protein, partial, hypoth, predict, isoform, low, outpath ), shell=True) #variable orgnaism, protein, partial, hypoth, predict, isoform, low, outpath are placed into {}s in that order 
#The sequences are stored in fasta format in 'analysis' directory, which is in a space specified by the user at the beginning

seq = open(outpath+"/analysis/sequences.fasta").read() #the dataset is read and stored in the variable 'seq'
count = seq.count('>') #counting number of sequences

#Checking that the number of sequences is less than the maximum number. 
if count > int(max_seq): #if the number of sequences is bigger than the maximum threshold set by the user, program exits and asks user to refine the dataset
  print('\n'+ "*****There are " +str(count)+ ' sequences in the dataset, which is more than the maximum, ' +str(max_seq) + ' sequences.'+'\n'+ '*****Please redefine the dataset. Program will exit.' )
  shutil.rmtree(outpath+'/analysis/') #The 'analysis' directory is removed 
  sys.exit() #The program is exited 
  
if count == 0: #if there are no sequences, program quits and asks user to redefine the search
  print('\n' +'*****There are no sequences in the dataset. Please redefine the search. Program will exit.')
  shutil.rmtree(outpath+'/analysis/') #The 'analysis' directory is removed 
  sys.exit() #The program is exited
  
else:
  print('\n'+'*****There are ' +str(count)+ ' sequences in the dataset.') #If the number is sequences is less than the max number, program continues.
  
##################################################
############### STEP 2 ###########################
##################################################  
  
#Asking user if s/he would like to display sequences per genus or species 
while True:
  ans = input('Would you like to display sequences per genus or per species? Please enter genus or species: ').lower().rstrip()
  if ans == 'species' or ans == 'genus':
    if ans == 'species': #if species, the program calls function species(). This function displays different species and number of sequences per species. It also asks user to enter the name of species s/he would like to remove from the dataset 
      seq_modified = species(seq) #the modified dataset is assigned to the variable seq_modified 
    if ans == 'genus': #if genus, the program calls function genus().  This function displays different genus and number of sequences per genus. It also asks user to enter the name of genus s/he would like to remove from the dataset 
      seq_modified = genus(seq) #The modified dataset is assigned to the variable seq_modified 
    break #The while loop is broken
  print('Incorrect input. Please enter \'species\' or \'genus\'. ' +'\n') #error message for the main outer loop  
  
#creating a file containing only user specified sequences 
outfile = open(outpath+"/analysis/sequence.fasta",'w') #creating a file called sequence.fasta 
for line in seq_modified: #All lines in the modified dataset are written 
  outfile.write(line) 
outfile.close() #close file
print('\n'+'*****The modified dataset containing sequences is at ' +outpath+'/analysis. The file is called \'sequence.fasta\' .' + '\n')

#removing the first dataset
os.remove(outpath+'/analysis/sequences.fasta')

#Asking user if s/he would like to continue
choice = yesno("Would you like to continue? If no, the program will quit. (y/n) ") 
if choice == True: #if yes, the program continues
  print('Great! Continuing to the next step...')
else: #If no, program exits 
  print('\n'+ '*****Very well. Thank you for using my program. All output files have been deleted. See you soon!')
  shutil.rmtree(outpath+'/analysis/') #if the user doesn't wish to continue, the outfile files are deleted and program is shut
  sys.exit() #The program is exited      
        
##################################################
################## STEP 3 ########################
##################################################

#Counting the number of sequences. if it's aobve 250, blastp is performed 
sequence = open(outpath+'/analysis/sequence.fasta').read() #Modified dataset is read and stored in the variable 'sequence' 
sequence_split = sequence.split('>') #The dataset is split into individual sequences 
sequence_split.pop(0) #First element which is blank, is removed 

if sequence.count('>') > 250: #if there are more than 250 sequences in the dataset
  print('\n'+'*****There are more than 250 sequences in the dataset. We will perform Clustal-Omega and BLASTp to identify 250 most similar sequences now.')  
  
  #Asking user to select the species s/he would like to align the sequences against
  while True:
    name = input('Please enter the name of the reference genus or species you would like to align the sequences: ').rstrip()
    if re.search(name, sequence): #if the specified species or genus name is in dataset
      target_name = name #The user-specified genus or species name is stored in the variable target_name
      break #The loop is broken
    print(str(name)+ ' is incorrect. Please enter the name again. Please make sure the first letter is in upper case. E.G. Nestor notabilis' +'\n') #Error message
  
  consensus = [] #List contains all sequences for the user-specified genus or species   
  for sequence in sequence_split: #For every sequence in the dataset:
    if target_name in sequence: #If the user-specified genus or species is present, the sequence is added to the list
      consensus.append('>'+sequence) #'>' is added to the beginning of every sequence added 
  
  if len(consensus) > 1: #If the number of sequence of the reference genus or species is greater than 1: 
    outfile = open(outpath+'/analysis/consensus_in.fasta','w') #File is created 
    for line in consensus: #All sequences of the reference genus or species is written to this new file 
      outfile.write(line)
    outfile.close() #File is closed 
    
    #Multiple sequence alignment using Clustal-Omega
    subprocess.call('clustalo -i {}/analysis/consensus_in.fasta -o {}/analysis/consensus_out.fasta --outfmt fasta --force' .format(outpath,outpath),shell=True)
    #Consensus is sequence from the multiple sequence alignment is found using cons
    subprocess.call('cons -sequence {}/analysis/consensus_out.fasta -outseq {}/analysis/consensus.fasta -osformat2 fasta'.format(outpath,outpath),shell=True)
    #Files are deleted 
    os.remove(outpath+'/analysis/consensus_in.fasta')
    os.remove(outpath+'/analysis/consensus_out.fasta')
    
  else: #if the number of sequence of the reference genus or species is 1:
    outfile = open(outpath+'/analysis/consensus.fasta','w') #output file is created
    for line in consensus: #One sequence is written to the file 
      outfile.write(line)
    outfile.close() #file is closed 
    
  #blast database is made using the original number of sequences 
  subprocess.call('makeblastdb -in {}/analysis/sequence.fasta -dbtype prot -out {}/analysis/organism'.format(outpath,outpath), shell=True)
  #blastp performed, the output shows top 250 matches  
  print('\n'+'*****Performing BLASTp now. Please wait a few moments...')
  subprocess.call('blastp -db {}/analysis/organism -query {}/analysis/consensus.fasta -outfmt 6 -num_alignments 250  > {}/analysis/blastoutput.out'.format(outpath, outpath, outpath) , shell = True)
  
  #removing blast databases
  os.remove(outpath+'/analysis/organism.phr')
  os.remove(outpath+'/analysis/organism.pin')
  os.remove(outpath+'/analysis/organism.psq')
  
  #Finding 250 most similar sequences 
  blast = open(outpath+'/analysis/blastoutput.out').read().replace('>','') #blast output is opened and read 
  blast_split = blast.split('\n') #The output is split into indexes
  blast_split.pop(-1) #The last element which is blank is removed 

  accession250_list=[]#For every blast output, the accession number is stored in this list 
  for sequence in blast_split: 
    accession250_list.append(sequence.split('\t')[1]) #The subject accession number is in the second column

  seq250 = [] #this list will contain 250 most similar sequences 
  for n in range(0,250): #Only first 250 accession numbers are used 
    accession = accession250_list[n] 
    for sequence in sequence_split:
      if accession in sequence: #if the accession number is found in the sequence:
        seq250.append('>'+sequence)  #the sequence is added to the list seq250
  seq_modified = ''.join(seq250) #The modified dataset is made into string 

  #Creating new fasta file containing only top 250 most similar sequence 
  outfile = open(outpath+"/analysis/sequence.fasta",'w') #new file is created and opened 
  for line in seq_modified: 
    outfile.write(line) #All lines in the seq_modified are written to the new file 
  outfile.close()
  
  print('\n'+'*****The dataset has been updated to contain 250 sequences. This file can be found at: ' +outpath+'/analysis/sequence.fasta.')
  print('\n'+'*****We will now perform Clustal-Omega to align 250 sequences.' +'\n')
  
  #Multiple sequence alignment of 250 sequences
  subprocess.call('clustalo -i {}/analysis/sequence.fasta -o {}/analysis/clustal.out --force' .format(outpath, outpath) , shell=True) 
  
else: #If the number of sequences in the dataset is less than 250, those sequences are aligned using Clustal-Omega. The output is ordered from most to least similar
  print('\n'+ '*****Currently performing multiple sequence alignment using Clustal-Omega. Please wait a few moments...')
  subprocess.call('clustalo -i {}/analysis/sequence.fasta -o {}/analysis/clustal.out --output-order=tree-order --force' .format(outpath, outpath) , shell=True) #Multiple sequence alignment of 250 sequences

clustal = open(outpath+'/analysis/clustal.out').read() #The clustalo output is opened and read 
clustal_split = clustal.split('>') #clustalo output is split into sequences 
clustal_split.pop(0) #first element, which is blank is removed 

accession250_list = [] #contains accession number of sequences from most to least similar 
for sequence in clustal_split: #For every sequence in clustalo output
  header = sequence.split('\n')[0] #the header is the first element when split by new line
  accession = header.split(' ')[0] #the accession number is first element when the header is split by a space 
  accession250_list.append(accession) #the accession number is added to the list in order. This will be used for next step when the user specifies number of sequences to scan for motifs 
  
  
#########################################################
################## STEP 4 ###############################
#########################################################

#Plotting the similarity of aligned sequences 
#plotcon plots the similarity of aligned sequences. The plot is saved 
print('\n'+'*****Now we will perform some analysis on these aligned sequences. First we will plot the similarity of aligned sequences.')  
subprocess.call('plotcon -sequence {}/analysis/clustal.out -winsize 5 -graph svg -goutfile {}/analysis/similarity' .format(outpath, outpath) ,shell=True)

#Asking user if s/he would like to view the plot 
choice = yesno('Would you like to view and save the plot? (y/n) ')
if choice == True: #If yes, the plot is displayed onto the screen 
  subprocess.call('display {}/analysis/similarity.svg' .format(outpath) , shell=True ) 
  print('\n'+'*****The similarity plot of aligned sequences can be found in ' +outpath+'/analysis. The plot is called similarity.svg'+'\n')
else: #if no, the plot is removed 
  os.remove(outpath+'/analysis/similarity.svg')
  
#########################################################

#Determining level of conservation by calculating the proportion of unique amino acids in each position
#Proportion of unique amino acids in each column
print('\n'+'*****The level of conservation between protein sequences will be determined by calculating the proportion of unique amino acids in each position of the sequence. The proportion will be plotted.' +'\n')
aligned = open(outpath+'/analysis/clustal.out').read() #the clustal output is opened and read 
aligned_split = re.split(r'>',aligned) #sequences are split into indexes using '>'
aligned_split.pop(0) #First two elements, which are empty are removed
aligned_split.pop(0)

noheader_seq = [] #this list contains one sequence per line without header
for seq in aligned_split:
  aa=seq.split('\n')[1:] #fist element is the header. Excluding the header, the sequences are stored in the variable
  noheader_seq.append(''.join(aa)) #the sequences are joint into one string and added to the list 

#Calculating the proportion of unique amino acids in each column 
unique_aa = [] #All 
for col in range(len(noheader_seq[0])): #for every column in the sequence:
  aa = []
  for seq in noheader_seq: #For every sequence:
    if seq[col] != '-': #if gap is found, it is ignored 
      aa.append(seq[col]) #all amino acids in one column are added to the list
  unique_aa.append(len(set(aa))) #for each column, non-rebundant list of amino acids are added to the list unique_aa

window = 10 #sliding window methd 
proportion = []
for i in range(len(unique_aa)-window):
  win = unique_aa[i:i+window] #window slides across the sequence
  fraction = sum(win)/len(win) #fraction at each position is calculcated 
  proportion.append(fraction) #this fraction is added to the list proportion 

plt.figure(figsize=(15,8)) #figure size is set to 15 by 8
plt.plot(proportion) #the proportion at each position is plotted 
plt.title('Proportion of unique amino acids in each position') #title
plt.xlabel('Residue position') #x-axis label
plt.ylabel('Proportion of unique amino acid') #y-axis label

#Asking user if s/he would like to view the plot 
print('\n') #pritning new line
choice = yesno('Would you like to view and save the proportion of unique amino acids? (y/n) ')
if choice == True: #if yes, the plot is saved and opened 
  plt.savefig(outpath+'/analysis/unique_aa', transparent=True)
  print('\n' +'The plot is saved in: ' +outpath+'/analysis/unique_aa.png')
  plt.show() #opening connection 
  plt.close() #closing connection 
else:
  pass

#Obtaining basic information about the multiple sequence alignment using infoalign 
#By default, the consensus sequence is calculated and all sequences are compared to this reference sequence 
print('\n'+'*****For more information on the multiple sequence alignment, please see \'alignmentInfo.txt\' in: ' +outpath+'/analysis.')
subprocess.call('infoalign {}/analysis/clustal.out -only -name -seqlength -alignlength -gaps -idcount -simcount -diffcount -change -outfile {}/analysis/alignmentInfo.txt'.format(outpath,outpath),shell=True)

#####################################################
############## STEP 5 ###############################
#####################################################

sequence = open(outpath+'/analysis/sequence.fasta').read() #The dataset is opened and read 
sequence_split = re.split(r'>',sequence) #The dataset is split into individual sequences
sequence_split.pop(0) #The first element which is blank is removed

#To peform further analyses, the user is given option to select a subset of sequences 
print('\n'+ '*****Now let\'s perform some more analyses.')
print('\n'+ '*****There are currently ' +str(sequence.count('>'))+ ' sequences in the dataset.') #reminds user the number of sequences in the current dataset
choice = yesno('Would you like to select a subset of sequences to perform further analyses? (y/n) ') #Gives an option to select a subset of sequences 
print('\n') #prints new line 
if choice == False: #if no, the analyses are done on original dataset
  seq_modified = sequence
if choice == True: #if yes, the program presents different ways of selecting a subset of sequences

  while True:
    ans = input('If you would like to select sequences based on genus, enter \'genus\' and for species, enter \'species\'. If you would like to select most similar sequences based on the previous multiple sequence alignment, please enter \'alignment\': ').lower().rstrip()
    if ans == 'species' or ans == 'genus' or ans == 'alignment':
      if ans == 'species': #if user selects species, program calls function species(). This function presents different species and number of sequences per species and asks user to enter the name of species s/he would like to remove from the dataset.   
        seq_modified = species(sequence) #The modified dataset is saved in the variable seq_modified.
      if ans == 'genus': #if user selects genus, program calls function genus(). This function presents different genus and number of sequences per genus and asks user to enter the name of genus s/he would like to remove from the dataset. .  
        seq_modified = genus(sequence) #The modified dataset is saved in the variable seq_modified
      if ans == 'alignment': #If user selects alignment, program calles function topseq(). This function selects user-specified number of most similar sequences from the previous multiple sequence alignment. 
        seq_modified = topseq('How many most similar sequences would you like to use? : ') #The modified dataset is saved in the variable seq_modified
      break
    print('Incorrect input. Please enter \'species\' or \'genus\' or \'alignment\'. ' +'\n') #Error message for the while loop

#Creating a file containing the user-specified subset of sequences
outfile = open(outpath+"/analysis/subset_sequence.fasta",'w') #creating a file called subset_sequence.fasta 
for line in seq_modified: #All lines in the modified dataset are written 
  outfile.write(line) 
outfile.close() #file closed

#Multiple sequence alignment of subset of sequences is performed using Clustal-Omega for later use 
print('\n' +'*****Great. Now, the subset of sequences will be aligned using Clustal-Omega for later use. Please wait a few moments... ' +'\n')
subprocess.call('clustalo -i {}/analysis/subset_sequence.fasta -o {}/analysis/subset_clustal.out --force' .format(outpath,outpath) ,shell =True)

#Outlining the different analyses to the user
print('\n' + '*****Great! In this program, you will have the option to perform 4 different analyses. These are: ' +'\n\t'+ '1. Scan for known motifs.' +'\n\t'+ '2. Scan for a specific pattern.' +'\n\t'+ '3. Predict and plot transmembrane segments.' +'\n\t'+ '4. Plot hydropathy.')

#Opening and reading the file containing the subset of sequences for further analyses 
subset = open(outpath+'/analysis/subset_sequence.fasta').read()
subset_split = subset.split('>') #The subset is split into individual sequences 
subset_split.pop(0) #The first element, which is empty, is removed 

####################################################
#Four analyses are performed
#ANALYSIS 1 : scanning sequences for any known motifs 
print('\n') #print new line 
choice = yesno('*****Would you like to scan these sequences for any known motifs? (y/n): ') #Gives user option to scan sequences for any known motifs
if choice == False: #If no, this analysis is skipped
  pass
else: #If yes, the analysis is performed 
  print('\n'+ '*****Great. Let\'s scan for any known motifs.')
  motif_dict={} #all outputs will be saved in this dictionary 
  for sequence in subset_split: #for every sequence in the dataset:
    sequence = '>'+sequence #'>' is added to the beginning of the sequence since patmatmotifs requires this symbol
    out_seq = open(outpath+'/analysis/in_seq.fasta','w') #File is created 
    for line in sequence: #One sequence is written to the file 
      out_seq.write(line)
    out_seq.close() #File closed 
    
    #patmatmotifs scans the file containing one sequence for any known motifs. It creates an output file called motifs.out
    subprocess.call('patmatmotifs -sequence {}/analysis/in_seq.fasta -outfile {}/analysis/motifs.out'.format(outpath,outpath),shell=True) 
    
    motif_out = open(outpath+'/analysis/motifs.out').read() #The output from patmatmotifs is opened and read
    name =  re.findall(r'# Sequence: .+',motif_out) #finds the sequence name from the output file 
    all_motif = re.findall(r'Motif = .+',motif_out) #finds the motif name from the output file 
    for motif in all_motif:
      motif_dict[str(name)] = str(motif) #the name is stored as the key, and the motif name is stored as the value of that key
  
    os.remove(outpath+'/analysis/in_seq.fasta') #input files for patmatmotifs are deleted 
    os.remove(outpath+'/analysis/motifs.out') #output files for patmatmotifs are deleted 
    
  print('\n') #prints a new line
  if len(motif_dict) > 0: #If motifs were detected:
    for motif in motif_dict: #For every element in the dictionary:
      print(motif.replace('[\'#','').replace('\']','') +', has '+ motif_dict[motif]) #the sequence name and the motif name are printed onto the screen 
      
    choice=yesno('Would you like to save these results? (y/n)') #Gives user the option to save the results 
    if choice == True:
      outfile = open(outpath+'/analysis/motif.out','w') #If yes, file is created 
      for motif in motif_dict: #every element in the dictionary is saved onto the file
        outfile.write(motif.replace('[\'#','').replace('\']','') +', has '+ motif_dict[motif]+'\n')
      outfile.close() #file closed 
      print('\n' + '*****The results are stored in ' +outpath+ '/analysis/motif.out')
    else: #if not, the file is not saved 
      pass
      
  else: #outer conditional: if no motifs were identified, prints the messgage to the screen 
    print('*****No motifs were identified.' +'\n')
                                                                                              
#####################################################
#ANALYSIS2: Pattern searching 
print('\n')
choice = yesno('*****Would you like to search for a particular pattern in the subset of sequences? (y/n): ') #Gives user an option to search for patterns in sequences 
if choice == False: #if no, this analysis is skipped 
  pass
else: #If yes, the analysis is performed 

  print('\n' +'*****The search pattern is written in PROSITE style.') 
  choice = yesno('Would you like to open the manual on pattern syntax? (y/n): ') #Asks user if s/he would like to see the manual on pattern syntax
  if choice == True: #If yes, program opens the PROSITE website, which has the manual 
    subprocess.call('firefox https://prosite.expasy.org/scanprosite/scanprosite_doc.html',shell=True) 
  else: #If no, skipped to the next step
    pass
    
  pattern = input('Please enter the search pattern in PROSITE style: ').upper() #User specifies the pattern
  choice = yesno('Is ' +pattern+ ' correct? (y/n): ') #Double checks the pattern with the user 
  if choice == False: #if the pattern is incorrect, asks user to enter the pattern again 
    pattern = input('Please enter the search pattern in PROSITE style: ')
  else:
    pass
  #fuzzpro is used to search for the user-specified pattern in the sequences. The output file is called fuzzpro.out. 
  subprocess.call('fuzzpro -sequence {}/analysis/subset_sequence.fasta -pattern {} -outfile {}/analysis/fuzzpro.out' .format(outpath,pattern,outpath) ,shell = True)
  
  fuzzpro = open(outpath+'/analysis/fuzzpro.out').read() #The fuzzpro output file is opened and read 
  hitcount = re.search(r'Reported_hitcount: \d+',fuzzpro).group().replace('Reported_hitcount: ','') #the hit count is saved in the variable 
  
  if hitcount == '0': #if the hitcount is 0, the fuzzpro outfile is removed 
    print('\n' +'*****Pattern not found.') #message printed to the screen 
    os.remove(outpath+'/analysis/fuzzpro.out') #output file removed 
    
  else: #if the hitcount is not 0:

    end_position =re.findall(r'\d+ p',fuzzpro) #The end residue numbers of the pattern match are found 
    hitcount = re.search(r'Reported_hitcount: \d+',fuzzpro).group().replace('Reported_hitcount: ','') #hitnumber is found 
    total_seq = re.search(r'Total_sequences: \d+',fuzzpro).group().replace('Total_sequences: ','') #total number of sequence is found
    total_length = re.search(r'Total_length: \d+',fuzzpro).group().replace('Total_length: ','') #total length of the sequences is found 
    pattern= re.search(r'\. .+',fuzzpro).group().replace(' ','').replace('.','') #the pattern is found 
    seq_len = int(total_length)/int(total_seq) #the average length of the sequence is found by total length of sequences / total number of sequences
    
    start_end = [] #the start and end residue numbers for every pattern match are stored in the list 
    for line in end_position:
      end = line.replace(' ','').replace('p','') #spaces and 'p' are removed so 'end' variable only contains the residue number 
      start = int(end) - len(pattern) #the start residue number is found by end residue number - length of the pattern
      both = str(start)+','+str(end) #both variable has both the start and end residue number separated by a comma
      start_end.append(both) #This is added to the list start_end 
    
    count = []
    for n in range(len(start_end)): #For every element in the list start_end:
      start = start_end[n].split(',')[0] #the start residue number is the first element 
      end = start_end[n].split(',')[1] #the end residue vnumber is the second element
      pattern_range = range(int(start),int(end)+1,1) #fills in the numbers between start and end residue number e.g. if start=1 and end=4, pattern_range = 1,2,3,4
      for n in pattern_range: #Each residue number where pattern was found is added to the list count 
        count.append(n)
    occurence = Counter(count) #Counter function counts the number of occurence in each residue position. The key of occurrence dictionary contains residue numbers and the value contains the corresponding number of occurrence.
    
    dict0 ={} #dict0 is empty 
    for n in range(int(seq_len)): #the key of dict0 contains residue number from 0 to the average length of a sequence
      dict0[n] = 0 #for now, all keys have value '0'
    
    for x in dict0.keys(): #for every key in the dict0:
      if x in occurence.keys(): #If the residue number matches the key of the occurrence dictioinary:
        value = occurence[x] #the variable value contains the value of the occurrence dictionary
        dict0[x] = value/int(hitcount) #This value is divided by the total hitcount and stored in the corresponding residue number 
    
    plt.figure(figsize=(15,8)) #plotting results. The figure has size 15 by 8
    plt.plot(list(dict0.values())) #The proportion of occurrences in each residue position is plotted 
    plt.title('Occurrence of pattern in each position') #title
    plt.xlabel('Residue position') #x-axis label
    plt.ylabel('Occurence of pattern') #y-axis label
    
    choice = yesno('Would you like to view and save the plot? (y/n) :') #gives user option to view and the plot 
    if choice == True: 
      plt.show() #if yes, shown onto the screen 
      plt.savefig(outpath+'/analysis/pattern.png', transparent=True) #plot is saved 
      plt.close() #closing connection 
      print('\n' +'*****The plot is saved at: ' +outpath+ '/analysis/pattern_plot.png')
    else:
      pass #If not, passed 
    
    choice = yesno('Would you like to save these results? (y/n) ') #gives user option to save the fuzzpro output 
    if choice == True:#if yes, the message is printed 
      print('The fuzzpro output is saved at :' +outpath+'/analysis/fuzzpro.out') 
    else:#if not, the output is removed 
      os.remove(outpath+'/analysis/fuzzpro.out') 
      
#####################################################
#ANALYSIS3: predicting and plotting transmembrane segments in protein sequences 
print('\n')
choice = yesno('*****Would you like to predict and plot transmembrane segments for aligned sequences? (y/n): ') #Gives user an option to plot transmembrane segments  
if choice == False: #if no, this analysis is skipped 
  pass
else: #If yes, the analysis is performed 

  #tmap uses the clustalo output to predict and plot transmembrane segments in protein sequences 
  subprocess.call('tmap {}/analysis/subset_clustal.out -out {}/analysis/tmap.out -graph svg -goutfile {}/analysis/tmap' .format(outpath,outpath,outpath) ,shell =True)
  
  #gives user option to view and save the plot 
  choice= yesno('Would you like to view and save the plot? (y/n): ') 
  if choice == False: #if no, the plot is deleted 
    os.remove(outpath+'/analysis/tmap.svg')
  else: #if yes, the plot is displayed onto the sceren and message is printed 
    subprocess.call('display {}/analysis/tmap.svg' .format(outpath) , shell=True )
    print('\n'+'*****The plot is saved at :' +outpath+'/analysis/tmap.svg')
    
  #gives user option to save the tmap output 
  choice = yesno('Would you like to save the output from tmap? (y/n): ')
  if choice == True: #if yes, message displayed 
    print('\n'+ '*****The output is saved at: ' +outpath+'/analysis/tmap.out')
  else: #if no, the output is removed 
    os.remove(outpath+'/analysis/tmap.out')
    
#####################################################
#ANALYSIS4: plotting hydropathy for aligned sequences 
print('\n')
choice = yesno('*****Would you like to plot the hydropathy of aligned sequences? (y/n): ') #Gives user an option to plot hydropathy  
if choice == False: #if no, this analysis is skipped 
  pass
else: #If yes, the analysis is performed 

  subprocess.call('pepwindowall {}/analysis/subset_clustal.out -graph svg -goutfile {}/analysis/hydropathy' .format(outpath,outpath) ,shell =True)
  
  #gives user option to view and save the plot 
  choice= yesno('Would you like to view and save the plot? (y/n): ') 
  if choice == False: #if no, the plot is deleted 
    os.remove(outpath+'/analysis/hydropathy.svg')
  else: #if yes, the plot is displayed onto the sceren and message is printed 
    subprocess.call('display {}/analysis/hydropathy.svg' .format(outpath) , shell=True )
    print('\n'+'*****The plot is saved at :' +outpath+'/analysis/hydropathy.svg')
    
####################################################
####################################################

#Gives option to remove the sequence and alignemnt files 
print('\n' + '*****Sadly, we have come to the end.')
choice = yesno('Would you like to remove all sequence and alignment files generated from this analysis? The analysis outputs will not be deleted. (y/n): ')
if choice == True: #if yes, the files are all removed 
  os.remove(outpath+'/analysis/clustal.out')
  os.remove(outpath+'/analysis/subset_clustal.out')
  os.remove(outpath+'/analysis/sequence.fasta')
  os.remove(outpath+'/analysis/subset_sequence.fasta')
else: #else, the files are kept 
  pass

#End of the program
print('\n' +'*****Thank you for using my program! All outputs are saved in: ' +outpath+'/analysis.'+ ' See you soon!')
    
