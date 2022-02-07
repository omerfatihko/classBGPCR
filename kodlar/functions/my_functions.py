import pandas as pd
from collections import Counter
import collections
import numpy as np
from typing import List, Dict

class Node:
    """define a node class. A node has at most 2 children (left and right), a name,
    and a data. A node with no children is called a leaf"""
    def __init__(self, name, left = None, right = None, rawdata = None, data = None, metadata = None):
        self.rawdata = rawdata
        self.data = data
        self.metadata = metadata
        self.name = name
        self.left = left
        self.right = right

def to_dictionary(nodele: Node) -> Dict:
    pass

def is_leaf(node: Node) -> bool:
    if(not node.left and not node.right):
        return True
    else:
        return False

def get_leaf_rawdata(root: Node) -> Dict:
    """returns a dictionary where keys are leaf names and values are leaf raw data"""
    # If node is null, return
    leafdict = {}
    if (not root):
        return leafdict
    
    # If node is leaf node,
    # key --> leaf name, value --> leaf data
    if (not root.left and not root.right):
        leafdict[root.name] = root.rawdata
        #leaflist.append(root.name)
        #print(root.data,
        #	end = " ")
        return leafdict

    # If left child exists,
    # check for leaf recursively
    if root.left:
        temp = get_leaf_rawdata(root.left)
        leafdict.update(temp)

    # If right child exists,
    # check for leaf recursively
    if root.right:
        temp = get_leaf_rawdata(root.right)
        leafdict.update(temp)
    
    return leafdict

def get_leaf_data(root: Node) -> Dict:
    """returns a dictionary where keys are leaf names and values are leaf data"""
    # If node is null, return
    leafdict = {}
    if (not root):
        return leafdict
    
    # If node is leaf node,
    # key --> leaf name, value --> leaf data
    if (not root.left and not root.right):
        leafdict[root.name] = root.data
        #leaflist.append(root.name)
        #print(root.data,
        #	end = " ")
        return leafdict

    # If left child exists,
    # check for leaf recursively
    if root.left:
        temp = get_leaf_data(root.left)
        leafdict.update(temp)

    # If right child exists,
    # check for leaf recursively
    if root.right:
        temp = get_leaf_data(root.right)
        leafdict.update(temp)
    
    return leafdict

def get_leaf_list(root: Node) -> List:
    """returns a list of node objects (leaves)"""
    leaflist = []
    # If node is null, return
    if (not root):
        return leaflist

    # If node is leaf node,
    # get the node
    if (not root.left and not root.right):
        leaflist.append(root)
        #leaflist.append(root.name)
        #print(root.data,
        #	end = " ")
        return leaflist

    # If left child exists,
    # check for leaf recursively
    if root.left:
        leaflist.extend(get_leaf_list(root.left))

    # If right child exists,
    # check for leaf recursively
    if root.right:
        leaflist.extend(get_leaf_list(root.right))

    return leaflist

def get_node_list(root: Node) -> List:
    "returns a list of node objects (in order traversal)"
    #if node is null, return
    nodelist = []
    if (not root):
        return nodelist

    #left -> root -> right
    else:
        nodelist.extend(get_node_list(root.left))
        nodelist.append(root)
        nodelist.extend(get_node_list(root.right))
        return nodelist

def get_inner_node_list(root: Node) -> List:
    """get inner nodes (nodes that are not leaves) of the binary tree"""
    nodelist = get_node_list(root)
    leaflist = get_leaf_list(root)
    innernodelist = [k for k in nodelist if k not in leaflist]
    return innernodelist

def df_same(df1: pd.DataFrame, ind1: int, df2: pd.DataFrame, ind2: int) -> List:
    """return indexes of rows of two columns where two elements of the columns at the
    row are equal"""
    colnames1 = list(df1.columns)
    colnames2 = list(df2.columns)
    newdf = pd.concat([df1, df2], axis=1) #concatanate dfs to subset them together
    subdf = newdf.loc[newdf[colnames1[ind1]] == newdf[colnames2[ind2]]]
    idx = list(subdf.index)
    return idx

def df_equalto(df: pd.DataFrame, idx: int, ele) -> List:
    """return indexes of rows of a column where the element of the column at the
    row is equal to given element"""
    colnames = list(df.columns)
    subdf = df.loc[df[colnames[idx]] == ele]
    idx = list(subdf.index)
    return idx    
    pass

def df_diff(df1: pd.DataFrame, ind1: int, df2: pd.DataFrame, ind2:int) -> List:
    """return indexes of rows of two columns where two elements of the columns at the
    row are not equal"""
    colnames1 = list(df1.columns)
    colnames2 = list(df2.columns)
    newdf = pd.concat([df1, df2], axis=1) #concatanate dfs to subset them together
    subdf = newdf.loc[newdf[colnames1[ind1]] != newdf[colnames2[ind2]]]
    idx = list(subdf.index)
    return idx

def inner_check(datadict: Dict) -> List:
    """this function will take a dictionary of dataframes. Columns of dataframes are
    0 - sequence, 1 - residue, 2 - frequency, 3 - status. This function will test
    1. residues are the same
    2. status is "C" or conserved (frequency above 0.9).
    We will iterate over keys 2 at a time and check for the given conditions.
    Function returns the indexes of said residues."""
    #get the keys of the dictionary
    keys = sorted(datadict.keys())
    #create a list to hold indexes that satisfy given conditions
    listofidx = []
    #if there is only single dataframe, return indexes of the df
    if len(keys) == 1:
        dfsingle = datadict[keys[0]]
        idx = df_equalto(dfsingle, 3, "C") #list(dfsingle.index)
        idx.sort()
        return idx
    else: #more than one dataframe
        for i in range(len(keys) - 1): #iterate over two elements at a time
            keypair = keys[i:i+2]
            df1 = datadict[keypair[0]] #first df
            df2 = datadict[keypair[1]] #second df
            temp =[df_same(df1, 1, df2, 1), df_equalto(df1, 3, "C"), df_equalto(df2, 3, "C")]
            idx = set(temp[0]).intersection(*temp)
            listofidx.append(idx)
        #return the indexes that satisfy the conditions for all dataframes at the same time
        commonels = list(set.intersection(*[set(x) for x in listofidx]))
        commonels.sort()
        return commonels

def divergence(datadict1: Dict, datadict2: Dict) -> List:
    """this function will take two dictionaries of dataframes. Columns of dataframes are
    0 - sequence, 1 - residue, 2 - frequency, 3 - status. This function will test
    1. residues are not the same
    2. status is "C" or conserved (frequency above 0.9).
    We will iterate over keys 2 at a time and check for the given conditions.
    Function returns the indexes of said residues."""
    listofidx = []
    for value1 in datadict1.values(): #iterate over two dicts
        for value2 in datadict2.values():
            temp = [df_diff(value1, 1, value2, 1), df_equalto(value1, 3, "C"), df_equalto(value2, 3, "C")]
            idx = set(temp[0]).intersection(*temp)
            listofidx.append(idx)
    #return the indexes that satisfy the conditions for all dataframes at the same time
    commonels = list(set.intersection(*[set(x) for x in listofidx]))
    commonels.sort()
    return commonels

def emergence(datadict1: Dict, datadict2: Dict) -> List:
    """this function will take two dictionaries of dataframes. Columns of dataframes are
    0 - sequence, 1 - residue, 2 - frequency, 3 - status. This function will test
    1. residues are not the same
    2. status is "C" or conserved for first dict and "NC" or nonconserved
    for second dict. (frequency above 0.9). We will iterate over keys 2 at a time and 
    check for the given conditions. Function returns the indexes of said residues."""
    listofidx = []
    for value1 in datadict1.values():
        for value2 in datadict2.values():
            temp = [df_diff(value1, 1, value2, 1), df_equalto(value1, 3, "C"), df_equalto(value2, 3, "NC")]
            idx = set(temp[0]).intersection(*temp)
            listofidx.append(idx)
    #return the indexes that satisfy the conditions for all dataframes at the same time
    commonels = list(set.intersection(*[set(x) for x in listofidx]))
    commonels.sort()
    return commonels

def outer_check(datadict1: Dict, datadict2: Dict) -> List:
    """this function will take two dictionaries of dataframes. Columns of dataframes are
    0 - sequence, 1 - residue, 2 - frequency, 3 - status. This function will test
    1. residues are not the same
    2. status is "C" or conserved for first dict (frequency above 0.9).
    We will iterate over keys 2 at a time and 
    check for the given conditions. Function returns the indexes of said residues."""
    listofidx = []
    for value1 in datadict1.values():
        for value2 in datadict2.values():
            temp = [df_diff(value1, 1, value2, 1), df_equalto(value1, 3, "C")]
            idx = set(temp[0]).intersection(*temp)
            listofidx.append(idx)
    #return the indexes that satisfy the conditions for all dataframes at the same time
    commonels = list(set.intersection(*[set(x) for x in listofidx]))
    commonels.sort()
    return commonels

def chisquare(observed, expected):
    """This function takes two lists and perform chi square test using lists,
    for info on chi-square test, visit https://www.statisticshowto.com/probability-and-statistics/chi-square/.
    Lists should have same length, categories should share indexes, degree of freedom is equal to the length of lists.
    returns a tuple (first element chi value, second degree of freedom)."""
    chival = 0
    for i in range(len(observed)):
        a = ((observed[i] - expected[i])**2)/expected[i]
        chival += a 
    return (chival, len(observed) - 1)



def combine_csv(path1,name1, path2, name2):
    """This function combines csv files that contain ortholog protein identifiers.
    Combined orthologs will be used to perform MSA."""

    #read the csv files from specified paths
    df1 = pd.read_csv(path1, header = None)
    df2 = pd.read_csv(path2, header = None)

    #combine the files
    frames = [df1, df2]
    results = pd.concat(frames)

    #write them to the path. File name includes names of both target proteins
    path_name = "/cta/users/ofkonar/work/results/" + name1 + "_" + name2 + "_orthologs.csv"
    results.to_csv(path_name, index = False, header = False)

def select_csv(path, column):
    """This function takes the csv file that contains ortholog identifiers
    and returns the column that has unique identifier as a list. column 
    is an index"""

    #read the csv file from specified path
    df = pd.read_csv(path, header = None)

    #get the specified columns with index
    selected = df.iloc[:, column].tolist()

    return selected

def filter_fasta(selected, unified_fasta_path):
    """This function accepts a list of unique protein ids, search them
    in a fasta database and returns the resulted fasta sequences as a list"""
    #open a list to append fasta sequences
    fasta_seqs = []

    #get fasta sequences sequence by sequence
    with open(unified_fasta_path, 'r') as file:
        counter = 0
        temp = ""
        for line in file:
            if line[0] == '>':
                counter = counter + 1
                if counter == 2:
                    fasta_seq = FilterSeq(temp, selected)
                    if fasta_seq != None:
                        fasta_seqs.append(fasta_seq)
                        if len(fasta_seqs) == len(selected):
                            break
                    temp = ""
                    counter = 1				
                temp += line
            elif line[0] != ">":
                temp += line
        fasta_seq = FilterSeq(temp, selected)
        if fasta_seq != None:
            if fasta_seq not in fasta_seqs:
                fasta_seqs.append(fasta_seq)

    return fasta_seqs

def read_fasta(fasta_path):
    """This function reads the fasta file from specified path
    and returns all the fastas as a list"""

    fasta_seqs = []

    with open(fasta_path, 'r') as file:
        counter = 0
        temp = ""
        for line in file:
            if line[0] == '>':
                counter = counter + 1
                if counter == 2:
                    fasta_seq = temp
                    if fasta_seq != "":
                        fasta_seqs.append(fasta_seq)
                    temp = ""
                    counter = 1				
                temp += line
            elif line[0] != ">":
                temp += line
        fasta_seq = temp
        if fasta_seq != "":
            fasta_seqs.append(fasta_seq)

    return fasta_seqs

def write_fasta(write_path, fasta_seqs):
    """This function takes a list of fastas and write them to 
    the specified path"""

    fastas = open(write_path, "w")

    for element in fasta_seqs:
        fastas.write(element)
    fastas.close()


def FilterSeq(body, substr):
    """check if a substring from a list of substrings
    exist in a string. If true, returns the string, else
    returns None"""
    cond = any(sub in body for sub in substr)
    if cond:
         return body
    else:
         return None

def Pruner(body, indexes):
    """This function takes the body of fasta and removes
    specified indexes"""
    indexes = set(indexes)
    #join indexes where it is not specified index
    return "".join([char for idx, char in enumerate(body) if idx not in indexes])

def fasta_pruner(fasta_seqs):
    """This function takes fastas as a list, removes indexes that are gaps for 
    all human proteins and returns a list of fasta sequences"""

    #get human sequences
    human_seqs = []

    for seq in fasta_seqs:
        if "HUMAN" in seq:
            human_seqs.append(seq)

    #create a dictionary where header and body of fasta are split
    list_of_fastas = []

    for seq in fasta_seqs:
        #get seperate lines
        split_seq = seq.split("\n")
        #get header (first line)
        header = split_seq[0]
        #combine the body
        body = "".join(map(str,split_seq[1:]))
        #assign header and body
        list_of_fastas.append({"header": header, "body": body})

    #create same dictionary for only human sequences
    list_of_humans = []

    for seq in human_seqs:
        split_seq = seq.split("\n")
        header = split_seq[0]
        body = "".join(map(str,split_seq[1:]))
        list_of_humans.append({"header": header, "body": body})

    #get indices of double gaps
    doublegapindices = []
    #15li
    if len(list_of_humans) == 15:
        for i in range(len(list_of_humans[0]["body"])):
            if (list_of_humans[0]["body"][i] == "-") & (list_of_humans[1]["body"][i] == "-") & (list_of_humans[2]["body"][i] == "-") & (list_of_humans[3]["body"][i] == "-") & (list_of_humans[4]["body"][i] == "-") & (list_of_humans[5]["body"][i] == "-") & (list_of_humans[6]["body"][i] == "-") & (list_of_humans[7]["body"][i] == "-") & (list_of_humans[8]["body"][i] == "-") & (list_of_humans[9]["body"][i] == "-") & (list_of_humans[10]["body"][i] == "-") & (list_of_humans[11]["body"][i] == "-") & (list_of_humans[12]["body"][i] == "-") & (list_of_humans[13]["body"][i] == "-") & (list_of_humans[14]["body"][i] == "-"):
                doublegapindices.append(i)
    
    elif len(list_of_humans) == 2:
        for i in range(len(list_of_humans[0]["body"])):
            if (list_of_humans[0]["body"][i] == "-") & (list_of_humans[1]["body"][i] == "-"):
                doublegapindices.append(i)
    #if there is only one human sequence in the list
    elif len(list_of_humans) == 1:
        for i in range(len(list_of_humans[0]["body"])):
            if (list_of_humans[0]["body"][i] == "-"):
                doublegapindices.append(i)

    #prune the body of fastas
    for i in range(len(list_of_fastas)):
        #get body
        temp = list_of_fastas[i]["body"]
        #remove double gap indices
        pruned_body = Pruner(temp, doublegapindices)
        #reassign pruned body to body
        list_of_fastas[i]["body"] = pruned_body

    #return the fastas to list state
    fasta_like = []

    for element in list_of_fastas:
        temp = element["header"] + "\n"
        counter = (-(-len(element["body"])//60))-1
        temp2 = element["body"]
        for i in range(counter):
            temp += (temp2[0:60] + "\n")
            temp3 = temp2[60:]
            temp2 = temp3
        temp += (temp2 + "\n")
        fasta_like.append(temp)

    return fasta_like

def fasta_splitter(fasta_seqs,path1, path2):
    """This function takes a list of fasta sequences and 2 csv files,
    as paths, then splits the list according to csv files"""
    
    colnames = ["sseqid", "stitle", "ox"]

    ortholog_list1 = pd.read_csv(path1, header = None, names = colnames)
    ortholog_list2 = pd.read_csv(path2, header = None, names = colnames)
    sseqids1 = ortholog_list1["sseqid"]
    sseqids2 = ortholog_list2["sseqid"]

    fasta_list1 = []

    for fasta in fasta_seqs:
        fasta_seq = FilterSeq(fasta, sseqids1)
        if fasta_seq != None:
            fasta_list1.append(fasta_seq)

    fasta_list2 = []

    for fasta in fasta_seqs:
        fasta_seq = FilterSeq(fasta, sseqids2)
        if fasta_seq != None:
            fasta_list2.append(fasta_seq)
    
    fasta_dictionary = {"fasta_list1":fasta_list1, "fasta_list2": fasta_list2}

    return fasta_dictionary

def fasta_to_dataframe(fasta_seqs):
    """This function takes a list of fasta sequences and
    turn them into a pandas dataframe where columns represent sequences and rows
    contain amino acid sequences at that position. It must be that every sequence
    in the list is in equal length"""

    list_of_fastas = []
    for seq in fasta_seqs:
        split_seq = seq.split("\n")
        header = split_seq[0]
        body = "".join(map(str,split_seq[1:]))
        body = list(body)
        list_of_fastas.append({"header": header, "body": body})

    df = pd.DataFrame()
    for fasta in list_of_fastas:
        df[fasta["header"]] = fasta["body"]

    return df

def Most_Common(lst):
    """this function returns the most common element in a list"""

    #count every element with counter
    data = Counter(lst)
    return data.most_common(1)[0][0]

def consensus(df, limit):
    """This function takes a pandas dataframe that has the multiple sequence
    alignment fasta data. It defines consensus as any amino acid that appears 
    above given limit in the position.It does not take gaps	into consideration.
    IMPORTANT: ADD A CONDITION FOR TOO MANY GAPS. For now if gaps are above 50%"""

    #create an empty dataframe to add consensus sequence with the frequency of the residues 
    consensus_seq = pd.DataFrame(columns = ["residue", "frequency", "status"])
    #get the dimensions of original dataframe to iterate over rows
    dims = df.shape

    for i in range(dims[0]):
        #get the values of the row (each row has the residue at that position for every protein)
        x = df.iloc[i].astype(str).tolist()
        #count the gaps
        gapcount = x.count("-")
        #count how many proteins there are
        poslength = len(x)
        #if more than 50% is gaps, there is no consensus
        if (gapcount/poslength) >= 0.5:
            consensus_seq = consensus_seq.append({"residue" : "-", "frequency" : 0.0, "status" : "NC"}, ignore_index = True)
        else:
            #get rid of the gaps
            y = [ elem for elem in x if elem != "-"]
            #get most common residue, possible consensus
            mostcommonres = Most_Common(y)
            #count how many times most common residue occurs
            mostcommonrescount = y.count(mostcommonres)
            #we take residues above 50% frequency so we can later use them (I CHANGED THE LIMIT TO 0.34 FROM 0.5) -- reverted back the change 
            if (mostcommonrescount/len(y)) > (0.5):
                #if most common residue is below limit it is not conserved (NC)
                if (mostcommonrescount/len(y)) < (limit):
                    consensus_seq = consensus_seq.append({"residue" : mostcommonres, "frequency" : mostcommonrescount/len(y), "status" : "NC"}, ignore_index = True)
                #if most common residue equal or above limit it is conserved
                else:
                    consensus_seq = consensus_seq.append({"residue" : mostcommonres, "frequency" : mostcommonrescount/len(y), "status" : "C"}, ignore_index = True)
            #if residues is below 50% it is not conserved
            else:
                consensus_seq = consensus_seq.append({"residue" : "-", "frequency": 0.0, "status" : "NC"}, ignore_index = True)

    return consensus_seq

def residuecontent(df):
    """Returns the most common 3 residues and their frequency for every position in the table of
    aligned proteins, gaps are considered"""
    #get the dimensions of original dataframe to iterate over rows
    dims = df.shape
    #create an empty dataframe to add consensus sequence with the frequency of the residues 
    contenttable = pd.DataFrame(columns = ["residue1", "frequency1", "residue2", "frequency2", "residue3", "frequency3"], index=list(range(dims[0])))

    for i in range(dims[0]):
        #get the values of the row (each row has the residue at that position for every protein)
        x = df.iloc[i].astype(str).tolist()
        residuefrequencies = {}
        for j in set(x):
            residuefrequencies[j] = x.count(j)/len(x)
        orderedfrequencies = {k: v for k, v in sorted(residuefrequencies.items(), key=lambda item: item[1], reverse=True)}
        first3pairs = {k: orderedfrequencies[k] for k in list(orderedfrequencies)[:3]}
        keylist = list(first3pairs.keys())
        for ii in range(len(keylist)):
            contenttable.at[i,"residue" + str(ii+1)] = keylist[ii]
            contenttable.at[i, "frequency" + str(ii+1)] = first3pairs[keylist[ii]]
            # contenttable = contenttable.append({"residue1" :keylist[0], "residue2" :keylist[1], "residue3" :keylist[2], "frequency1": first3pairs[keylist[0]], "frequency2": first3pairs[keylist[1]], "frequency3": first3pairs[keylist[2]]}, ignore_index=True) 
    return contenttable

def consensus_posadjusted(df_ort, df):
    """This function adjusts the positions of the residues from
    function"""

    #position adjustment
    #df_ort will include all ortholog sequences (return of fasta_to_dataframe function)
    #take the column names which include all sequence names
    cnames = list(df_ort.columns)
    cname = ""
    for element in cnames:
        #find human sequence
        if element.find("HUMAN") != -1:
            cname = element	
    #take human sequence and turn it into a list
    human_seq = df_ort[cname].tolist()
    indexes = []
    #take the indexes of gaps
    for i in range(len(human_seq)):
        if human_seq[i] == "-":
            indexes.append(i)
    #take the indexes we would like to keep
    indexes_to_keep = set(range(df.shape[0])) - set(indexes)
    df_posadjusted = df.take(list(indexes_to_keep))
    df_posadjusted.reset_index(drop = True, inplace = True)

    return df_posadjusted

def lookout(family, name, indexes):
    """This function looks up the positions specified in tables that 
    contain sequences of proteins of Class B1 GPCRs glucagon and secretin subfamily
    that are gathered from all orthologs msa, consensus residues from all orthologs msa and 
    frequencies of residues for every protein at every position. Takes subfamily names 
    (either "glucagon" or "secretin"), protein names, and positions to look for in that protein.
    Returns a dataframe with appropriate information. Protein name should be in all capitals.
    Possible names for glucagon: GIPR, GLP1R, GLP2R, GLR; for secretin: GHRHR, PACR, SCTR, VIPR1, VIPR2."""
    if family == "glucagon":
        glucagon_table = pd.read_csv("/cta/users/ofkonar/work/results/all_orthologs_pruned_consensus_withhuman.csv", header = 0)
        column = name + "_index"
        looked = glucagon_table.loc[glucagon_table[column].isin(indexes)]
    elif family == "secretin":
        secretin_table = pd.read_csv("/cta/users/ofkonar/work/results/secretin_subfamily/secretin_subfamily_orthologs_consensus_withhuman.csv", header = 0)
        column = name + "_index"
        looked = secretin_table.loc[secretin_table[column].isin(map(str, indexes))]		
    return looked

def special_residues(df1, df2):
    """This function takes two dataframe of consensus sequences with frequencies. It decides on the speciality of the residue
    based on 2 criteria by comparing df1 to df2: 
    1. df1 residue should not be a gap and
    2. df1 residue has to be different than df2	residue 
    !!!MAKE SURE THAT DATAFRAMES SHARE SAME DIMENSIONS AND HAD THE SAME LIMIT (consensus function limit)!!!
    NO RANKING IS DONE! Results are returned as a dataframe"""

    #read blosum80 table
    blosum80 = pd.read_csv("/cta/users/ofkonar/work/resources/blosum80.csv", index_col = 0)
    
    #create an empty dataframe to add the results
    special_residues = pd.DataFrame(columns = ["residue", "frequency", "blosum80_score"])

    #get the dimensions of the dataframes so we can iterate over rows
    dims = df1.shape[0]

    #non-specific residues are assigned blosum score of 100 so they are sorted at the bottom of the list(regular blosum ranges between -6 and 11) (I CHANGED IT TO "-")
    for i in range(dims):
        #the residue is a gap
        if df1.iloc[i]["residue"] == "-":
            special_residues = special_residues.append({"residue" : "-", "frequency" : 0.0, "blosum80_score": "-"}, ignore_index = True)
        #the residue is not frequent enough (not conserved)
        elif df1.iloc[i]["status"] == "NC":
            special_residues = special_residues.append({"residue" : "-", "frequency" : 0.0, "blosum80_score": "-"}, ignore_index = True)
        #residues are different
        elif df1.iloc[i]["residue"] != df2.iloc[i]["residue"]:
            residue1 = df1.loc[i]["residue"]
            residue2 = df2.loc[i]["residue"]
            #if second residue is a gap blosum score is given as -10 to indicate this situation (blosum scores range from -6 to 11) 
            if residue2 == "-":
                blosum80_score = -10
                special_residues = special_residues.append({"residue" : residue1, "frequency" : df1.iloc[i]["frequency"], "blosum80_score": blosum80_score}, ignore_index = True)
            #if second residue is not gap, blosum score is provided
            else:
                blosum80_score = blosum80.loc[residue1, residue2]
                special_residues = special_residues.append({"residue" : residue1, "frequency" : df1.iloc[i]["frequency"], "blosum80_score": blosum80_score}, ignore_index = True)
        #residues are the same so it is not specific
        else:
            special_residues = special_residues.append({"residue" : "-", "frequency" : 0.0, "blosum80_score": "-"}, ignore_index = True)
    
    return special_residues

def dendogram_seq(df_ort, special, newick):
    """This function will take ortholog data (df_ort, dataframe),
    special sequences data (special, dataframe) and newick data (string) 
    with clade info. It will add special sequence pos info to the 
    tree so we can see which residue resides within special positions
    of all orthologous proteins """

    #position adjustment, we get rid of positions of gaps on human sequence and all corresponding positions of orthologs
    #df_ort will include all ortholog sequences (return of fasta_to_dataframe function)
    #take the column names which include all sequence names
    cnames = list(df_ort.columns)
    cname = ""
    for element in cnames:
        #find human sequence
        if element.find("HUMAN") != -1:
            cname = element
    human_seq = df_ort[cname].tolist()
    indexes = []
    #take the indexes of gaps
    for i in range(len(human_seq)):
        if human_seq[i] == "-":
            indexes.append(i)
    #take the indexes we would like to keep
    indexes_to_keep = set(range(df_ort.shape[0])) - set(indexes)
    df_posadjusted = df_ort.take(list(indexes_to_keep))
    df_posadjusted.reset_index(drop = True, inplace = True)
    
    #get column names of special, residues should be under "residue column"
    special_residues = special["residue"].tolist()
    #get the positions of special residues
    idxs = []
    for res in special_residues:
        idx = int(''.join(x for x in res if x.isdigit())) - 1
        idxs.append(idx)
    #remove all non-special residues from dataframe
    idxs_to_keep = set(idxs)
    df_onlyspecific = df_posadjusted.take(list(idxs_to_keep))

    #get protein names and sequences
    names = []
    sequences = []
    for element in cnames:
        parts = element.split()
        names.append(parts[0].replace(">", ""))
        seq = df_onlyspecific[element].tolist()
        str1 = ""
        str2 = str1.join(seq)
        sequences.append(str2)

    for i in range(len(names)):
        new_label = sequences[i] + " " + names[i]
        newick = newick.replace(names[i], new_label)

    return newick

def finalize(df_ort, df, prt_name):
    """This function
    1. adjust the positions of residues to human sequence
    2. add Wootten numbering scheme to the data and
    3. eleminate not specific residues
    4. locate the residue on the domains (7TMs, ICLs, ECLs or N and C terminals)"""

    #position adjustment
    #df_ort will include all ortholog sequences (return of fasta_to_dataframe function)
    #take the column names which include all sequence names
    cnames = list(df_ort.columns)
    cname = ""
    for element in cnames:
        #find human sequence
        if element.find("HUMAN") != -1:
            cname = element	
    #take human sequence and turn it into a list
    human_seq = df_ort[cname].tolist()
    indexes = []
    #take the indexes of gaps
    for i in range(len(human_seq)):
        if human_seq[i] == "-":
            indexes.append(i)
    #take the indexes we would like to keep
    indexes_to_keep = set(range(df.shape[0])) - set(indexes)
    df_posadjusted = df.take(list(indexes_to_keep))
    df_posadjusted.reset_index(drop = True, inplace = True)

    #add Wootten numbers
    #Column names : 'CALCR', 'CALCRL', 'CRHR1', 'CRHR2', 'GHRHR', 'GIPR',
    #'GLP1R', 'GLP2R', 'GCGR', 'SCTR', 'PTH1R', 'PTH2R', 'PACR','VIPR1', 'VIPR2'
    Wootten = []
    name = []
    #read the Wootten file
    Wootten_table = pd.read_csv("/cta/users/ofkonar/work/resources/residue_table.csv")
    #get the required protein column
    my_column = Wootten_table[prt_name]
    for i in range(df_posadjusted.shape[0]):
        #get the residue name 
        res = df_posadjusted.iloc[i]["residue"] + str(i+1)
        name.append(res)
        #check if the residue has an Wootten number
        if res in set(my_column):
            #get the index of Wootten 
            windex = my_column[my_column == res].index[0]
            #append the Wootten numbering
            Wootten.append(str(Wootten_table.iloc[windex]["Wootten"]))
        #if residue has no Wootten numbering, append "-"
        else:
            Wootten.append("-")

    if prt_name == "GCGR":
        domains = ["N-ter"]*130 + ["TM1"]*36 + ["ICL1"]*4 + ["TM2"]*33 + ["ECL1"]*16 + ["TM3"]*36 + ["ICL2"]*6 + ["TM4"]*29 + ["ECL2"]*13 + ["TM5"]*32 + ["ICL3"]*6 + ["TM6"]*29 + ["ECL3"]*4 + ["TM7"]*29 + ["H8"]*26 + ["C-ter"]*48
    elif prt_name == "GIPR":
        domains = ["N-ter"]*126 + ["TM1"]*36 + ["ICL1"]*4 + ["TM2"]*28 + ["ECL1"]*19 + ["TM3"]*34 + ["ICL2"]*6 + ["TM4"]*29 + ["ECL2"]*13 + ["TM5"]*32 + ["ICL3"]*6 + ["TM6"]*29 + ["ECL3"]*4 + ["TM7"]*29 + ["H8"]*26 + ["C-ter"]*45
    elif prt_name == "GLP1R":
        domains = ["N-ter"]*137 + ["TM1"]*32 + ["ICL1"]*4 + ["TM2"]*28 + ["ECL1"]*22 + ["TM3"]*34 + ["ICL2"]*4 + ["TM4"]*31 + ["ECL2"]*9 + ["TM5"]*38 + ["ICL3"]*1 + ["TM6"]*31 + ["ECL3"]*5 + ["TM7"]*29 + ["H8"]*26 + ["C-ter"]*32
    elif prt_name == "GLP2R":
        domains = ["N-ter"]*167 + ["TM1"]*36 + ["ICL1"]*4 + ["TM2"]*28 + ["ECL1"]*22 + ["TM3"]*34 + ["ICL2"]*6 + ["TM4"]*29 + ["ECL2"]*13 + ["TM5"]*32 + ["ICL3"]*6 + ["TM6"]*29 + ["ECL3"]*4 + ["TM7"]*29 + ["H8"]*26 + ["C-ter"]*88
    elif prt_name == "GHRHR":
        domains = ["N-ter"]*115 + ["TM1"]*39 + ["ICL1"]*4 + ["TM2"]*29 + ["ECL1"]*11 + ["TM3"]*35 + ["ICL2"]*5 + ["TM4"]*30 + ["ECL2"]*13 + ["TM5"]*30 + ["ICL3"]*7 + ["TM6"]*30 + ["ECL3"]*4 + ["TM7"]*27 + ["H8"]*26 + ["C-ter"]*18
    elif prt_name == "PACR":
        domains = ["N-ter"]*138 + ["TM1"]*40 + ["ICL1"]*4 + ["TM2"]*28 + ["ECL1"]*12 + ["TM3"]*35 + ["ICL2"]*6 + ["TM4"]*29 + ["ECL2"]*9 + ["TM5"]*34 + ["ICL3"]*7 + ["TM6"]*30 + ["ECL3"]*5 + ["TM7"]*26 + ["H8"]*26 + ["C-ter"]*39
    elif prt_name == "SCTR":
        domains = ["N-ter"]*130 + ["TM1"]*37 + ["ICL1"]*4 + ["TM2"]*28 + ["ECL1"]*12 + ["TM3"]*35 + ["ICL2"]*3 + ["TM4"]*32 + ["ECL2"]*9 + ["TM5"]*36 + ["ICL3"]*5 + ["TM6"]*30 + ["ECL3"]*3 + ["TM7"]*27 + ["H8"]*26 + ["C-ter"]*23
    elif prt_name == "VIPR1":
        domains = ["N-ter"]*128 + ["TM1"]*39 + ["ICL1"]*4 + ["TM2"]*28 + ["ECL1"]*12 + ["TM3"]*35 + ["ICL2"]*3 + ["TM4"]*32 + ["ECL2"]*8 + ["TM5"]*34 + ["ICL3"]*7 + ["TM6"]*30 + ["ECL3"]*5 + ["TM7"]*26 + ["H8"]*26 + ["C-ter"]*40
    elif prt_name == "VIPR2":
        domains = ["N-ter"]*111 + ["TM1"]*40 + ["ICL1"]*4 + ["TM2"]*28 + ["ECL1"]*15 + ["TM3"]*35 + ["ICL2"]*5 + ["TM4"]*29 + ["ECL2"]*9 + ["TM5"]*34 + ["ICL3"]*7 + ["TM6"]*30 + ["ECL3"]*5 + ["TM7"]*26 + ["H8"]*26 + ["C-ter"]*34


    df_posadjusted.insert(3, "Wootten", Wootten, True)
    df_posadjusted.insert(4, "Name", name, True)
    df_posadjusted.insert(5, "domains", domains, True)
    

    #eleminate not specific residues
    df_onlyspecific = df_posadjusted[df_posadjusted.residue != "-"]
    df_onlyspecific.reset_index(drop = True, inplace = True)

    return df_onlyspecific 

def newick_to_beauty(string):
    """This function transforms newick to a more readable style"""

    new_string = list(string)
    beauty_string = []
    tab_counter = 0
    for i in range(len(new_string)):
        if new_string[i] == "(":
            beauty_string.append(new_string[i])
            beauty_string.append("\n")
            tab_counter += 1
            beauty_string += ["\t"]*tab_counter
        elif new_string[i] == ",":
            beauty_string.append(new_string[i])
            beauty_string.append("\n")
            beauty_string += ["\t"]*tab_counter
        elif new_string[i] == ")":
            if new_string[i-1] == ",":
                beauty_string.append("\n")
                beauty_string += ["\t"]*tab_counter
                beauty_string.append(new_string[i])
            else:
                beauty_string.append("\n")
                tab_counter -= 1
                beauty_string += ["\t"]*tab_counter
                beauty_string.append(new_string[i])

        else:
            beauty_string.append(new_string[i])

    return beauty_string

################################################################################################

################################################################################################

################################################################################################

#BURADAN AsAgISI UNDER CONSTRUCTION

################################################################################################

################################################################################################



aadata = {"aminoacid" : ["alanine", "arginine", "asparagine", "aspartic_acid", "asparagine_or_aspartic_acid", "cysteine", "glutamic_acid", "glutamine", "glutamine_or_glutamic_acid", "glycine", "histidine", "isoleucine", "leucine", "lysine", "methionine", "phenylalanine", "proline", "serine", "threonine", "tryptophan", "tyrosine", "valine"],
"three_letter_code" : ["ala", "arg", "asn", "asp", "asx", "cys", "glu", "gln", "glx", "gly", "his", "ile", "leu", "lys", "met", "phe", "pro", "ser", "thr", "trp", "tyr", "val"],
"one_letter_code" : ["A", "R", "N", "D", "B", "C", "E", "Q", "Z", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]}
aminoacid_table = pd.DataFrame.from_dict(aadata)
one_letter_code = ["A", "R", "N", "D", "B", "C", "E", "Q", "Z", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]