from Bio.Seq import translate
from typing import Tuple
import logging
import requests
import json
import sys
import re
import argparse

#Contains server addresses as values and the name of the
#databases as keys. Set as global variable since the these
#values are used constantly throughout the script
SERVER = {
    "mygene":"https://mygene.info/v3/",
    "ensembl":"https://rest.ensembl.org/"
    }

#Contains server specific endpoints along with content-type interactions 
#if applicable as values within in a list. If a content-type is not needed
#for the requests.get() to complete the task then it is an empty string 
#the name of the search focus are the keys. Set as global variable since the these
#values are used constantly throughout the script
#NOTE all query templates have the species set to human.
SEARCH_TYPE = {"gene":["query?q={}&species=human",""],
                "entrez_id":["gene/{}?species=human",""],
                "sequence":["sequence/id/{}","text/plain"],
                "homology":["homology/id/human/{}","application/json"]
                }


def fetch_endpoint_POST(server_address:str, content_type:str)->requests.get:
    #takes in a string with the the full web url and the endpoints
    #returns the Response object for further manipulation
    logging.info(f"Searching {server_address}")
    results = requests.get(server_address, headers={"Content-Type":content_type})
    #checks to see if the results address is ok which is web code 200
    if not results.ok:
        #if the code is not 200, the script exists.
        logging.warning(f"Something went wrong while fetching the results!!!")
        results.raise_for_status()
        sys.exit(1)
    return results

def convert_to_json(request_results:requests.get)->json:
    #simple one line function that takes in a requests object
    #and converts it to a json.
    return request_results.json()

def create_query_string(query_item:str, server_address:dict, server:str, search_type:dict, search_by:str)->Tuple[str,str]:
    #this function is used to weave togther the URL string and
    #returns what the content type should be.
    #a server url is selected based on the server name
    server_to_search = server_address[server]
    #the query template and output content are extracted from the search type dict
    query, output_content = search_type[search_by]
    #the query string is updated to include the item that is being queried
    query = query.format(query_item)
    #the final output could look like this if MC1R, mygene, and gene are
    #the query_item, server, and by terms
    #"https://mygene.info/v3/query?q=MC1R&species=human
    server_query = server_to_search + query
    return server_query, output_content


def mygene_search(gene_name:str, search_by:str="gene", search_type:dict = SEARCH_TYPE, server:str="mygene", server_address:dict = SERVER)->str:
    #this function takes in the gene name (the only required argument) as a string 
    #the necesseary server url and end points are extrapolated using the default values associated with
    #search_by and server
    #the output is a json from the MyGene database containing information such as the entrez ID
    logging.info("Constructing search string...")
    server_query, output_content = create_query_string(gene_name, server_address,
                                    server, search_type, search_by)
    logging.info(f"Looking for {gene_name} within MyGene database.")
    results = fetch_endpoint_POST(server_query, output_content)
    json_results = convert_to_json(results)
    return json_results

def entrez_to_ensembl(gene_name:str ,mygene_json:json, search_by:str = "entrez_id", search_type:dict = SEARCH_TYPE, server:str="mygene", server_address:dict = SERVER)->str:
    #this function takes in the gene name (required argument) as a string and 
    #the associated json generated from the MyGene database (second requried argument)
    #the necesseary server url and end points are extrapolated using the default values associated with
    #search_by and server
    #the output is a string reflecting the ensembl ID
    logging.info("Extracting information from provided JSON")
    mygene_results = mygene_json['hits']
    try:
        mygene_result = [hit for hit in mygene_results if hit['symbol'] == gene_name][0]
    except IndexError as error:
        logging.error(f"No hits were fround for {gene_name}")
        logging.warning(f"Result JSON from MyGene Database: {mygene_json}")
        logging.warning(f"Exiting Workflow")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An error was found! Exiting!")
        logging.error(e)
        sys.exit(1)
    entrez_id= mygene_result['entrezgene']
    logging.info(f"Found the following Entrenz ID: {entrez_id}")
    server_query, output_content = create_query_string(entrez_id, server_address, server, search_type, search_by)
    results = fetch_endpoint_POST(server_query, output_content)
    json_results = convert_to_json(results)
    ensembl_id = json_results['ensembl']['gene']
    logging.info(f"Found the following Ensembl ID: {ensembl_id}")
    return ensembl_id


def find_ensemble_sequence(en_id:str, search_by:str = "sequence", search_type:dict=SEARCH_TYPE, server:str = "ensembl", server_address:dict = SERVER)->str:
    #this function takes in the ensembl id (the only required argument) as a string
    #the necesseary server url and end points are extrapolated using the default values associated with
    #search_by and server
    #the output is the DNA sequence associated with the gene
    logging.info(f"Searching for DNA sequence using Ensembl ID: {en_id}")
    server_query, output_content = create_query_string(en_id, server_address, server, search_type, search_by)
    results = fetch_endpoint_POST(server_query, output_content)
    dna_sequence = results.text
    return str(dna_sequence)

def write_to_file(content:str, filename:str, writemode:str="write")->None:
    #this function outputs a text file using variabel content for the internal
    #content of the text file, filename will determine the name and extension of the
    #output file, and writemode by default is write, however the user can append
    #content to a file by using append.
    mode = {"write":"w", "append":"a"}
    #dictionary mode is used what allows this function to be flexible
    logging.info(f"Using mode {writemode} on {filename}")
    try:
        #this try and except block is used to catch an issues that arise
        #when attempting to write or append to the file.
        with open(filename, mode[writemode]) as f:
            f.write(content)
        #if the operation is successful the user will be informed, otherwise
        #the exception will be caught and presented to the user. The script
        #will exit out with error code 1
        logging.info(f"{writemode[0].upper()}{writemode[1:]} mode on {filename} successful")
    except Exception as e:
        logging.warning("An unexpected error, exiting the workflow!!!")
        logging.error(e)
        sys.exit(1)


def find_longest_orf(dna_sequence:str)->str:
    #this functing takes in a dna sequence as a string
    #and outputs the longest orf as a string
    orfs = []
    #an empty list orfs is initialized and the start
    #and stop codons are initialized as strings
    start_codon = "ATG"
    #the stop codoon has pipe | included for use with 
    #regex's find iter function.
    stop_codon = "TAA|TGA|TAG"
    if start_codon in dna_sequence:
        #if our stop codon exists at all in the string dna_sequence
        #we itterate through each match found by finditer()
        for starting_match in re.finditer(start_codon,dna_sequence):
            #the starting index is extracted using the start() method on 
            #our match object
            remaining_sequence = dna_sequence[starting_match.start():]
            #the same style of itteration and use of re.finditer is used
            for stopping_match in re.finditer(stop_codon,remaining_sequence):
                #this stop match should be complimentary to the start match from the start
                potential_orf = remaining_sequence[:stopping_match.end()]
                if len(potential_orf) % 3 == 0:
                    #we makes sure that our sequence is in frame
                    #if not we will not see an asterisk (*) when we translate
                    #the DNA to Amino Acid space indicating an inccorect translation
                    orfs.append(potential_orf)
                    #the orf is appended, and we no longer need to look at the 
                    #remaining sequence using our first start codon as the refrence
                    #and we move onto the next one.
                    break
    #the list orf is sorted based on the length of the strings
    #and is sorted in descending order
    orfs.sort(key=len, reverse=True)
    #implying that our largest orf should be the 0th index
    found_longest_orf = orfs[0]
    return found_longest_orf


def dna_to_aa(dna_sequence:str)->str:
    #this function takes in a DNA string and uses
    #Bio pythons translate function to convert the
    #DNA sequence to an Amino Acid sequence
    aa_sequence = translate(dna_sequence)
    return aa_sequence


def find_homology(en_id:str, search_by:str = "homology", search_type:dict=SEARCH_TYPE, 
                  server:str = "ensembl", server_address:dict = SERVER)->json:
    #this function takes in the ensemble ID, and uses a default value to extract the end point
    #for the ensembl database.
    logging.info(f"Searching homologs for: {en_id}")
    #the url with the endpoint and content output is established
    query_string, output_content = create_query_string(en_id, server_address, server, search_type, search_by)
    #the results are a results object that is converted into a JSON for parsing
    results = fetch_endpoint_POST(query_string, output_content)
    json_results = convert_to_json(results)
    logging.info("Completed homolog search")
    return json_results

def sort_homology(homology_json:json)->str:
    #this function takes in a json containing homology results
    #from the ensembl database.
    logging.info("Sorting homology results...")
    #an empty lists is created, and we set homology results to
    #be associated with the value of key homologies, as it contains
    #a list of the resulting homologies
    homology_species = []
    homology_results = homology_json['data'][0]["homologies"]
    logging.info(f"Number of homologs found: {len(homology_results)}")
    #we itterate through the list of homology results
    for result in homology_results:
        species = result['target']['species']
        #before doing anything we check to make sure what the target
        #species is.
        if species == "homo_sapiens":
            #if the target species is homo_sapien it is skipped 
            #and the for loop continues
            continue
        #if the species is not homo sapiens we append the resulting species
        #to our instantiated list homology_species.
        homology_species.append(species+"\n")
    logging.info(f"Number of non-unique species found: {len(homology_species)}")
    #a set is generated from the list of homology_species, as this is a quick
    #way to eliminate duplicates within the list. The set is converted back into
    #a list in order to use join to create a long string. That will be used
    #to generatre a text file.
    unique_species = list(set(homology_species))
    logging.info(f"Number of unique species found: {len(unique_species)}")
    unique_species = ''.join(unique_species)
    return unique_species

def arg_parser()-> argparse.ArgumentParser:
    #this function encapsulates argparse arguments 
    #this line describes what the program does
    parser = argparse.ArgumentParser(prog='Gene2FASTA',
                                    description='Outputs a FASTA file based on a given gene name using MyGene and Ensembl databases via REST API',
                                    epilog='Example use case: python3 main.py "TTN"')
    #add arguments adds additional inputs
    #gene name is a required positional argument
    parser.add_argument('gene_name', type=str, help='gene name of interest ex: HER2, TTN, MC1R')
    #--logging is an optional argument and not calling it keeps logging off.
    parser.add_argument('-l', '--logging', action='store_true', help='enable verbose logging.')
    #the parser object is returned
    return parser

def main(arguments:argparse.ArgumentParser)->None:
    #main orchestrates the fucntions above and uses
    #the values passed into the command line to direct
    #the workflow
    set_logging = args.logging
    gene_name = args.gene_name.upper()
    #if --logging was used in the command line,
    #the logging basic config is kicked off and enables
    #logging at the INFO level.
    if args.logging is True:
        logging.basicConfig(level=logging.INFO, 
        format='[%(levelname)s]-%(asctime)s:::%(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
    #a gene sarch using MyGene database is kicked off.
    output_json = mygene_search(gene_name)
    #using the outputted json, the ensembl ID is generated using
    #the entrez ID from the json.
    en_out = entrez_to_ensembl(gene_name, output_json)
    #the ensembl ID is used to find the sequence from
    #the ensembl database
    seq_out = find_ensemble_sequence(en_out)
    #a fasta header is generated using the gene name
    #ensembl ID and the sequence is added on after a new line
    fasta = f">{gene_name} {en_out}\n"+seq_out
    #this very long strong is used to generate the fasta file
    write_to_file(fasta,f"{gene_name.lower()}.fasta")
    #the dna sequence found from the ensembl database
    longest_orf = find_longest_orf(seq_out)
    #the amino acid sequence is created using the longest orf
    aa_seq = dna_to_aa(longest_orf)
    #a fasta header for the protein sequence is created 
    #using the gene name with the phrase "longest Orf" added on 
    aa_header = f"\n>{gene_name.lower()} Longest Orf\n"
    #the amino acid header and sequence are concatinated into
    #one large string which is passed to the write file function
    aa_fasta = aa_header + aa_seq
    #append is used so that it can be added to our original DNA
    #fasta
    write_to_file(aa_fasta,f"{gene_name.lower()}.fasta","append")
    #using the ensembl ID the homologs are found
    found_homologs = find_homology(en_out)
    #the homologs of unique species are kept
    unique_species = sort_homology(found_homologs)
    #the unique species are written to a text file with *_homology_list.txt
    write_to_file(unique_species, f"{gene_name.lower()}_homology_list.txt")
    logging.info(f"Job complete for gene: {gene_name}")

if __name__ == "__main__":
    parser = arg_parser()
    args = parser.parse_args()
    main(args)