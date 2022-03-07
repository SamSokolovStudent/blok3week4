from aminodict import amino_counter, amino_name_dict
from collections import Counter
import os


def main():
    files = file_acceptor()
    file_opener(files)


def file_acceptor():
    """Accepts files from command line, returns all command line
    arguments as "files" in a list.

    :return file_list"""
    file_extension = ".FASTA"
    file_list = []
    for file in os.listdir():
        if file.endswith(file_extension):
            file_list.append(file)
    return file_list


def file_opener(file_list):
    """Opens file in current scope using with to automatically close
    when processing on file is done."""
    for file in file_list:
        # Opens file in current scope.
        with open(file, 'r') as open_file:
            sequence_dict = fasta_reader(open_file)
            envelop_dict, internal_dict = dict_sorter(sequence_dict)
            data_output(file, envelop_dict, internal_dict)


def fasta_reader(open_file):
    """Reads all the lines in a file and makes a structured dictionary
    out of them that uses the header of the file as key and the sequence
    of that protein as the value.

    :param open_file
    :return sequence_dict"""
    # Make empty header and dictionary that can be used to store all
    # sequences in an organized manner.
    header = ""
    sequence_dict = {}
    for line in open_file:
        if line.startswith(">"):
            # Makes a list of all elements of the header line split by
            # spaces. Then gets the 1st / 0th element from that list.
            header = line.split(" ")
            header = header[0], header[1]
            header = ''.join(map(str, header))
            # Makes a key out of the header with an empty value.
            sequence_dict[header] = ""
        else:
            # Removes newline character at end of sentence and adds the
            # correct sequence as value to the dictionary.
            sequence_dict[header] += line.rstrip()
    return sequence_dict


def dict_sorter(sequence_dict):
    """"Imports a protein sequence dictionary and determines if the
    protein is an envelope protein by the code that is appended to the
    string. Returns two dictionaries with all envelope proteins and
    all internal proteins.

    :param sequence_dict
    :return env_dict, internal_dict"""
    env_dict = {}
    internal_dict = {}
    for key, value in sequence_dict.items():
        if "=env" in key:
            env_dict[key] = value
        else:
            internal_dict[key] = value
    return env_dict, internal_dict


def data_output(file_name, envelop, internal):
    """"Outputs the data from all envelope and internal proteins per
    file. Loops over dictionaries and gives output per key in order to
    ensure all data is covered.

    :param file_name
    :param envelop
    :param internal"""
    # Prints in bold text for easier oversight.
    print('\033[1m' + f"Output for {file_name}." + '\033[0m')
    if envelop:
        print("Envelop Protein Data")
        print("-" * 30)
        for key, value in envelop.items():
            data_processor(key, value)
    print("Internal Protein Data\n")
    print("-" * 30)
    for key, value in internal.items():
        data_processor(key, value)
    print("\n"*2)


def data_processor(key, value):
    """"Takes in sequence dictionaries from previous loop and outputs
    desired data per key. Prints to console in formatted manner and
    resets global import dict when done.

    :param key
    :param value"""
    # For every letter in sequence add a point to specific letter in
    # global imported dict.
    for s in value:
        amino_counter[s] += 1
    # Averages the data inside global dict.
    for s in amino_counter:
        amino_counter[s] = round(amino_counter[s] / len(value) * 100, 2)
    # Sets cysteine and tryptophan counts.
    cys_count = amino_counter["C"]
    trypt_count = amino_counter["W"]
    # Sets hydrophobic and hydrophilic counts.
    hydrophobic = (
        amino_counter["F"], amino_counter["W"], amino_counter["I"],
        amino_counter["L"], amino_counter["M"], amino_counter["V"],
        amino_counter["A"], amino_counter["C"])
    hydrophilic = (
        amino_counter["K"], amino_counter["R"], amino_counter["E"],
        amino_counter["D"], amino_counter["Q"], amino_counter["N"])
    # Counts amounts of amino acids in dict and sets most/least common.
    counted_dict = Counter(amino_counter)
    most_common = counted_dict.most_common(3)
    least_common = counted_dict.most_common()[:-3-1:-1]
    print('\033[1m' + key + '\033[0m' + '\n')
    print(f"Percentage Cysteine: {cys_count}%")
    print(f"Percentage Tryptophan: {trypt_count}%\n")
    print(f"Hydrophobic amino acids: {round(sum(hydrophobic, 2))}%")
    print(f"Hydrophilic amino acids: {round(sum(hydrophilic, 2))}%\n")
    print(f"Most common amino acids")
    for element in most_common:
        print(f"{amino_name_dict[element[0]]} : {element[1]}%")
    print(f"\nLeast common amino acids")
    for element in least_common:
        print(f"{amino_name_dict[element[0]]} : {element[1]}%")
    print('\n')
    # Resets global dictionary.
    amino_counter.update({}.fromkeys(amino_counter, 0))


if __name__ == '__main__':
    main()
