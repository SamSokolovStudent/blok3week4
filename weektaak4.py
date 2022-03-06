from aminodict import amino_counter
import sys
from collections import Counter


def main():
    files = file_acceptor()
    file_opener(files)


def file_acceptor():
    """Accepts files from command line, returns all command line
    arguments as "files" in a list.

    :return file_list"""
    # All files from position 1 and on are put inside a list.
    file_list = sys.argv[1:]
    file_1 = "fumarasehydratase.FASTA"
    file_2 = "glucosetransporter.FASTA"
    file_3 = "gpcr.FASTA"
    file_4 = "receptorkinase.FASTA"
    file_5 = "hiv1cdsprotein.FASTA"
    file_6 = "hiv2cdsprotein.FASTA"
    file_7 = "sivcdsprotein.FASTA"
    file_8 = "sivmnd2cdsprotein.FASTA"
    file_list = [file_1, file_2, file_3, file_4,
                 file_5, file_6, file_7, file_8]
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
    # sequences in an oranized manner.
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
    env = "=env"
    env_dict = {}
    internal_dict = {}
    for key, value in sequence_dict.items():
        if env in key:
            env_dict[key] = value
        else:
            internal_dict[key] = value
    return env_dict, internal_dict


def data_output(file_name, envelop, internal):
    print(f"Output for {file_name}.")
    if envelop:
        print("Envelop Protein Data\n", "-" * 30)
        for value in envelop.values():
            data_processor(value)
    print("Internal Protein Data\n", "-" * 20)
    for value in internal.values():
        data_processor(value)
    print("\n"*2)


def data_processor(value):
    for s in value:
        amino_counter[s] += 1
    for s in amino_counter:
        amino_counter[s] = round(amino_counter[s] / len(value) * 100, 2)
    cys_count = amino_counter["C"]
    trypt_count = amino_counter["W"]
    hydrophobic = (
        amino_counter["F"], amino_counter["W"], amino_counter["I"],
        amino_counter["L"], amino_counter["M"], amino_counter["V"],
        amino_counter["A"], amino_counter["C"])
    hydrophilic = (
        amino_counter["K"], amino_counter["R"], amino_counter["E"],
        amino_counter["D"], amino_counter["Q"], amino_counter["N"])
    counted_dict = Counter(amino_counter)
    most_common = counted_dict.most_common(3)
    least_common = counted_dict.most_common()[:-3-1:-1]
    print(f"Percentage Cysteine: {cys_count}%")
    print(f"Percentage Tryptophan: {trypt_count}%")
    print(f"Hydrophobic amino acids: {round(sum(hydrophobic, 2))}%")
    print(f"Hydrophilic amino acids: {round(sum(hydrophilic, 2))}%")
    print(f"Most common amino acids")
    for keys in most_common:
        print(f"{keys[0]} : {keys[1]}%")
    print(f"Least common amino acids")
    for keys in least_common:
        print(f"{keys[0]} : {keys[1]}%")
    amino_counter.update({}.fromkeys(amino_counter, 0))


if __name__ == '__main__':
    main()
