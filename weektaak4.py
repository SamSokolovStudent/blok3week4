from codondict import codons
import sys


class ProteinData:
    def __init__(self, data_dict):
        self._file_name = self.file_name_setter(data_dict[""])
        self._cys_value = self.cys_setter()


    def file_name_setter(self, file_name):
        return self._file_name

    def cys_setter(self, cys_value):
        _cys_value = cys_value
        return _cys_value

    def trp_setter(self, trp_value):
        _trp_value = trp_value
        return _trp_value

    def hydrophobe_setter(self), :

    def hydrophile_setter(self):

    def most_common(self):

    def least_common(self):


def main():
    files = file_acceptor()
    file_opener(files)


def file_acceptor():
    """Accepts files from command line, returns all command line
    arguments as "files" in a list.

    :return file_list"""
    # All files from position 1 and on are put inside a list.
    file_list = sys.argv[1:]
    return file_list


def file_opener(file_list):
    """Opens file in current scope using with to automatically close
    when processing on file is done."""
    for file in file_list:
        # Opens file in current scope.
        with open(file, 'r') as open_file:
            sequence_dict = fasta_reader(open_file)
            split_sequence_dict = codon_splitter(sequence_dict)


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
            header = header[0]
            # Makes a key out of the header with an empty value.
            sequence_dict[header] = ""
        else:
            # Removes newline character at end of sentence and adds the
            # correct sequence as value to the dictionary.
            sequence_dict[header] += line.rstrip()
    return sequence_dict


def codon_splitter(sequence_dict):
    """Takes in the dictionary made in fasta_reader() and splits the
    sequences in segments of 3. Then creates a new dictionary using
    the same header and sets a list of all split codons as the value
    for this dictionary.

    :param sequence_dict
    :return splitting_dictionary"""
    splitting_dictionary = {}
    for header, sequence in sequence_dict.items():
        current_sequence = sequence
        # n is the constant which determines the amount of bases the
        # codons are split by
        n = 3
        # Splits sequence in segments of 3 and makes a list.
        split_codons = [sequence[i:i+n] for i in range(0, len(sequence), n)]
        splitting_dictionary[header] = split_codons
    return splitting_dictionary


if __name__ == '__main__':
    main()
