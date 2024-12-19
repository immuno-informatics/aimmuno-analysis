# %%
import numpy as np
import random
from unicodedata import decimal
import pandas as pd
import copy
import argparse
import logging
from sklearn.model_selection import train_test_split


logging.basicConfig(level=logging.DEBUG)
logging.debug("This will get logged")
logging.disable()
random.seed(42)
np.random.seed(42)

parser = argparse.ArgumentParser(
    description="Give the minimum and maxim lenghth of peptide and MHC moelcule"
)

parser.add_argument(
    "--sample_num",
    type=int,
    default=100,
    required=True,
    help="number of true samples to generate",
)
parser.add_argument( "--out_dir", type=str, required=True, help="Output directory")
parser.add_argument(
    "--decoys", type=int, default=99, required=False, help="number of decoys per sample"
)
parser.add_argument(
    "--peptide_min", type=int, default=7, required=False, help="peptide maximum length"
)
parser.add_argument(
    "--peptide_max", type=int, default=15, required=False, help="peptide maximum length"
)
parser.add_argument(
    "--mhc_min", type=int, default=15, required=False, help="MHC minimum length"
)
parser.add_argument(
    "--mhc_max", type=int, default=35, required=False, help="MHC maximum length"
)
parser.add_argument(
    "--random_decoy",
    action="store_true",
    default=False,
    help="If the decoy should be completely random",
)
# parser.add_argument(
#     "--num_letter_decoy",
#     default=random,
#     type=str,
#     help="number of letters to change in the decoy",
# )
args = parser.parse_args()


amino_acids = "G A L M F W K Q E S P V I C Y H R N D T"
amino_acid_letter = amino_acids.split()
amino_acid_letter = [
    "A",
    "R",
    "N",
    "D",
    "C",
    "Q",
    "E",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
]


class AminoAcids:
    def __init__(self, letters):
        """
        a class to create a dictionary of amino acid afinity from a given matrix
        """
        self.lettrs = letters
        self.binding_matrix = np.load("observed_expected_ratio.npy")

    def _matrix_to_dic(self):
        """
        Converts a matrix to a dictionary.

        Returns:
            dict: A dictionary where the keys are tuples representing the indices of the matrix,
                  and the values are the corresponding matrix values.
        """
        mat_to_dict = {}

        for a_acid_1, value_1 in zip(self.lettrs, self.binding_matrix):
            for a_acid_2, value_2 in zip(self.lettrs, value_1):
                mat_to_dict[(a_acid_1, a_acid_2)] = value_2

        return {
            k: v
            for k, v in sorted(
                mat_to_dict.items(), key=lambda item: item[1], reverse=True
            )
        }

    def k_pair(self, k=40):
        """
        Returns a list of the first `k` keys from a sorted dictionary.

        Parameters:
            k (int): The number of keys to return. Defaults to 40.

        Returns:
            list: A list of the first `k` keys(tuples) from the sorted dictionary.
        """
        sort_mat_to_dict = self._matrix_to_dic()
        list_sorted_dict = [i for i in sort_mat_to_dict.keys()]
        return list_sorted_dict[:k]


class PeptideMhcPair:

    def __init__(self, pairs):
        self.pairs = pairs
         # number of places to change the letter in the decoy. "random" or 1, 2, 3

    @staticmethod
    def random_letter_maker(N):
        """
        Generates random letter sequence length N
        """
        return "".join(random.choice(amino_acid_letter) for _ in range(N))

    def _peptide_mhc_sequence(self, peptide_part_len, binding_pairs, index_selection):
        """
        A function made for creating peptids OR MHC sequence
        args:
            peptide_part_len("tupl or list"): a tuple or list indiciating the number of random
            letters between the binding pairs e.g.: (2,4,3,4) ==> aa b aaaa b aaa b aaaa
            binding_pairs(list of tuples): a list of tuples of binding amino acid pair letter e.g. [("A", "G"), ("Y", "E"), ("Q", "N")]
            index_selection(int): a radnom number between 0,1 to select the direction of binding from amino acid pairs e.g. ("A", "G") ==> 0 =="A" , 1 == "G"
        returns:
            peptide_str: a string represting a peptide(or MHC)
            sequences_list_peptide: a list of the parts of the peptide(or MHC)
        """

        # print(binding_pairs)

        first_part_peptide = PeptideMhcPair.random_letter_maker(peptide_part_len[0])
        second_part_peptide = PeptideMhcPair.random_letter_maker(peptide_part_len[1])
        third_part_peptide = PeptideMhcPair.random_letter_maker(peptide_part_len[2])
        fourth_part_peptide = PeptideMhcPair.random_letter_maker(peptide_part_len[3])

        first_amino_peptide = binding_pairs[0][index_selection]
        second_amino_peptide = binding_pairs[1][index_selection]
        third_amino_peptide = binding_pairs[2][index_selection]

        sequences_list_peptide = [
            first_part_peptide,
            first_amino_peptide,
            second_part_peptide,
            second_amino_peptide,
            third_part_peptide,
            third_amino_peptide,
            fourth_part_peptide,
        ]

        peptide_str = "".join(sequences_list_peptide)

        return peptide_str, sequences_list_peptide

    def flank_maker(self, len_flank):
        """
        Makes random string of size len_flank
        """
        return PeptideMhcPair.random_letter_maker(len_flank)

    def peptide_maker(
        self,
        peptide_part_len=(2, 4, 2, 2),
        mhc_part_len=(4, 3, 3, 20),
        n_flank_len=15,
        c_flank_len=15,
        random_choice=True,
        num_of_pairs=3,
    ):
        """
        This function makes peptied and its corresponding mhc protein.

        :param amino_pairs: a list of tuples of pairing amino acids e.g. [("A","G"),("Y","E"),("Q","N")]
        :param peptide_part_len: a tuple of peptide length around each amino acid first_seq|peptide1|second_seq|peptide2|third_seq|peptide3|foruth_seq
        :param mhc_part_len: a tuple of peptide length around each amino acid first_seq|peptide1|second_seq|peptide2|third_seq|peptide3|foruth_seq
        :param repeat: if repate is true, two matching seq would occure in mhc (NOT IN THIS VERSION)
        :param random_choice: if True, the order of binding letters tuples will be random
        :param num_of_pairs: number of binding pairs

        :returns: amindo acid sequence string, amindo acid sequence list beofre joining, peptied sequence string, peptide sequence list beofre joining
        """
        amino_pairs = self.pairs  # [("A","G"),("Y","E"),...,("Q","N")]

        if random_choice:
            binding_pairs = random.sample(
                amino_pairs, k=num_of_pairs
            )  # choosing pairs randomly (3 pairs)

        else:
            binding_pairs = amino_pairs

        index_selection = random.choice([0, 1])
        pair_selection_mhc = 0 if index_selection == 1 else 1
        # print(index_selection,pair_selection_mhc)
        # random_bind_seq = random.choice([0, 1])

        # n_flank and c_flank
        n_flank = self.flank_maker(n_flank_len)  # Not using for now
        c_flank = self.flank_maker(c_flank_len)  # Not using for now

        peptidseq_str, peptidseq_list = self._peptide_mhc_sequence(
            peptide_part_len, binding_pairs, index_selection
        )
        Mhcseq_str, Mhcseq_list = self._peptide_mhc_sequence(
            mhc_part_len, binding_pairs, pair_selection_mhc
        )

        return peptidseq_str, peptidseq_list, Mhcseq_str, Mhcseq_list

    def _peptide_decoy_maker(
        self, a_list_main_seq, a_list_second_seq, random_letter
    ):
        """
        args:
            a_list_main_seq(list): a list of splitted peptides and amino acids, index 1,3,5 are the spots of binding
            random_letter(list, str): if it is a list it cho
        returns:
            seq(str): altered version of the originial sequence
        """
        # print("random_letter    ", random_letter)
        if random_letter == "random":
            number_of_letters = random.choice([1, 2, 3])
        elif int(random_letter) == 1:
            number_of_letters = 1
        elif int(random_letter) == 2:
            number_of_letters = 2
        elif int(random_letter) == 3:
            number_of_letters = 3

        a_list = copy.deepcopy(a_list_main_seq)
        counter_part_list = a_list_second_seq

        while True:
            if number_of_letters == 1:

                n = random.choice([1, 3, 5])
                d = random.choice(range(len(amino_acid_letter)))

                # if the letter pair and its pairs are not in the list of pairing amino acids
                if (amino_acid_letter[d], counter_part_list[n][0]) not in self.pairs:
                    # if a_list[n][0] != amino_acid_letter[d]:  # not to choose the same letter: Does not work (15/2)
                    replaced_text = a_list[n].replace(
                        a_list[n][0], amino_acid_letter[d]
                    )
                    a_list[n] = replaced_text

            elif number_of_letters == 2:
                n = random.choices([1, 3, 5], k=2)
                d_1 = random.choice(range(len(amino_acid_letter)))
                d_2 = random.choice(range(len(amino_acid_letter)))

                # not to choose the same letter
                # if a_list[n[0]][0] != amino_acids[d_1] and a_list[n[1]][0] != amino_acids[d_2]:
                if (
                    amino_acid_letter[d_1],
                    counter_part_list[n[0]],
                ) not in self.pairs and (
                    amino_acid_letter[d_2],
                    counter_part_list[n[1]],
                ) not in self.pairs:
                    # print("I am here 2???")
                    replaced_text_1 = a_list[n[0]].replace(
                        a_list[n[0]][0], amino_acid_letter[d_1]
                    )
                    replaced_text_2 = a_list[n[1]].replace(
                        a_list[n[1]][0], amino_acid_letter[d_2]
                    )
                    a_list[n[0]] = replaced_text_1
                    a_list[n[1]] = replaced_text_2

            elif number_of_letters == 3:
                n = [1, 3, 5]
                d_1 = random.choice(range(len(amino_acid_letter)))
                d_2 = random.choice(range(len(amino_acid_letter)))
                d_3 = random.choice(range(len(amino_acid_letter)))
                # not to choose the same letter
                if (
                    (amino_acid_letter[d_1], counter_part_list[n[0]]) not in self.pairs
                    and (amino_acid_letter[d_2], counter_part_list[n[1]])
                    not in self.pairs
                    and (amino_acid_letter[d_3], counter_part_list[n[2]])
                    not in self.pairs
                ):

                    replaced_text_1 = a_list[n[0]].replace(
                        a_list[n[0]][0], amino_acid_letter[d_1]
                    )
                    replaced_text_2 = a_list[n[1]].replace(
                        a_list[n[1]][0], amino_acid_letter[d_2]
                    )
                    replaced_text_3 = a_list[n[2]].replace(
                        a_list[n[2]][0], amino_acid_letter[d_2]
                    )
                    a_list[n[0]] = replaced_text_1
                    a_list[n[1]] = replaced_text_2
                    a_list[n[2]] = replaced_text_3

            if "".join(a_list_main_seq) != "".join(a_list):
                break

        seq = "".join(a_list)
        seq = seq.replace(" ", "")
        return seq, a_list

    def decoy_output(
        self,
        peptidseq_str,
        peptidseq_list,
        Mhcseq_str,
        Mhcseq_list,
        decoy_count=1,
        rand_let="random",
    ):
        """

        return:
            temp_set(list): a listed set of (peptide_decoy, mhc_decoy,"False"); "false" is the label of decoys
        """

        temp_set = set()

        while len(temp_set) < decoy_count:

            peptide_or_mhc_choice = random.choice(
                ["peptidseq_list", "Mhcseq_list", "both"]
            )
            peptide_decoy = None
            if peptide_or_mhc_choice == "peptidseq_list":
                # print("making peptide decoy")
                logging.info("making peptide decoy")
                peptide_decoy, _ = self._peptide_decoy_maker(
                    peptidseq_list, Mhcseq_list, random_letter=rand_let
                )
                # print("making peptide decoy 2")

                mhc_decoy = Mhcseq_str

            elif peptide_or_mhc_choice == "Mhcseq_list":
                # print("making mhc decoy")
                logging.info("making MHC decoy")
                mhc_decoy, _ = self._peptide_decoy_maker(
                    Mhcseq_list, peptidseq_list, random_letter=rand_let
                )
                peptide_decoy = peptidseq_str

            elif peptide_or_mhc_choice == "both" and rand_let == "random":
                # print("I am doing both")
                logging.info("making both peptide and MHC decoy")
                peptide_decoy, peptide_decoy_list = self._peptide_decoy_maker(
                    peptidseq_list, Mhcseq_list, random_letter=rand_let
                )
                mhc_decoy, _ = self._peptide_decoy_maker(
                    Mhcseq_list, peptide_decoy_list, random_letter=rand_let
                )

            if peptide_decoy is not None and mhc_decoy is not None:
                # print("do you ever come here?")
                temp_set.add((peptide_decoy, mhc_decoy, "False"))

        return list(temp_set)

    def random_decoys(self, peptid_str, mhc_str, decoy_count=10):
        """
        creates decoys as many as decoy_count for a give peptide_mhc pair
        """

        len_peptid = len(peptid_str)
        len_mhc = len(mhc_str)
        temp_set = set()

        while len(temp_set) < decoy_count:
            rand_peptide = PeptideMhcPair.random_letter_maker(len_peptid)
            rand_mhc = PeptideMhcPair.random_letter_maker(len_mhc)
            temp_set.add((rand_peptide, rand_mhc, "False")) #to ensure the decoys are not repeated

        return list(temp_set)

"""
#some test

amino_pairs = AminoAcids(
    amino_acid_letter,
)
pairs = amino_pairs.k_pair(k=40)
# print(pairs)
#
test = PeptideMhcPair(pairs=pairs)

# print("test.pairs", test.pairs)
# input()
peptidseq_str, peptidseq_list, Mhcseq_str, Mhcseq_list = test.peptide_maker()
# print("hhhh",peptidseq_str)
# a = test._peptide_decoy_maker(['DM', 'A', 'TLDI', 'Q', 'HI', 'Y', 'QK'],['DM', 'A', 'TLDI', 'Q', 'HI', 'Y', 'QK'])
# print(a)
# c = test.decoy_output('DMPTLDIRHIYQK', ['DM', 'P', 'TLDI', 'R', 'HI', 'Y', 'QK'],
#                       'DMPTLDIRHIYBB', ['DM', 'P', 'TLDI', 'R', 'HI', 'Y', 'BB'], decoy_count=10)
# # print(len(c))
# print("aaaa",c)
# d = test.random_decoys(10,20, decoy_count=10)
# print(d)
# print(len(d))
"""

class dataWrapper:
    def __init__(
        self,
        peptide_min,
        peptide_max,
        mhc_min,
        mhc_max,
        sample_count=1,
        decoy_count=9,
        pairs=6,
        binding_pos_count=3,
        random_decoys=False
        # decoy_let_num="random",
    ):

        self.sample_count = sample_count
        self.decoy_count = decoy_count
        # self.pairs = pairs  #self.random_pairs(pairs) if type(pairs) is int else pairs

        self.peptide_min = peptide_min
        self.peptide_max = peptide_max
        self.mhc_min = mhc_min
        self.mhc_max = mhc_max
        self.binding_pos_count = binding_pos_count
        self.random_decoys = random_decoys
        # self.decoy_let_num = decoy_let_num
        # print(pairs)
        # input()
        if type(pairs) == int:
            if pairs < 3:
                raise ValueError("the number of pairs should be 3 >=")
            self.pairs = self.random_pairs(length=pairs)
        else:
            if len(pairs) < 3:
                raise ValueError("the number of pairs given manually should be 3 >=")
            self.pairs = pairs

    def random_pairs(self, length=40, replacement=True):
        """
        This is a helper function to generate a list of random pairs to be chosen later for binding pairs
        args:
            length(int): number of random pairs for binding
        return:
            pair_list(list(tuples)): a list of tuples of binding amino acides
        """

        pairs_list = []
        # for _ in range(length):
        while len(pairs_list) < length:
            pair = random.sample(amino_acid_letter, k=2)

            if replacement:
                pairs_list.append(tuple(pair))

            elif not replacement:
                if pair not in pairs_list:  # this stop from replacement from the list of pairs
                    pairs_list.append(tuple(pair))

        return pairs_list

    def _min_max_checker(self, parts_len, min_len, max_len, mhc_or_peptide=None):
        """
        a funciton to check the pre_specified values for min and max length of peptide or mhc

        """

        peptide_lenghth = sum(parts_len) + self.binding_pos_count

        if peptide_lenghth < min_len or peptide_lenghth > max_len:
            raise ValueError(
                f"the sequence of {mhc_or_peptide} should have {min_len} minimum and {max_len} maximum length, {mhc_or_peptide} length adds up to {peptide_lenghth}"
            )

    def _random_indexes_letter(self, length=4):
        """
        args:
            length(int): number of spots between binding indeces. for 3 binding sites we have 4 spots
        """
        numbers = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        return random.choices(numbers, k=length)

    def _multiple_seq_maker(self, peptide_part_len=None, mhc_part_len=None):
        """
        args:
            peptide_part_len(list(int)): a list of indices of each binding index of peptide: [2, 2, 3, 4]
            mhc_part_len(list(int)): a list of indices of each binding index of mhc: [2, 2, 3, 4]
        return:
            peptide_str(str): a string of peptide, Mhc_str(str): a string of mhc, true value
            decoys(list(tuples)): a list of tuples of decoys


        """

        peptid_mhc_obj = PeptideMhcPair(self.pairs)

        peptidseq_str, peptidseq_list, Mhcseq_str, Mhcseq_list = (
            peptid_mhc_obj.peptide_maker(
                peptide_part_len=peptide_part_len, mhc_part_len=mhc_part_len
            )
        )
        if self.random_decoys:

            decoys = peptid_mhc_obj.random_decoys(
                peptidseq_str, Mhcseq_str, decoy_count=self.decoy_count
            )

        else:
            decoy_list = []  # for 1,2,3 or random substitute
            for i in [1, 2, 3, "random"]:
                decoys = peptid_mhc_obj.decoy_output(
                    peptidseq_str,
                    peptidseq_list,
                    Mhcseq_str,
                    Mhcseq_list,
                    decoy_count=self.decoy_count,
                    rand_let=i,
                )
                decoy_list.append(decoys)
        # print((peptidseq_str, Mhcseq_str,"True"), decoy_list)

        return (peptidseq_str, Mhcseq_str, "True"), decoy_list

    def print_to_files(
        self,
        true_counts,
        peptide_indeces,
        mhc_indeces,
        tail_peptide_min,
        tail_peptide_max,
        tail_mhc_min,
        tail_mhc_max,
        out_dir
    ):
        true_seqs = []
        decoys = []

        for _ in range(true_counts):
            peptide_indeces_ = [
                peptide_indeces[0],
                peptide_indeces[1],
                peptide_indeces[2],
                random.randint(tail_peptide_min, tail_peptide_max),
            ]
            mhc_indeces_ = [
                mhc_indeces[0],
                mhc_indeces[1],
                mhc_indeces[2],
                random.randint(tail_mhc_min, tail_mhc_max),
            ]

            self._min_max_checker(
                peptide_indeces[:3] + [tail_peptide_max],
                self.peptide_min,
                self.peptide_max,
                mhc_or_peptide="peptide",
            )
            self._min_max_checker(
                mhc_indeces[:3] + [tail_mhc_max],
                self.mhc_min,
                self.mhc_max,
                mhc_or_peptide="Mhc",
            )
            seqs, deqs = self._multiple_seq_maker(
                peptide_part_len=peptide_indeces_, mhc_part_len=mhc_indeces_
            )

            true_seqs.append(seqs)
            decoys.append(deqs)
        # each  corresponds with 1,2 or 3 places of change in decoys. Random is random places.
        peptide_decs_1 = []
        mhc_decs_1 = []
        peptide_decs_2 = []
        mhc_decs_2 = []
        peptide_decs_3 = []
        mhc_decs_3 = []
        peptide_decs_random = []
        mhc_decs_random = []
        labels_decs = []


        for each_true in decoys:
            for j in range(len(each_true)):

                for k in range(len(each_true[j])):

                    if j == 0:
                        peptide_decs_1.append(each_true[j][k][0])
                        mhc_decs_1.append(each_true[j][k][1])
                        labels_decs.append(each_true[j][k][2])
                    elif j == 1:
                        peptide_decs_2.append(each_true[j][k][0])
                        mhc_decs_2.append(each_true[j][k][1])
                    elif j == 2:
                        peptide_decs_3.append(each_true[j][k][0])
                        mhc_decs_3.append(each_true[j][k][1])
                    elif j == 3:
                        peptide_decs_random.append(each_true[j][k][0])
                        mhc_decs_random.append(each_true[j][k][1])


        df_decoys = pd.DataFrame(
            list(
                zip(
                    peptide_decs_1,
                    mhc_decs_1,
                    peptide_decs_2,
                    mhc_decs_2,
                    peptide_decs_3,
                    mhc_decs_3,
                    peptide_decs_random,
                    mhc_decs_random,
                    labels_decs,
                )
            ),
            columns=[
                "peptide_1",
                "mhc_1",
                "peptide_2",
                "mhc_2",
                "peptide_3",
                "mhc_3",
                "peptide_rand",
                "mhc_rand",
                "label",
            ],
        )
        df_not_decoys = pd.DataFrame(true_seqs, columns=["peptide", "mhc", "label"])

        print(len(df_not_decoys["peptide"]))
        print(len(df_not_decoys["peptide"].unique()))
        print(len(df_not_decoys["mhc"]))
        print(len(df_not_decoys["mhc"].unique()))

        train_decoy, test_decoy = train_test_split(
            df_decoys, test_size=0.2, shuffle=False
        )  # from scikit-learn
        train_true, test_true = train_test_split(
            df_not_decoys, test_size=0.2, shuffle=False
        )

        train_decoy.to_csv(f"{out_dir}/decoys_dataset_train.csv")
        test_decoy.to_csv(f"{out_dir}/decoys_dataset_test.csv")

        train_true.to_csv(f"{out_dir}/true_dataset_train.csv")
        test_true.to_csv(f"{out_dir}/true_dataset_test.csv")


def main():
    sample_num = args.sample_num
    decoys_ratio = args.decoys
    min_peptide_len = args.peptide_min
    max_peptide_len = args.peptide_max
    min_mhc_len = args.mhc_min
    max_mhc_len = args.mhc_max
    random_decoy_flag = args.random_decoy
    # number_of_let_in_decoy = args.num_letter_decoy
    out_dir = args.out_dir

    tail_peptide_min = 1
    tail_peptide_max = 2
    tail_mhc_min = 5
    tail_mhc_max = 20
    peptide_indeces = [2, 2, 3]  # for manual indexing
    mhc_indeces = [4, 2, 5]  # for manual indexing

    # mhc_indeces = [7,5,3, 30]

    # pairs_cosutm = [('K', 'V'), ('K', 'L'), ('C', 'M'), ('I', 'T'), ('P', 'W'), ('A', 'M'), ('E', 'D'), ('L', 'Q'), ('S', 'N'), ('A', 'L'),
    #          ('H', 'T'), ('L', 'F'), ('F', 'M'), ('G', 'T'), ('F', 'G'), ('R', 'N'), ('D', 'M'), ('E', 'S'), ('F', 'V'), ('I', 'F'),
    #          ('Y', 'V'), ('G', 'Q'), ('W', 'W'), ('I', 'H'), ('E', 'I'), ('F', 'P'), ('Q', 'H'), ('M', 'C'), ('C', 'A'), ('E', 'P'),
    #          ('K', 'D'), ('A', 'A'), ('D', 'E'), ('V', 'L'), ('R', 'V'), ('G', 'M'), ('E', 'Y'), ('G', 'I'), ('G', 'R'), ('L', 'E')]

    # pairs = [('K', 'V'), ('K', 'L'),  ('I', 'T'), ('P', 'W'), ('A', 'M'), ('E', 'D'), ('L', 'Q'), ('S', 'N'), ('A', 'L'),
    #          ('H', 'T'), ('G', 'T'), ('F', 'G'), ('R', 'N'), ('D', 'M'), ('E', 'S'),
    #          ('Y', 'V'), ('G', 'Q'), ('I', 'H'), ('E', 'I'), ('F', 'P'), ('Q', 'H'), ('C', 'A'), ('E', 'P'),
    #          ('K', 'D'), ('D', 'E'), ('R', 'V'), ('G', 'M'), ('E', 'Y'), ('G', 'I'), ('G', 'R'), ('L', 'E')] #to use for the hard case

    amino_pairs = AminoAcids(
        amino_acid_letter,
    )
    pairs = amino_pairs.k_pair(k=40)
    # pairs = pairs[30:]
    print(pairs)
    print(len(pairs))
    # input()
    # pairs = 40
    test2 = dataWrapper(
        min_peptide_len,
        max_peptide_len,
        min_mhc_len,
        max_mhc_len,
        sample_count=sample_num,
        decoy_count=decoys_ratio,
        pairs=pairs,
        random_decoys=random_decoy_flag,
    )

    # print(test2.pairs)

    # a, b =test2.multiple_seq_maker(peptide_part_len= [2,2,3, random.randint(1, 10)],mhc_part_len = [7,2,3, random.randint(9, 15)])
    test2.print_to_files(
        sample_num,
        peptide_indeces,
        mhc_indeces,
        tail_peptide_min,
        tail_peptide_max,
        tail_mhc_min,
        tail_mhc_max,
        out_dir
    )


if __name__ == "__main__":

    main()
