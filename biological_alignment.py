"""
This program implements the Needleman-Wunsch algorithm to perform biological alignment. It uses dynamic programming as
the base of workflow of the algorithm to create a series of matrices that score the mach, mismatch and gap between two
sequences.

Date: March 19, 2022

Author: Ricardo Y. Rodriguez Gonzalez
Student Number: 802-18-2754

Course: CIIC 4025

Professor: Wilfredo E. Lugo-Beauchamp
"""


# ________________________________________________Biological Alignment________________________________________________ #

#                                              Needleman-Wunsch Algorithm


# ----------------------------------------------Importing Modules----------------------------------------------------- #

import numpy as np
import sys
import csv

# -------------------------------------------------input.csv Reader--------------------------------------------------- #


'''
 This method has an initial empty string that is used to store the data that is being extracted 
 form the file the is inputted in csv format. This file is opened and read using the csv.reader() method. 
 
 * @name - data_extract()
 * @return - It returns the variable sequence_data that stores the extracted sequence as a string.
 
'''

def data_extract():
    sequence_data = ""
    try:
        file = open(sys.argv[1], 'r')
        sequence_data = csv.reader(file)
    except:
        print("Invalid file path or file.")
        exit()
    return sequence_data


# -------------------------------------Matrix with Needleman-Wunsch Algorithm----------------------------------------- #

'''
 This method is an implementation of the Needleman-Wunsch algorithm for sequence alignment. It is used to compare and 
 analyze the alignment of two biological sequences. 
 
 For example:
 
        Sequence 1: "ATCGT" ------> Aligned Sequence 1: "ATCGT"
        Sequence 2: "ACGT"  ------> Aligned Sequence 2: "A-CGT"
        
        Socre = 2

 * @name - needleman_wunsch()
 * @param - sequence_1, sequence_2
 * @type - sequence_1: str, sequence_2: str
 * @return - It returns a list of the alignment texts and the score values as integers. The expected return values 
            are a column of alignment text and another column of scores.

'''


def needleman_wunsch(sequence_1: str, sequence_2: str):
    # -------------------------------------------------Matrix Creation------------------------------------------------ #

    main_matrix = np.zeros((len(sequence_1) + 1, len(sequence_2) + 1))

    matching_matrix = np.zeros((len(sequence_1), len(sequence_2)))

    # -----------------------------Score assignment for match, mismatch, and gap penalty------------------------------ #

    reward_match = 1
    penalty_mismatch = -1
    gap_penalty = -2

    # ---------------------------------------Filling Match/Mismatch Matrix-------------------------------------------- #

    for i in range(len(sequence_1)):
        for j in range(len(sequence_2)):
            if sequence_1[i] == sequence_2[j]:
                matching_matrix[i][j] = reward_match
            else:
                matching_matrix[i][j] = penalty_mismatch

    # ---------------------------------------------------------------------------------------------------------------- #

    # Step 1: Initialize gap-penalty values in matrix

    for i in range(len(sequence_1) + 1):
        main_matrix[i][0] = i * gap_penalty
    for j in range(len(sequence_2) + 1):
        main_matrix[0][j] = j * gap_penalty

    # Step 2: Fill Matrix:

    for i in range(1, len(sequence_1) + 1):
        for j in range(1, len(sequence_2) + 1):
            main_matrix[i][j] = max(main_matrix[i - 1][j - 1] + matching_matrix[i - 1][j - 1],
                                    main_matrix[i - 1][j] + gap_penalty,
                                    main_matrix[i][j - 1] + gap_penalty)
    score = main_matrix[-1][-1]

    # Step 3: Backtracking:

    aligned_1 = ""
    aligned_2 = ""

    seq_i, seq_j = len(sequence_1), len(sequence_2)

    while seq_i > 0 or seq_j > 0:
        if seq_i > 0 and seq_j > 0 and main_matrix[seq_i][seq_j] == main_matrix[seq_i - 1][seq_j - 1] + \
                matching_matrix[seq_i - 1][seq_j - 1]:
            aligned_1 = sequence_1[seq_i - 1] + aligned_1
            aligned_2 = sequence_2[seq_j - 1] + aligned_2

            seq_i -= 1
            seq_j -= 1

        elif seq_i > 0 and main_matrix[seq_i][seq_j] == main_matrix[seq_i - 1][seq_j] + gap_penalty:
            aligned_1 = sequence_1[seq_i - 1] + aligned_1
            aligned_2 = "-" + aligned_2

            seq_i -= 1

        else:
            aligned_1 = "-" + aligned_1
            aligned_2 = sequence_2[seq_j - 1] + aligned_2

            seq_j -= 1

    return ["".join(i for i in [aligned_1, " ", aligned_2]), int(score)]


# ---------------------------------------------------Runnable--------------------------------------------------------- #
'''
This section is where the data is stored into a variable and used to write the first to columns 
of the output file that is being created and written named: results.csv. Text to identify each column is 
also added to the file. A for loop writes each row of the newly established columns [alignment text] and 
[alignment score]. The needleman_wunsch() method is called to fill each row with the proper values. The alignment text 
section receives the alignment as strings and the alignment score receives the integer values of each alignment. 
'''

data = data_extract()

with open('./results.csv', 'w', newline='') as output_csv:
    out = csv.writer(output_csv, delimiter=",")
    m = next(data) + ['alignment text', 'alignment score']
    out.writerow(m)

    for row in data:
        row += needleman_wunsch(row[0], row[1])
        out.writerow(row)
    output_csv.close()
