"""
Copyright 2018 Iris van de Pol and Shane Steinert-Threlkeld.

This file is part of quantifier-LZ-complexity.

quantifier-LZ-complexity is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

quantifier-LZ-complexity is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with learn-complex.  If not, see <https://www.gnu.org/licenses/>.
"""


import time
from bitarray import bitarray
import numpy as np


def merge_bitarrays_randomly(bitarray1,bitarray2):
    ''' Takes as input two bitarrays,
        it merges them in a random way
        in such a way that the order of elements within each of the two input bitarrays
        remains the same.

    Returns: 
        a bitarray
    '''
    guide_bitarray = bitarray(len(bitarray1))
    guide_bitarray.setall(True)
    guide_bitarray2 = bitarray(len(bitarray2))
    guide_bitarray2.setall(False)
    guide_bitarray.extend(guide_bitarray2)
    np.random.shuffle(guide_bitarray)

    merged_bitarray = bitarray(len(guide_bitarray))

    count_1 = 0
    count_2 = 0
    for idx in range(len(merged_bitarray)):
        if guide_bitarray[idx] == 1:
            merged_bitarray[idx] = bitarray1[count_1]
            count_1+=1
        elif guide_bitarray[idx] == 0:
            merged_bitarray[idx] = bitarray2[count_2]
            count_2+=1
    return merged_bitarray;

def generate_filename(quantifier, max_model_size, num_chars, order_number):
    """ Generates a filename, which is used for storing and retrieving auxilary files, 
        used for making quantifier representations based on random model orders.

    Args:
        quantifier: a bitarray
        max_model_size: an integer
        num_chars:
        order_number: an integer

    Returns: a string
    """

    filename = str(num_chars)+"char/"
    filename += "aux_"
    filename += str(quantifier)+"_"
    filename += "size"+str(max_model_size)+"_"
    filename += str(num_chars)+"char_"   
    filename += "order_num"+str(order_number)
    filename += ".ba"
    return filename

def store_bitarray(seq, filename):
    """ Stores a bitarray as a binary file in filename, 
        stores the length of the bitarray in filename.len
        (storing a bitarray adds trailing 0's to make a multiple of 8, enabling storage using bytes).

    Args:
        seq: a bitarray
        filename: a string
    """

    with open(filename,'wb') as output_file:
        seq.tofile(output_file)
    with open(filename+'.len','w') as output_file:
        output_file.write(str(len(seq)))

def load_bitarray(filename):
    """ Loads a bitarray stored as a binary file in filename, 
        using the length of the bitarray stored in filename.len.

    Args:
        filename: a string

    Returns:
        a bitarray
    """

    seq = bitarray()
    with open(filename,'rb') as input_file:
        seq.fromfile(input_file)
    with open(filename+'.len','r') as input_file:
        seq_len = int(input_file.read())
    seq = seq[0:seq_len]
    return seq







        







