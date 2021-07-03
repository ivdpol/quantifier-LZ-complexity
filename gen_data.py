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

import quantifiers
import util
from itertools import product
from itertools import permutations
from multiprocessing import Pool, set_start_method
import argparse
from bitarray import bitarray
from lempel_ziv_complexity import lempel_ziv_complexity
import pandas as pd
import gzip
import bz2
import time
from math import log2
import numpy as np
import os
import shutil

def quantifiers_reps_lex(quantifiers, max_model_size, num_perms, num_chars,min_model_size=1):
    """ Produces num_perms many binary representations of a given quantifier,
        based on a lexicographical order on the models.
        Such a representation is produced by checking for all quantifier models  
        of at most max_model_size (in lexicographic order), 
        if they model that quantifier and then placing the answers in a bitarray.
        A quantifier model is represented by a sequence of elements from (0,1,2,3), 
        representing (AB,AnotB,BnotA,neither) for num_chars = 4, 
        or elements from (0,1,3) representing (AB,AnotB,BnotA) for num_chars = 3,
        or elements from (0,1) representing (AB,AnotB) for num_chars = 2.
        The i'th element in the sequence represents the i'th object in the model.
        Each representation is taken over a different ordering of quantifier models. 
        For each ordering a different permutation over (0,1,2,3) is used as the base order
        for a lexicographical ordering over all quantifier models.

    Args:
        quantifiers: a list of quantifier objects
        max_model_size: an integer, determining the size (i.e., number of elements) of the quantifier models
        num_perms: an integer, determinig the number of permutations over (0,1,2,3);
                    i.e., the number of different base orders for the lexicographical ordering over the models
        num_chars: an integer, determining the number of characters used in the quantifier representation,
                    (the number of subareas—--AB, AnotB, BnotA, and neither---of the quantifier model that are considered)
                    needs to be either 2, 3, or 4
        min_model_size: an integer, 1 for the 3 and 4 char case, 2 for the 2 char case

    Returns:
        a list of list (one per quantifier) with num_perms many bitarrays, 
        in each bitarray the value at index i represents whether a given quantifier 
        holds in the i'th model in a lexicographical order of quantifier models,
        in each list, the different bitarray corresponds 
        with a different lexicographical order over models (based on a different base order over the character)
    """

    if not num_chars in [2,3,4]:
        raise ValueError('Supply num_chars=2, num_chars=3, or num_chars=4!')

    if num_chars == 2:
        # check conservativity
        for quantifier in quantifiers:
            if quantifier.cons == False: 
                raise ValueError('Supply a conservative quantifier!')
        if num_perms != 1:
            raise ValueError('Supply num_perms = 1')
        if max_model_size < 2:
            raise ValueError('For char = 2, supply max_model_size > 1')

    if num_chars == 3:
        if num_perms > 3 or num_perms < 1:
            raise ValueError('Supply a num_perms <= 3 and > 1')
        if max_model_size < 1:
            raise ValueError('Supply max_model_size > 0')

    if num_chars == 4:
        if num_perms > 12 or num_perms < 1:
            raise ValueError('Supply a num_perms <= 12 and > 1')
        if max_model_size < 1:
            raise ValueError('Supply max_model_size > 0')

    # generate a list of lists of dummy bitarrays of the correct length in which the output values will be stored 
    quants_reps = [[bitarray(sum([num_chars**x for x in range(min_model_size, max_model_size+1)])) for _ in range(num_perms)] for _ in range(len(quantifiers))]

    # generate the unique (meaning excluding reversed ones) base orders (perms) for the lexicographical orderings 
    # use the property that the first symbol < the last symbol (or vice versa) for all orders that are not the reverse of each other
    perms = [order for order in list(permutations(list(range(num_chars)))) if order[0] < order[-1]][:num_perms]

    # generate all quantifier models of size at most max_model_size 
    # which are all of the sequences of length at most max_model_size over alphabet {0,...,num_chars}
    # in a lexicographical order, given some permutation over (0,...,num_chars ) as a base
    # for each quantifier model and for each quantifier, check whether it models the given quantifier and put the answer in a bitarray 
    
    for order_idx in range(len(perms)):
        order = perms[order_idx] 
        start_of_block = 0
        for model_size in range(min_model_size,max_model_size+1):
            models = product(order, repeat=model_size)
            for model_idx in range(num_chars**model_size):
                # check whether the quantifier holds in this model 
                # and iff so put a 1 at that model_idx in the quant_rep
                model = next(models)
                for quant_idx in range(len(quantifiers)):
                    quantifier = quantifiers[quant_idx]
                    quants_reps[quant_idx][order_idx][start_of_block + model_idx] = quantifier(model)
            start_of_block += num_chars**model_size

    return quants_reps


def get_uniformity(bin_seq):
    """Computes the ratio of 0's in a binary sequence (or the ratio of 1's, if that is larger)

    Args: a binary sequence

    Returns: an integer
    """
    avg = sum(bin_seq) / len(bin_seq)
    return max(avg, 1-avg)


def get_row_of_table(quant_order_rep):
    """ A top level function, used in pool.map for parallelizaiton
        i. can only take one argument
        ii. need to be at top-level, so picklable
        Computes LZ, bzip2, and gzip complexities of a quantifier representation 
        and returns it along with other metrics about the quantifier,
        so that a row can be added to a pandas data frame.

    Args: a 7-tuple, containing a quantifier representation and information about it

    Returns: a dictionary, to be used for appending to a data frame
    """

    max_model_size, quantifier, order_num, order_type, rep, uniformity, num_chars = quant_order_rep
    
    rev_rep = rep[::-1]
    rep_len = len(rep)
    LZ_complexity = lempel_ziv_complexity(rep)
    LZ_rev_complexity = lempel_ziv_complexity(rev_rep)
    gzip_complexity = len(gzip.compress(rep.tobytes()))
    gzip_rev_complexity = len(gzip.compress(rev_rep.tobytes()))
    bzip_complexity = len(bz2.compress(rep.tobytes()))
    bzip_rev_complexity = len(bz2.compress(rev_rep.tobytes()))

    return {'max_model_size': max_model_size,
            'quantifier': quantifier.name,
            'order_num': order_num,
            'order_type': order_type,
            'LZ': LZ_complexity,
            'LZ_rev': LZ_rev_complexity,
            'C_LZ': log2(rep_len)*(LZ_complexity + LZ_rev_complexity)/2,
            'gzip': gzip_complexity,
            'gzip_rev': gzip_rev_complexity,
            'C_gzip': (gzip_complexity + gzip_rev_complexity)/2,
            'bzip': bzip_complexity,
            'bzip_rev': bzip_rev_complexity,
            'C_bzip': (bzip_complexity + bzip_rev_complexity)/2,
            'uniformity': uniformity,
            'monotonicity': quantifier.mono,
            'quantity': quantifier.quan,
            'conservativity': quantifier.cons,
            'num_chars': num_chars
            }


def LZ_dist_quantifier(quantifiers, max_model_size, num_model_orders, num_chars=4, 
                            num_procs=1, out_file=None, order_type = 'lex'):
    """ Generates num_model_orders many different lexicographical- or random-order quantifier representations for each quantifier in quantifiers,
        and computes the LZ, bzip2, and gzip complexities for all of those quantifier representations

    Args:
        quantifiers: a list of quantifier objects
        max_model_size: an integer, determining the size (i.e., number of elements) of the quantifier models
        num_model_orders: an integer, determing the number of different model order permutations
        num_procs: an integer, determining the number of prcessors used for parallel computing
        num_chars: an integer, determining the number of characters used in the quantifier representation,
                    (the number of subareas—--AB, AnotB, BnotA, and neither---of the quantifier model that are considered)
                    needs to be either 2, 3, or 4
        out_file: a string, determining the name of the csv file in which the results are stored
        order_type: either 'lex' (for lexicographic model orders) or 'rand' (for random model orders)

    Returns:
        a pandas data frame which contains the complexity values and other metrics of the quantifiers
    """

    time1 = time.time()
    data = pd.DataFrame()
    pool = Pool(num_procs)

    if (order_type == 'lex'):

        # get the quantifier representations for the different lexicographical orders for each quantifiers
        quantifiers_reps = quantifiers_reps_lex(quantifiers,max_model_size,num_model_orders,num_chars)

        # get uniformities per quantifier
        uniformities = [get_uniformity(reps[0]) for reps in quantifiers_reps]

        # get quantifier representations for each quantifier per model_order (e.g., quants_reps_per_order = [Q1_rep1, Q2_rep1, ... ])
        for order_num in range(num_model_orders):
            quants_reps_per_order = [reps[order_num] for reps in quantifiers_reps]

            # store all the LZ complexities and other metrics in a pandas data frame
            data = data.append(pool.map(
                    get_row_of_table,
                    zip([max_model_size]*len(quantifiers),quantifiers, [order_num]*len(quantifiers),[order_type]*len(quantifiers), quants_reps_per_order,
                        uniformities, [num_chars]*len(quantifiers))),
                ignore_index=True)

    if order_type == 'rand':
        
        # make folder for storing the quantifier representations so that they can be retrieved and merged 
        # when computing the representation for larger sizes (with another function call)
        if not os.path.exists("aux"):
           os.mkdir("aux")
        if not os.path.exists("aux/"+str(num_chars)+"char"):
           os.mkdir("aux/"+str(num_chars)+"char")

        # get the quantifier representations for one lexicographical orders for each quantifiers
        # only include the maximum model size in the representations
        min_model_size = max_model_size
        if num_chars == 2  and max_model_size == 2:
            min_model_size = 1
        current_reps = [q_rep for sublist in quantifiers_reps_lex(quantifiers, max_model_size, 1, num_chars, min_model_size) for q_rep in sublist]

        min_model_size = 1  
        # for the 2 char case model size 1 gives a bug, this is related to a bug in storing a bitarray of size 2 as a binary file
        if num_chars == 2:
            min_model_size = 2

        for order_num in range(num_model_orders):

            # make a (new) random-order representation for each quantifier by shuffling
            for quant in current_reps:
                np.random.shuffle(quant) 

            # when this function call is for the smallest model size, store the representations so that they can be 
            # retrieved and merged with a larger model size in a future function call
            if max_model_size == min_model_size:
                quants_reps = current_reps
                for quant_idx in range(len(quantifiers)):
                    util.store_bitarray(quants_reps[quant_idx], 
                        "aux/"+util.generate_filename(quantifiers[quant_idx], max_model_size, num_chars, order_num))

            elif max_model_size > min_model_size:
                # merge the quant representations of the current size with the representation of the previous size
                # merge them in a random way, in such a way that the order of the elements within 
                # the current size en previous size representation remain the same
                quants_reps = [None for _ in range(len(quantifiers))]
                for quant_idx in range(len(quantifiers)):
                    aux_filename = "aux/"+util.generate_filename(quantifiers[quant_idx], max_model_size-1, num_chars, order_num)
                    if os.path.exists(aux_filename):
                       prev_rep = util.load_bitarray(aux_filename)
                    else:
                       raise Exception("Needed temporary file not found: " + aux_filename);
                    quants_reps[quant_idx] = util.merge_bitarrays_randomly(prev_rep,current_reps[quant_idx])
                    # store the merged representation so that they can be retrieved and merged
                    # with a larger model size in a future function call
                    util.store_bitarray(quants_reps[quant_idx], 
                        "aux/"+util.generate_filename(quantifiers[quant_idx], max_model_size, num_chars, order_num))

            # get uniformities per quantifier
            uniformities = [get_uniformity(rep) for rep in quants_reps]

            # store all the LZ complexities and other metrics in a pandas data frame       
            data = data.append(pool.map(
                    get_row_of_table,
                    zip([max_model_size]*len(quantifiers), quantifiers, [order_num]*len(quantifiers),[order_type]*len(quantifiers), quants_reps,
                        uniformities, [num_chars]*len(quantifiers))),
                ignore_index=True) 

    # store the data in a csv file
    if out_file:
        if not os.path.exists("data"):
            os.mkdir("data")
        if not os.path.exists("data/"+str(num_chars)+"char/"):
            os.mkdir("data/"+str(num_chars)+"char/")
        if not os.path.exists("data/"+str(num_chars)+"char/"+order_type+"/"):
            os.mkdir("data/"+str(num_chars)+"char/"+order_type+"/")
        data.to_csv("data/"+str(num_chars)+"char/"+order_type+"/"+out_file)

    time2 = time.time()
    print('max_model_size =', max_model_size, ',', time2-time1, '\n')
    return data


def LZ_dist_quantifier_cumulative(quantifiers, max_model_size, num_model_orders, num_chars=4,
                            num_procs=1, out_file='Complexity_dist_results.csv', order_type=None):
    """ Incrementally builds a pandas data frame with the distribution of quantifier complexities over num_model_orders
        many different lexicographical- or random-order quantifier representations for each quantifier in quantifiers,
        and for each value for max_model_size, from 1 (or from 2, in the 2 char case) till max_model_size.
        A quantifier representation for some maximum model size, includes all previous sizes.

    Args:
        quantifiers: a list of quantifier objects
        max_model_size: an integer, determining the size (i.e., number of elements) of the quantifier models
        num_model_orders: an integer, determing the number of different model order permutations
        num_chars: an integer, determining the number of characters used in the quantifier representation,
                    (the number of subareas—--AB, AnotB, BnotA, and neither---of the quantifier model that are considered)
                    needs to be either 2, 3, or 4
        num_procs: an integer, determining the number of processors used for parallel computing
        out_file: a string, determining the name of the csv file in which the results are stored
        order_type: either 'lex' (for lexicographic model orders) or 'rand' (for random model orders)

    Returns:
        a pandas data frame which contains the complexity values and other metrics of the quantifiers
    """

    # delete the auxilirary files with supporting information for building the random-order representations, 
    # that were made by a previous function call, if any
    if order_type == 'rand':
        if os.path.exists("aux/"+str(num_chars)+"char/"):
            shutil.rmtree("aux/"+str(num_chars)+"char/")

    big_data_table = pd.DataFrame()

    min_model_size = 1  
    # for the 2 char case model size 1 gives a bug, this is related to a bug in storing a bitarray of size 2 as a binary file
    if num_chars == 2:
        min_model_size = 2

    # incrementally build the data frame
    for model_size in range(min_model_size,max_model_size+1):
        big_data_table = big_data_table.append(LZ_dist_quantifier(quantifiers, model_size, num_model_orders, 
                           num_chars, num_procs, None, order_type), ignore_index = True)

    # store the data in a csv file
    if not os.path.exists("data"):
        os.mkdir("data")
    if not os.path.exists("data/"+str(num_chars)+"char/"):
        os.mkdir("data/"+str(num_chars)+"char/")
    if not os.path.exists("data/"+str(num_chars)+"char/"+order_type+"/"):
        os.mkdir("data/"+str(num_chars)+"char/"+order_type+"/")
    big_data_table.to_csv("data/"+str(num_chars)+"char/"+order_type+"/"+out_file)
    return big_data_table


if __name__ == '__main__':

    # NOTE: this is to avoid running out of memory in large Pools
    set_start_method('spawn')

    start = time.time()

    all_cons_quantifiers = quantifiers.get_all_cons_quantifiers()
    all_quantifiers = quantifiers.get_all_quantifiers()

    parser = argparse.ArgumentParser()
    parser.add_argument('--max_model_size', type=int, default=6)
    parser.add_argument('--num_model_orders', type=int, default=12)
    parser.add_argument('--out_file', type=str, default='Complexity_results.csv')
    parser.add_argument('--num_procs', type=int, default=1)
    parser.add_argument('--num_chars', type=int, default=4)
    parser.add_argument('--quantifiers', type=str, default='all_quantifiers')
    parser.add_argument('--jobtype', type=str, default='cumulative')
    parser.add_argument('--order_type', type=str, default='quantity_lex')
    args = parser.parse_args()

    print('*********** JOB ***********')
    print('LZ_dist_quantifier_'+args.jobtype)
    print('max_model_size =', args.max_model_size)
    print('num_chars =', args.num_chars)
    print('order_type =', args.order_type)
    print('num_model_orders =', args.num_model_orders)
    print('num_procs =', args.num_procs)
    print('quantifiers =', args.quantifiers)
    print('***************************\n\n')
    
    chosen_quantifiers = None;
    if args.quantifiers == 'all_cons_quantifiers':
        chosen_quantifiers = all_cons_quantifiers;
    elif args.quantifiers == 'all_quantifiers':
        chosen_quantifiers = all_quantifiers;


    if args.jobtype == 'cumulative':
        print(LZ_dist_quantifier_cumulative(chosen_quantifiers,
                       args.max_model_size,
                       args.num_model_orders,
                       args.num_chars,
                       args.num_procs,
                       args.out_file,
                        args.order_type))

    elif args.jobtype == 'one_size':
        print(LZ_dist_quantifier(chosen_quantifiers,
                       args.max_model_size,
                       args.num_model_orders,
                       args.num_chars,
                       args.num_procs,
                       args.out_file,
                       args.order_type))

    end = time.time()

    print('running_time_this_program =', end-start)






















