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

import gc


class Quantifier:

    def __init__(self, name, mono=None, quan=None, cons=None, fn=None):

        if mono is None:
            raise ValueError("supply a value for monotonicity!")
        if quan is None:
            raise ValueError("supply a value for quantity!")
        if cons is None:
            raise ValueError("supply a value for conservativity!")
        if fn is None:
            raise ValueError("supply a function for verifying a quantifier!")

        self.name = name
        self.mono = mono
        self.quan = quan
        self.cons = cons
        self.verify = fn

    def __call__(self, seq):
        return self.verify(seq)

    def __repr__(self):
       return(self.name);


# a quantifier model is a model mathcal{M} = <M,A,B>, with A,B subseteq M,
# mathcal{M} has four different area's, A cap B, A - B, B - A, M - (A cup B)
# we represent an object in the model by an integer that signals in which area that object is 
# a model is represented by a sequence of integers
# the objects in the model are ordered, the index in the sequence represents the order 
# for readability, we name the integers:
AB = 0 #A cap B
AnotB = 1 #A - B
BnotA = 2 #B - A 
neither = 3 #M - (A cup B)


def models_at_least_4(model):
    """ Verifies whether at least 4 As are Bs in a sequence.

    Args:
        model:  an n-tuple, a sequence of n elements from (0,1,2,3), representing (AB,AnotB,BnotA,neither),
                the i'th element in the sequence represents the i'th object in the model

    Returns:
        True iff at least 4 elements in model are AB
    """
    count_AB = 0
    for symbol in model:
        if symbol == AB:
            count_AB += 1
            if count_AB == 4:
                return True
    return False


# make quantifier object 
at_least_4 = Quantifier("at_least_4",
                        mono="up", quan=True, cons=True,
                        fn=models_at_least_4)


# def models_at_least_6(model):
#     """ Verifies whether at least 6 As are Bs in a sequence.

#     Args:
#         model:  a sequence of elements from (0,1,2,3), representing (AB,AnotB,BnotA,neither),
#                 the i'th element in the sequence represents the i'th object in the model

#     Returns:
#         True iff at least 4 elements in model are AB
#     """
#     count_AB = 0
#     for symbol in model:
#         if symbol == AB:
#             count_AB += 1
#             if count_AB == 6:
#                 return True
#     return False


# # make quantifier object 
# at_least_6 = Quantifier("at_least_6",
#                         mono="up", quan=True, cons=True,
#                         fn=models_at_least_6)


def models_at_least_6_or_at_most_2(model):
    """ Verifies whether at least 6 or at most 2 As are Bs in a sequence.

    Args:
        model:  a sequence of elements from (0,1,2,3), representing (AB,AnotB,BnotA,neither),
                the i'th element in the sequence represents the i'th object in the model

    Returns:
        True iff at least 6 or at most 2 elements in model are AB
    """
    count_AB = 0
    for symbol in model:
        if symbol == AB:
            count_AB += 1
            if count_AB == 6:
                return True
    return count_AB <= 2 


# make quantifier object 
at_least_6_or_at_most_2 = Quantifier("at_least_6_or_at_most_2",
                                     mono="non", quan=True, cons=True,
                                     fn=models_at_least_6_or_at_most_2)


def models_even(model):
    """ Verifies whether an even number of As are Bs in a sequence.

    Args:
        model:  a sequence of elements from (0,1,2,3), representing (AB,AnotB,BnotA,neither),
                the i'th element in the sequence represents the i'th object in the model

    Returns:
        True iff an even number of elements in model are AB
    """
    count_AB = 0
    for symbol in model:
        if symbol == AB:
            count_AB += 1
    return count_AB % 2 == 0


# make quantifier object 
even = Quantifier("even",
                  mono="non", quan=True, cons=True,
                  fn=models_even)


def models_at_most_3(model):
    """ Verifies whether at most 3 As are Bs in a sequence.

    Args:
        model:  a sequence of elements from (0,1,2,3), representing (AB,AnotB,BnotA,neither),
                the i'th element in the sequence represents the i'th object in the model

    Returns:
        True iff at most 3 elements in model are AB
    """
    count_AB = 0
    for symbol in model:
        if symbol == AB:
            count_AB += 1
            if count_AB == 4:
                return False
    return True


# make quantifier object 
at_most_3 = Quantifier("at_most_3",
                       mono="down", quan=True, cons=True,
                       fn=models_at_most_3)


def models_at_least_3(model):
    """ Verifies whether at least 3 As are Bs in a sequence.

    Args:
        model:  a sequence of elements from (0,1,2,3), representing (AB,AnotB,BnotA,neither),
                the i'th element in the sequence represents the i'th object in the model

    Returns:
        True iff at least 3 elements in model are AB
    """
    count_AB = 0
    for symbol in model:
        if symbol == AB:
            count_AB += 1
            if count_AB == 3:
                return True
    return False


# make quantifier object 
at_least_3 = Quantifier("at_least_3",
                        mono="up", quan=True, cons=True,
                        fn=models_at_least_3)


def models_first_3(model):
    """ Verifies whether the first 3 As in a sequence are Bs.

    Args:
        model:  a sequence of elements from (0,1,2,3), representing (AB,AnotB,BnotA,neither),
                the i'th element in the sequence represents the i'th object in the model

    Returns:
        True iff the first 3 elements of model that are either AB or AnotB are in fact AB.
        False if either model has length less than n, or if on of the first 3 elements of model that are AB or AnotB are in fact AnotB,
        or if fewer than 3 elements in model are As
    """
    if len(model) < 3:
        return False
    count_AB = 0
    for symbol in model:
        if symbol == AB:
            count_AB += 1
            if count_AB == 3:
                return True
        if symbol == AnotB:
            return False
    # there are less than 3 As
    return False


# make quantifier object 
first_3 = Quantifier("first_3",
                     mono="up", quan=False, cons=True,
                     fn=models_first_3)


def models_last_3(model):
    """ Verifies whether the last 3 As in a sequence are Bs.

    Args:
        model:  a sequence of elements from (0,1,2,3), representing (AB,AnotB,BnotA,neither),
                the i'th element in the sequence represents the i'th object in the model

    Returns:
        True iff the last 3 elements of model that are either AB or AnotB are in fact AB.
        False if either model has length less than n, or if the last 3 elements of model that are either AB or AnotB are in fact AB,
        or if fewer than 3 elements in model are As
    """
    return models_first_3(tuple(reversed(model)))


# make quantifier object 
last_3 = Quantifier("last_3",
                    mono="up", quan=False, cons=True,
                    fn=models_last_3)


def models_not_all(model):
    """ Verifies whether not all As are Bs in a sequence.

    Args:
        model:  a sequence of elements from (0,1,2,3), representing (AB,AnotB,BnotA,neither),
                the i'th element in the sequence represents the i'th object in the model

    Returns:
        True iff not all elements in model are AB
    """
    for symbol in model:
        if symbol == AnotB:
            return True
    return False


# make quantifier object 
not_all = Quantifier("not_all",
                     mono="down", quan=True, cons=True,
                     fn=models_not_all)


def models_not_only(model):
    """ Verifies whether not only As are Bs in a sequence, i.e., whether there is a BnotA in a sequence.

    Args:
        model:  a sequence of elements from (0,1,2,3), representing (AB,AnotB,BnotA,neither),
                the i'th element in the sequence represents the i'th object in the model

    Returns:
        True iff there is a BnotA in model
    """
    for symbol in model:
        if symbol == BnotA:
            return True
    return False


# make quantifier object 
not_only = Quantifier("not_only",
                      mono="up", quan=True, cons=False,
                      fn=models_not_only)


def models_most(model):
    """ Verifies whether more than half As are Bs in a sequence.

    Args:
        model:  a sequence of elements from (0,1,2,3), representing (AB,AnotB,BnotA,neither),
                the i'th element in the sequence represents the i'th object in the model

    Returns:
        True iff more than half of the elements in model are AB
    """
    count_AB = 0
    count_AnotB = 0
    for symbol in model:
        if symbol == AB:
            count_AB += 1
        if symbol == AnotB:
            count_AnotB += 1
    return count_AB > count_AnotB


# make quantifier object 
most = Quantifier("most",
                  mono="up", quan=True, cons=True,
                  fn=models_most)


def models_M(model):
    """ Verifies whether |A| > |B|

    Args:
        model:  a sequence of elements from (0,1,2,3), representing (AB,AnotB,BnotA,neither),
                the i'th element in the sequence represents the i'th object in the model

    Returns:
        True iff there are more AnotB's than BnotA's elements in model
    """
    count_AnotB = 0
    count_BnotA = 0
    for symbol in model:
        if symbol == AnotB:
            count_AnotB += 1
        if symbol == BnotA:
            count_BnotA += 1
    return count_AnotB > count_BnotA


# make quantifier object 
M = Quantifier("M",
               mono="up", quan=True, cons=False,
               fn=models_most)


def get_all_quantifiers():
    """
    Returns: a list of all Quantifiers that have been created so far.
    """
    return [quant for quant in gc.get_objects()
            if isinstance(quant, Quantifier)]


def get_all_cons_quantifiers():
    """
    Returns: a list of all cons=True Quantifiers that have been created so far.
    """
    return [quant for quant in gc.get_objects()
            if isinstance(quant, Quantifier) and quant.cons]


all_minimal_pairs = [('at_least_4','at_least_6_or_at_most_2'), # monotonicity
                        ('at_most_3', 'at_least_6_or_at_most_2'), # monotonicity
                        ('at_least_3', 'first_3'), # quantity
                        ('at_least_3', 'last_3'), # quantity
                        ('not_all', 'not_only'), # conservativity
                        ('most', 'M') # conservativity
                        ]


all_cons_minimal_pairs = [('at_least_4','at_least_6_or_at_most_2'), # monotonicity
                            ('at_most_3', 'at_least_6_or_at_most_2'), # monotonicity
                            ('at_least_3', 'first_3'), # quantity
                            ('at_least_3', 'last_3') # quantity
                            ]



