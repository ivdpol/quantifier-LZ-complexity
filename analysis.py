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

import argparse
import pandas as pd
from scipy import stats
import statsmodels.formula.api as sm
from plotnine import *
import quantifiers
import numpy as np


def plot(in_file, quants, max_model_size, out_file='test', orders=None,
            pairs_plot_per_order=False, summary_pair_plots = False,pairs_uniformity_plot=False,
            dist_plot=False, num_chars=None,order_type=None):
    """ Makes plots of complexity data, stores them in png files.

    Args:
        in_file: a string, the name of an existing csv file (without .csv at the end)
                    from which to load a pandas data frame
        quants: a list of quantifier names
        max_model_size: an integer
        out_file: a string, the name of a new csv file, not used in current code
        orders: integer, the num of model orders for which graphs are being made
        pairs_plot_per_order: True or False
        summary_pair_plots: True or False
        pairs_uniformity_plots: True or False
        dist_plot: True or False
        num_chars: an integer, 2, 3, or 4
        order_type: 'lex' or 'rand'
    """
    data = pd.read_csv(in_file + '.csv')
    data['order_num'] = data['order_num'].astype('category')

    if pairs_uniformity_plot:

        for quantifier_pair in quantifiers.all_minimal_pairs:
            for order in range(orders):
                make_uniformity_per_size_plot(data, order, quantifier_pair,
                                out_file = 'graphs/' + str(num_chars) + 'char/' + order_type +'/uniformity/'+ 
                                quantifier_pair[0] + '_&_' + quantifier_pair[1])

    if pairs_plot_per_order:

        for quantifier_pair in quantifiers.all_minimal_pairs:
            for order in range(orders):
                make_complexity_per_size_plot(data, order, quantifier_pair,
                                out_file = 'graphs/' + str(num_chars) + 'char/' + order_type +'/per_order/LZ/'+ 
                                quantifier_pair[0] + '_&_' + quantifier_pair[1] +'_order_' + str(order),
                                add_measure = True)

    if summary_pair_plots:

        for quantifier_pair in quantifiers.all_minimal_pairs:
                make_summary_complexity_per_size_plot(data, quantifier_pair,
                                out_file = 'graphs/' + str(num_chars) + 'char/' + order_type +'/summary/LZ/'+ 
                                quantifier_pair[0] + '_&_' + quantifier_pair[1],
                                add_measure = True)

    if dist_plot:
        make_combo_dist_plot(data, max_model_size, quants, 'dist_test')
 

def make_combo_dist_plot(data,max_model_size,quants=None,out_file='test_dist', field='C_LZ'):
    """ Makes a density plot of the distribution of complexity 
        of different (representations of) quantifiers over different model-orders, 
        for a given maximum model size.

    Args:
        data: a pandas data frame
        max_model_size: an integer
        quants: list of quantifier names
        out_file: a string, name for a csv file (including .csv)
        fiels: a string, the name of the complexity meassure to ben plotted
    """
    plot = (ggplot(
        data[(data['quantifier'].isin(quants)) &
                (data['max_model_size']== max_model_size)],
                    aes(field, colour='quantifier', fill='quantifier')) 
            + geom_density(alpha=0.2))

    plot.save(out_file)


def make_raw_plot(data,max_model_size,out_file='test_raw', field='C_LZ'):
    """ Makes line plots of the distribution of complexity 
        of different (representations of) quantifiers for a given maximum model size.

    Args:
        data: a pandas data frame
        max_model_size: an integer
        out_file: a string, name for a csv file (including .csv)
        fields: a string, the name of the complexity meassure to be plotted
    """
    data = data[data['max_model_size'] == max_model_size]
    plot = (ggplot(data,
                   aes(x='quantifier', y=field)) 
            + geom_point(aes(colour='order_num')) 
            + geom_line(aes(group='order_num', colour='order_type')) 
            + theme(axis_text_x  = element_text(angle = 300, hjust = 0))
            )

    plot.save(out_file)


def t_test_two_quantifiers(data, quant1, quant2, max_model_size,
                           min_model_size, measure='C_LZ'):
    """ Performs a t test bewteen complexity values for 2 quantifiers, 
        taking data points from different model sizes,
        and prints the results.

    Args:
        data: a pandas data frame
        quant1: a quantifier name
        quant2: a quantifier name
        max_model_size: an integer >= min_models_size
        min_model_size: an integer <= max_model_size
        measure: a string, name of a complexity measure
    """
    q1_data = data[(data['quantifier'] == quant1) &
                   (data['max_model_size'] <= max_model_size) &
                   (data['max_model_size'] >= min_model_size)]
    q2_data = data[(data['quantifier'] == quant2) &
                   (data['max_model_size'] <= max_model_size) &
                   (data['max_model_size'] >= min_model_size)]
    print()
    print('{}, {}:'.format(quant1, measure))
    print(q1_data[measure].describe())
    print()
    print('{}, {}:'.format(quant2, measure))
    print(q2_data[measure].describe())
    print()
    print('{} vs. {}:'.format(quant1, quant2))
    print(stats.ttest_rel(q1_data[measure], q2_data[measure]))

def t_test_min_pairs(data,max_model_size,cons=False,
                        measures = ('C_LZ','C_gzip','C_bzip')):
    """ Does multiple calls for t tests bewteen complexity values for pairs of quantifiers, 
        taking data points from different model sizes

    Args:
        data: a pandas data frame
        quant1: a quantifier name
        quant2: a quantifier name
        max_model_size: an integer >= min_models_size
        min_model_size: an integer <= max_model_size
        measure: a string, name of a complexity measure
    """
    if cons:
        pairs = quantifiers.cons_minimal_pairs       
    else:
        pairs = quantifiers.all_minimal_pairs

    for min_pair in pairs:
        quant1,quant2 = min_pair
        for measure in measures:
            t_test_two_quantifiers(data, quant1, quant2, max_model_size, measure=measure)          


def multiple_regression(data, dependent_var, independent_vars):
    """ Performs a simple multiple regression, and prints the results.

    Args:
        data: a pandas data frame
        dependent_var: a string, name of a complexity measure
        independent_var: a list of strings, names of quantifier properties
    """
    formula = dependent_var + ' ~ ' + ' + '.join(independent_vars)
    result = sm.ols(formula=formula, data=data).fit()
    print(result.summary())


def make_complexity_per_size_plot(data, order, quants=None, out_file='quantifier_pair', add_measure=False):
    """ Makes line plots of the complexity per quantifier in quants (y-axis), 
        per model size (x-axis)

    Args:
        data: a pandas data frame
        order: an integer, the name of the model order of the quantifier representation
        quants: a list of quantifier names
        out_file: a string, name for a csv file (including .csv)
        add_measure: True or False, set to True when using facet_wrap
    """
    data2 = data[data['quantifier'].isin(quants)]
    # "order" the quantifier column by the order specified in quants, to ensure
    # consistent coloring
    data2['quantifier'] = pd.Categorical(data2['quantifier'], ordered=True,
                                         categories=quants)
    if order > -1:
        data2 = data2[data2['order_num'] == order]

    # reshaping the data frame, adding a collumn called 'measure'
    # to plot different graphs per measure using facet_wrap
    if add_measure:
        data2 = pd.melt(data2, id_vars=['max_model_size','quantifier','order_num'],
                        #value_vars=['C_LZ','C_bzip','C_gzip'],
                        #var_name='measure', value_name='complexity')
                        value_vars=['C_LZ'],
                        var_name='measure', value_name='LZ complexity')

    plot = (ggplot(aes(x='max_model_size'), data=data2)            
            # + geom_line(aes(y='complexity', group='quantifier',
            #                colour='quantifier'))
            # + geom_point(aes(y='complexity', group='quantifier',
            #                colour='quantifier'))
            + geom_line(aes(y='LZ complexity', group='quantifier',
                           colour='quantifier'))
            + geom_point(aes(y='LZ complexity', group='quantifier',
                           colour='quantifier'))
            # + facet_wrap('measure', scales = "free")
            + scale_color_manual(values=['deepskyblue','red'])
            # take selection of x-axis:
            #+ xlim(3, 8)
           )

    plot.save(out_file)

def make_uniformity_per_size_plot(data, order,quants=None, out_file='quantifier_pair'):
    """ Makes line plots of the uniformity per quantifier in quants (y-axis), 
        per model size (x-axis)

    Args:
        data: a pandas data frame
        order: an integer, the name of the model order of the quantifier representation
        quants: a list of quantifier names
        out_file: a string, name for a csv file (including .csv)
    """
    data2 = data[data['quantifier'].isin(quants)]
    # "order" the quantifier column by the order specified in quants, to ensure
    # consistent coloring
    data2['quantifier'] = pd.Categorical(data2['quantifier'], ordered=True,
                                         categories=quants)
    if order > -1:
        data2 = data2[data2['order_num'] == order]

    plot = (ggplot(aes(x='max_model_size'), data=data2) 
            + geom_line(aes(y='uniformity', group='quantifier',
                           colour='quantifier'))
            + geom_point(aes(y='uniformity', group='quantifier',
                           colour='quantifier'))
            + scale_color_manual(values=['deepskyblue','red'])
           )

    plot.save(out_file)

def make_summary_complexity_per_size_plot(data, quants=None, out_file='quantifier_pair', add_measure=False):
    """ Makes line plots with confidence interval zone
        of the average complexity per quantifier in quants over the different model orders (y-axis), 
        per model size (x-axis)

    Args:
        data: a pandas data frame
        order: an integer, the name of the model order of the quantifier representation
        quants: a list of quantifier names
        out_file: a string, name for a csv file (including .csv)
        add_measure: True or False, set to True when using facet_wrap
    """
    pd.set_option('mode.chained_assignment',None)
    data2 = data[data['quantifier'].isin(quants)]
    # "order" the quantifier column by the order specified in quants, to ensure
    # consistent coloring
    data2['quantifier'] = pd.Categorical(data2['quantifier'], ordered=True,
                                         categories=quants)

    # reshaping the data frame, adding a collumn called 'measure'
    # to plot different graphs per measure using facet_wrap
    if add_measure:
        data2 = pd.melt(data2, id_vars=['max_model_size','quantifier','order_num'],
                        #value_vars=['C_LZ','C_bzip','C_gzip'],
                        #var_name='measure', value_name='complexity')
                        value_vars=['C_LZ'],
                        var_name='measure', value_name='LZ complexity')

    plot = (ggplot(aes(x='max_model_size'), data=data2)
            ## mean and confidence interval ##
            # + stat_summary(aes(y='complexity', group='quantifier',
            #              colour='quantifier'), geom='ribbon', fun_data='mean_cl_boot', fill='lightgray', size=0)
            # + stat_summary(aes(y='complexity', group='quantifier',
            #              colour='quantifier'), geom='line', fun_y=np.mean, size=0.75)
            + stat_summary(aes(y='LZ complexity', group='quantifier',
                         colour='quantifier'), geom='ribbon', fun_data='mean_cl_boot', fill='lightgray', size=0)
            + stat_summary(aes(y='LZ complexity', group='quantifier',
                         colour='quantifier'), geom='line', fun_y=np.mean, size=0.75)
           # + facet_wrap('measure', scales = "free", labeller=as_labeller({'C_bzip': 'C_Bzip2'}))
            + theme(panel_spacing=0.4,legend_position="top",legend_title=element_blank())
            + scale_x_continuous(name="maximum model size")
            + scale_color_manual(values=['deepskyblue','red'])
           )

    plot.save(out_file, width=2, height=4)


def concat_data(in_files,out_file):
    """ Concatinates pandas data frames stored in csv files into one big data frame, 
        and outputs it as a new csv file

    Args:
        in_files = a list of strings, names of exisiting csv files in current folder
        out_file = a string, name for new csv file
    """

    frames = [None for _ in range(len(in_files))]

    for idx in range(len(in_files)):
        frames[idx] = pd.read_csv(in_files[idx])
   
    out_data = pd.concat(frames, ignore_index=True)
    out_data.to_csv(out_file)

def describe_complexity_difference_per_pair(in_file='test'):
    """ Loads complexity data from csv file, pivots the data around quantifiers (i.e., makes quantifiers as columns),
        adds columns per quant pair with the difference between those quants, 
        outputs this to csv file (with only the column 'measure', all quantifier columns, and all difference columns),
        describes the statistics per measure of the difference between quantifier pairs,
        prints this description and saves it to a file per measure

    Args:
        in_file: a string, the name of an existing csv file (without .csv at the end)
                    from which to load a pandas data frame
    """

    data = pd.read_csv(in_file + '.csv')
    data = data[data['quantifier'] != 'at_least_6']
    data = pd.melt(data, id_vars=['max_model_size','quantifier','order_num'],
                            value_vars=['C_LZ','C_bzip','C_gzip'],
                            var_name='measure', value_name='complexity')
    data = data.pivot_table(index=['measure','max_model_size','order_num'], columns='quantifier', values='complexity')
    data.columns = data.columns.rename(None)
    data = data.reset_index()

    for quantifier_pair in quantifiers.all_minimal_pairs:
        data['diff_' + quantifier_pair[1] + '_&_' + 
            quantifier_pair[0]] = data[quantifier_pair[1]] - data[quantifier_pair[0]]

    data.drop(all_quantifiers,axis=1,inplace=True)
    data.to_csv(in_file + '_diff.csv')
    data.drop(['max_model_size', 'order_num'],axis=1,inplace=True) 

    for measure in ['C_LZ','C_bzip','C_gzip']:

        data1 = data[data['measure'] == measure]
        data_descr = data1.describe()
        data_descr.to_csv(in_file + '_diff_' + measure +'_descr.csv' )
        print(data_descr)


if __name__ == '__main__':

    all_cons_quantifiers = [quant.name for quant in quantifiers.get_all_cons_quantifiers()]
    all_quantifiers = [quant.name for quant in quantifiers.get_all_quantifiers()]

    parser = argparse.ArgumentParser()
    parser.add_argument('--in_file', type=str, default='data/4char/lex/C_4char_size10_lex_CogSci_1feb')
    parser.add_argument('--out_file', type=str, default='new_test_graphs')
    parser.add_argument('--orders', type=int, default=12)
    parser.add_argument('--quants', type=str, default='all_quantifiers')
    parser.add_argument('--max_model_size',type=int, default=10)
    parser.add_argument('--num_chars', type=int, default = 4)
    parser.add_argument('--order_type', type=str, default = 'lex')
    args = parser.parse_args()

    plot(args.in_file,
            quantifiers.all_minimal_pairs,
            args.max_model_size,
            args.out_file,
            args.orders,
            pairs_plot_per_order=False,
            summary_pair_plots=True,
            dist_plot=False,
            pairs_uniformity_plot=False,
            num_chars = args.num_chars,
            order_type = args.order_type)



