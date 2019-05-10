# Complexity and learnability in the explanation of semantic universals of quantifiers

This repository accompanies the following paper:
* Iris van de Pol, Shane Steinert-Threlkeld, and Jakub Szymanik, "Complexity and learnability in the explanation of semantic universals of quantifiers", Proceedings of the 41st Annual Meeting of the Cognitive Science Society (CogSci 2019).

This project explores the complexity of quantifiers in the explanation of semantic universals, and compares these to earlier results on learnability in the explanation of semantic universals (see [Steinert-Threlkeld and Szymanik (in press)](https://semanticsarchive.net/Archive/mQ2Y2Y2Z/LearnabilitySemanticUniversals.pdf)). In particular, we compute the Lempel-Ziv (LZ) complexity (as an approximation to Kolmogorov complexity) of several quantifiers, some satisfying proposed universals (monotonicity, quantity, and conservativity) and some not. 

Our results indicate that the monotonicity universal can be explained by complexity while the conservativity universal cannot. For quantity we did not find a robust result. We also found that learnability and complexity pattern together in the monotonicity and conservativity cases that we consider, while that pattern is less robust in the quantity cases

This repository contains all the code needed to replicate the data reported in that paper.  It also contains the data and figures that are reported therein.

If you have any questions and/or want to extend the code-base and/or generate your own quantifier complexity data, feel free to get in touch!

## Getting Started

### Requirements

Python 3, [lempel_ziv_complexity](https://github.com/Naereen/Lempel-Ziv_Complexity), Pandas, plotnine, scipy, statsmodels

### Generating data

The core method for generating the data is `LZ_dist_quantifier_cumulative` in *gen_data.py*. 
The `LZ_dist_quantifier_cumulative` method accepts the following parameters:

- `quantifiers`: a list of quantifier objects
- `max_model_size`: an integer, determining the size (i.e., number of elements) of the quantifier models
- `num_model_orders`: an integer, determing the number of different model-order permutations
- `num_chars`: an integer, determining the number of characters used in the quantifier representation,
			(the number of subareas—--AB, AnotB, BnotA, and neither---of the quantifier model that are considered)
            needs to be either 2, 3, or 4
- `num_procs`: an integer, determining the number of processors used for parallel computing
- `out_file`: a string, determining the name of the csv file in which the results are stored
- `order_type`: either 'lex' (for lexicographic model orders) or 'rand' (for random model orders)

Example of a run from the command-line:

```
python gen_data.py --num_procs 4 --max_model_size 10 --jobtype cumulative --num_chars 4 --quantifiers all_quantifiers  --num_model_orders 12 --order_type lex --out_file output.csv

```

`LZ_dist_quantifier_cumulative` stores the results in an .csv file.

Note that in addition to the LZ complexity of quantifiers, this code also computes the complexity based on two other compression methods, gzip (based on LZ compression), and bzip2 (a block sorting compressor). Plots showing LZ and and bzip2 complexities can be found in this repository.

### Analyzing and plotting the data

The plots in the data report average complexity with confidence interval zones and were created with method `make_summary_complexity_per_size_plot` in *analysis.py*.
The `make_summary_complexity_per_size_plot` method accepts the following parameters:

- `data`: a pandas data frame
- `order`: an integer, the name of the model order of the quantifier representation
- `quants`: a list of quantifier names
- `out_file`: a string, name for a csv file (including .csv)
- `add_measure`: True or False, set to True when using facet_wrap	

The `make_summary_complexity_per_size_plot` method is called by the `plot` method. The plots are stored as .png files.

## Defining New Quantifiers

All of the quantifiers are defined in `quantifiers.py`.  The core is a verification method, named `models_quantifier_name(model)`.  This takes in a sequence of elements from (0,1,2,3) and outputs True or False.  This method implements the core semantics of a quantifier. 

The verification method is then wrapped inside a `Quantifier` object, which provides a name for the quantifier and assigns it certain semantic properties.  The `Quantifier` object calls the `models_` function in its `__call__` method, so that you can treat it as a method.

## Authors

* **Iris van de Pol** 
* **Shane Steinert-Threlkeld** 
* **Jakub Szymanik** 
