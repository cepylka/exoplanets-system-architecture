# Exoplanets system architecture analysis

<!-- MarkdownTOC -->

- Preparing the sample
- Analyses of the sample
    - System architecture: planet parameters similarities
- Retrieving stellar parameters from GAIA Archive

<!-- /MarkdownTOC -->

## Preparing the sample

For retrieving information from databases initially we need a list of potentially interesting systems, which can be formed using the `1-get-systems-with-more-than-planets.py` script.

This script allows to retrieve star names for the systems, where at least two planets have both mass and radius defined, by fetching parameter `star_name` from NASA database.

``` sh
$ python ./1-get-systems-with-more-than-planets.py
```

The results will be output to `systems.txt` file in the current directory, which then could be edited manually, if you already know what system you want to exclude to save some processing time.

For retrieving the partial information about the system in interest, you can run `2-task-manager-mass-radius.sh`. It will work through the planetary system list in `systems.txt` using `2-mass-radius-only-from-NASA.py` script and will write a pickle file for each system to `./data/star_pickles/`:

``` sh
$ chmod +x ./2-task-manager-mass-radius.sh
$ ./2-task-manager-mass-radius.sh
```

The next script `3-copy-those-with-mass-radius.py` goes through produced pickles, selects those that have at least 2 planets with known both mass and radius values and copies them to `./data/star_pickles/minimum2withmassandradius/`. To run this script you need to specify the folder with pickles from the previous step:

``` sh
$ python ./3-copy-those-with-mass-radius.py ./data/star_pickles
```

To collect all other parameters you'd need a script `4-get-all-other-parameters.py`. You can run `task-manager.sh` for this task with the name of the script and input folder from which to take the pickles (*in this case `./data/minimum2withmassandradius/`*):

``` sh
$ ./task-manager.sh -i ./data/star_pickles/minimum2withmassandradius -s ./4-get-all-other-parameters.py
```

The output will be saved to `./data/star_pickles_enriched/`.

In order to further investigate the data you need to concatenate all chosen system files into one, which certain Jupyter notebooks will be working with. This combined file in `.pkl` format will also be needed for creating working files for internal systems structure analyses. You can do the concatenation by running the `5-merge-edges.py` script (*after putting the `.pkl` files of interest into the `./merging/` folder*):

``` sh
$ mkdir ./merging
$ cp ./data/star_pickles_enriched/*.pkl ./merging/
$ python ./5-merge-edges.py
```

Results will be saved to `./data/all_my_systems.pkl`.

You can additionally check if at least two planets having masses and radii determined are consecutive. The code for that is in the Jupyter notebook `./jupyter/two-consecutive-with-mass-and-radius.ipynb`.

In aforementioned folders (*`./data`, `./merging`*) I put example planetary system `.pkl` files for checking how the scripts work. After executing the routine from the begining, you will obtain such files of the new systems that meet your selection criteria. You can also make your own file tree for your convenience, but do not forget to change the file paths if you do.

## Analyses of the sample

### System architecture: planet parameters similarities

For analysing the "peas in the pod" tendency among masses and radii of adjacent planets in the system, you can create a new `.pkl` file from existing sample file `all_my_systems.pkl` by using `6-pairs_MRP_stellar.py`, which you can run in the terminal.

``` sh
$ python ./6-pairs_MRP_stellar.py
```

It will write a new file `./data/MRP_data_sets.pkl`. This file contains paired data for two adjacent planets in the system as one row in the pandas data frame.

For analyses of similarity in parameters of adjacent planet pairs and its depending on stellar parameters you can use the following script: `MRD_adjacent_planets.py`.

``` sh
$ python ./MRD_adjacent_planets.py
```

First, this script will create general plots for parameters: mass, radius, and density, and for two populations: P<sub>i</sub> for the inner planet in the pair of adjacent planets, P<sub>i+1</sub> for the outer planet in this pair, where P is one of the parameters and i=1,2 for first and second planet in the system by ascending of their semi-major axis. Note, that one system can produce multiple pairs, if it has more that 2 planets with well-defined parameters. The Pearson correlation test will be conducted and R- and P-values will be calculated, and the result along with the ordering in pairs calculations will be output to stdout.

Then, the script calculates the Pearson coefficient distributions for the Monte Carlo random uniform sampling populations from error intervals (10<sup>5</sup> attempts), and resulting figures will represent R- and P-values distributions, the normal distribution function fitting. And the mean and median coefficient values and the standard deviation $\sigma$ value will be printed in every figure.

In the third part, the script will make "moving window" test for different stellar parameters for sub-samples of x=40 data points with n=2 step. It will produce graphs of R- and P-value dynamics in stellar mass, radius, metallicity, T<sub>eff</sub> and age ranges.


For analysing periods ratio, you can create a new `.pkl` file from existing sample file `all_my_systems.pkl` by using `triples_MRPD.py`, which you can run in the terminal.

``` sh
$ python ./triples_MRPD.py
```

It will write a new file `triples_MR.pkl` in the `data` folder. The file consists now data for three adjacent planets in the system as one row in the pandas data frame.

For analyses of similarity in parameters of adjacent planet triples, and for assessing how these similarity trends depend on stellar parameters, you can use the following script: `P_adjacent_planets.py`.

``` sh
$ python ./P_adjacent_planets.py
```

Analogously with the script for masses, radii and densities (`MRD_adjacent_planets.py`), the `P_adjacent_planets.py` script calculates R- and P-values and output to stdout, creates general plots for period ratios, and for two populations: P<sub>i+1</sub>/P<sub>i</sub> for the inner planet pair of adjacent planets, P<sub>i+2</sub>/P<sub>i+1</sub>  for the outer planet pair, where P is the orbital period and i=1,2,3 for first, second and third planet in the system by ascending of their semi-major axis. Note, that one system can produce multiple triples, if it has more that 3 planets with well-defined parameters.

Second, the script calculates the Pearson coefficient distributions from random uniform simulations from error intervals (10<sup>5</sup> attempts), and plots figures of R- and P-values distributions, the normal distribution function fitting along with the mean and median coefficient values and the standard deviation $\sigma$ value.

After that, the script will make "moving window" test for different stellar parameters for sub-samples of x=40 data points with n=2 step. It will produce graphs of R- and P-value dynamics in stellar mass, radius, metallicity, T<sub>eff</sub> and age ranges.

## Retrieving stellar parameters from GAIA Archive

<!-- ### System architecture: planet parameters similarities -->
The python module `uio._task` \href{https://github.com/retifrav/uio-exoplanet-group}{
uio-exoplanet-group} is used for retrieving certain stellar parameters from GAIA DR3 release [Fouesneau et al., 2022](#Fouesneau).

```
from uio._tasks import reconfirming_stellar_parameters
```

In this case, the module works with the exisitng file `all_my_systems.pkl` and, from the table `gaiadr3.astrophysical_parameters`, it retrieves `radius_flame` and `radius_gspphot` parameters with corresponding upper and lower limits.

```
tbl = reconfirming_stellar_parameters.lookForParametersInGaia(
    # "./data/all_my_systems.pkl",
    "./data/your_new_file.pkl",
    "gaiadr3.astrophysical_parameters",
    [
        "radius_flame",
        "radius_flame_lower",
        "radius_flame_upper",
        "radius_gspphot",
        "radius_gspphot_lower",
        "radius_gspphot_upper"
    ]
)
```
For explanation of this parameters we refer to [Andrae et al., 2022](#Andrae). Then the module writes the retrieved data in the same table, which will be saved in `your_new_file.pkl`.


References:

<a name="Andrae"></a>Andrae, R. et al. (2022) ‘Gaia Data Release 3: Analysis of the Gaia BP/RP spectra using the General Stellar Parameterizer from Photometry’, Astronomy & Astrophysics [Preprint]. Available at: https://doi.org/10.1051/0004-6361/202243462.

<a name="Fouesneau"></a>Fouesneau, M. et al. (2022) ‘Gaia Data Release 3. Apsis II: Stellar parameters’, Astronomy & Astrophysics [Preprint]. Available at: https://doi.org/10.1051/0004-6361/202243919.

<a name="Kopparapu"></a> Kopparapu, R. K., Ramirez, R. M., SchottelKotte, J., Kasting, J. F., Domagal-Goldman, S. and Eymet, V. (2014). "Habitable zones around main-sequence stars: dependence
on planetary mass". In: The Astrophysical Journal Letters vol. 787, no. 2, p. L29.

<a name="Zeng"></a> Zeng, L., Jacobsen, S. B., Sasselov, D. D., Petaev, M. I., Vanderburg, A., Lopez-Morales, M., Perez-Mercader, J., Mattsson, T. R., Li, G., Heising, M. Z. et al. (2019). "Growth model interpretation of planet size distribution". In: Proceedings of the National Academy of Sciences vol. 116, no. 20, pp. 9723–9728.
