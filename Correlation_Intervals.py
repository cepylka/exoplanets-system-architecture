import numpy
from tabulate import tabulate
import pandas
import warnings
import matplotlib
warnings.filterwarnings("ignore", category=UserWarning)
pandas.options.mode.use_inf_as_na = True
import math

pandas.options.mode.chained_assignment = None  # default='warn'

import seaborn

import matplotlib.pyplot as plt
from astropy import units
from astropy import constants
from astropy.visualization import quantity_support
from matplotlib_inline.backend_inline import set_matplotlib_formats
from matplotlib import cm
import matplotlib.colors
import matplotlib.ticker
from scipy import stats # For in-built method to get PCC
import scipy.optimize
from scipy import interpolate
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.integrate import trapz

from scipy.stats import norm

from matplotlib.colors import ListedColormap
from matplotlib.colors import LogNorm
from matplotlib.colors import LinearSegmentedColormap
from collections import OrderedDict
from matplotlib.lines import Line2D
import itertools
from collections import OrderedDict
from functools import partial
import matplotlib.ticker as mticker
from cycler import cycler
mass_coeff = constants.M_jup / constants.M_earth
radius_coeff = constants.R_jup / constants.R_earth
pandas.options.mode.use_inf_as_na = True

quantity_support()
set_matplotlib_formats('svg')
plt.rc('legend', frameon=False)
plt.rc('figure', figsize=(10, 10))
plt.rc('font', size=12)


workingTableExoplanets1 = pandas.read_pickle("./data/MRP_data_sets.pkl")

workingTableExoplanets1["mass_ratio"] = workingTableExoplanets1["mass-2"]/workingTableExoplanets1["mass-1"]
workingTableExoplanets1 = workingTableExoplanets1[~workingTableExoplanets1['mass_detection_type'].isin(['Theoretical'])]
workingTableExoplanets1 = workingTableExoplanets1.query("`mass-1`.notna() & `mass-2`.notna()")
workingTableExoplanets1 = workingTableExoplanets1.query("`mass_error_min-1`.notna() & `mass_error_max-1`.notna() & `mass_error_min-2`.notna() & `mass_error_max-2`")
workingTableExoplanets1 = workingTableExoplanets1.query("`pl_massjlim`.isnull()")
workingTableExoplanets1 = workingTableExoplanets1.query("`ima_flag`.isnull()")
mass_median = numpy.median(workingTableExoplanets1["mass_ratio"])
mass_std = numpy.std(workingTableExoplanets1["mass_ratio"])
print(f"Mass median = {mass_median}")
print(f"Mass std = {mass_std}")
a = numpy.linspace(0.02,0.5, 300)
interval_value = None

for an in a:
    # while True:
    m1 = mass_median - an * mass_median #mass_std/20
    m2 = mass_median + an * mass_median #mass_std/20
    # print(an)

    Exoplanets0 = workingTableExoplanets1.query("`mass_ratio` < @m2 & `mass_ratio` > @m1")


    pearson_coef, p_value = stats.pearsonr(Exoplanets0["mass-1"] * mass_coeff, Exoplanets0["mass-2"] * mass_coeff)
    if pearson_coef > 0.95 and pearson_coef < 0.96:
        interval_value = an
        break
if not interval_value:
    raise SystemExit("[ERROR] interval_value is None")
print(f"Mass correlated sample size = {len(Exoplanets0)}")
granule_1u = numpy.unique(Exoplanets0["granule_uid-1"])
granule_2u = numpy.unique(Exoplanets0["granule_uid-2"])
granule = list(granule_1u)
granule.extend(x for x in granule_2u if x not in granule)
print(f"mass total planets:{len(granule)}")
sys_len =len(numpy.unique(Exoplanets0["star_name"]))
print(f"mass total systems:{sys_len}")
Exoplanets0.to_pickle("./data/corellated_masses_new.pkl")
z = numpy.arctanh(pearson_coef)

sigma = (1/((len(Exoplanets0)-3)**0.5))


cint = z + numpy.array([-1, 1]) * sigma * stats.norm.ppf((1+0.95)/2)


confidence_intervals =  numpy.tanh(cint)
print("Pearson Correlation Coefficient for mass data: ", pearson_coef, "and P-value of:", p_value, ", and the intervals +/-", interval_value, "confidence_intervals =",confidence_intervals) # Results


workingTableExoplanets2 = pandas.read_pickle("./data/MRP_data_sets.pkl")

workingTableExoplanets2 = workingTableExoplanets2.query("`radius-1`.notna() & `radius-2`.notna()")
workingTableExoplanets2 = workingTableExoplanets2.query("`radius_error_min-1`.notna() & `radius_error_max-1`.notna() & `radius_error_min-2`.notna() & `radius_error_max-2`")
workingTableExoplanets2 = workingTableExoplanets2.query("`pl_radjlim`.isnull()")
workingTableExoplanets2 = workingTableExoplanets2.query("`ima_flag`.isnull()")

workingTableExoplanets2["radius_ratio"] = workingTableExoplanets2["radius-2"]/workingTableExoplanets2["radius-1"]

radius_median = numpy.median(workingTableExoplanets2["radius_ratio"])

radius_std = numpy.std(workingTableExoplanets2["radius_ratio"])

print(f"Radius median = {radius_median}")
print(f"Radius std = {radius_std}")
interval_value = None
for an in a:

    r1 = radius_median - an * radius_median
    r2 = radius_median + an * radius_median


    Exoplanets1 = workingTableExoplanets2.query("`radius_ratio` < @r2 & `radius_ratio` > @r1")
    pearson_coef1, p_value1 = stats.pearsonr(Exoplanets1["radius-1"] * radius_coeff, Exoplanets1["radius-2"] * radius_coeff) #define the columns to perform calculations on

    if pearson_coef1 > 0.95 and pearson_coef1 < 0.96:
        interval_value = an
        break
if not interval_value:
    raise SystemExit("[ERROR] interval_value is None")
print(f"Radius correlated sample size = {len(Exoplanets1)}")
granule_1u = numpy.unique(Exoplanets1["granule_uid-1"])
granule_2u = numpy.unique(Exoplanets1["granule_uid-2"])
granule = list(granule_1u)
granule.extend(x for x in granule_2u if x not in granule)
print(f"radius total planets:{len(granule)}")
sys_len =len(numpy.unique(Exoplanets1["star_name"]))
print(f"radius total systems:{sys_len}")
Exoplanets1.to_pickle("./data/corellated_radii_new.pkl")
z = numpy.arctanh(pearson_coef1)

sigma = (1/((len(Exoplanets1)-3)**0.5))


cint = z + numpy.array([-1, 1]) * sigma * stats.norm.ppf((1+0.95)/2)


confidence_intervals1 =  numpy.tanh(cint)
print("Pearson Correlation Coefficient for radius data: ", pearson_coef1, "and a P-value of:", p_value1, ", and the intervals +/-", interval_value, "confidence_intervals =",confidence_intervals1) # Results

workingTableExoplanets4 = pandas.read_pickle("./data/MRP_data_sets.pkl")

workingTableExoplanets4['density-1'] = (workingTableExoplanets4["mass-1"] * mass_coeff * constants.M_earth) * 3 /(4 * numpy.pi*numpy.power(workingTableExoplanets4["radius-1"] * radius_coeff * constants.R_earth,3))
workingTableExoplanets4['density-2'] = (workingTableExoplanets4["mass-2"] * mass_coeff * constants.M_earth) * 3 /(4 * numpy.pi*numpy.power(workingTableExoplanets4["radius-2"] * radius_coeff * constants.R_earth,3))

workingTableExoplanets4['density_error_min-1'] = numpy.sqrt((workingTableExoplanets4["mass_error_min-1"]/workingTableExoplanets4["mass-1"]) ** 2 + (3 * workingTableExoplanets4["radius_error_min-1"]/workingTableExoplanets4["radius-1"]) ** 2) * workingTableExoplanets4['density-1']# assigned to a column
workingTableExoplanets4['density_error_max-1'] = numpy.sqrt((workingTableExoplanets4["mass_error_max-1"]/workingTableExoplanets4["mass-1"]) ** 2 + (3 * workingTableExoplanets4["radius_error_max-1"]/workingTableExoplanets4["radius-1"]) ** 2) * workingTableExoplanets4['density-1']
workingTableExoplanets4['density_error_min-2'] = numpy.sqrt((workingTableExoplanets4["mass_error_min-2"]/workingTableExoplanets4["mass-2"]) ** 2 + (3 * workingTableExoplanets4["radius_error_min-2"]/workingTableExoplanets4["radius-2"]) ** 2) * workingTableExoplanets4['density-2']# assigned to a column
workingTableExoplanets4['density_error_max-2'] = numpy.sqrt((workingTableExoplanets4["mass_error_max-2"]/workingTableExoplanets4["mass-2"]) ** 2 + (3 * workingTableExoplanets4["radius_error_max-2"]/workingTableExoplanets4["radius-2"]) ** 2) * workingTableExoplanets4['density-2']

workingTableExoplanets4 = workingTableExoplanets4.query("`density-1`.notna() & `density-2`.notna()")
workingTableExoplanets4 = workingTableExoplanets4.query("`density_error_min-1`.notna() & `density_error_max-1`.notna() & `density_error_min-2`.notna() & `density_error_max-2`")
workingTableExoplanets4 = workingTableExoplanets4.query("`pl_radjlim`.isnull() & `pl_massjlim`.isnull() ")
workingTableExoplanets4 = workingTableExoplanets4.query("`ima_flag`.isnull()")
print(len(workingTableExoplanets4))
workingTableExoplanets4["density_ratio"] = workingTableExoplanets4["density-2"]/workingTableExoplanets4["density-1"]
density_median = numpy.mean(workingTableExoplanets4["density_ratio"])

density_std = numpy.std(workingTableExoplanets4["density_ratio"])

print(f"Density median = {density_median}")
print(f"Density std = {density_std}")

interval_value = None
for an in a:
    d1 = density_median - an * density_median
    d2 = density_median + an * density_median


    Exoplanets4 = workingTableExoplanets4.query("`density_ratio` < @d2 & `density_ratio` > @d1")
    pearson_coef2, p_value2 = stats.pearsonr(Exoplanets4["density-1"], Exoplanets4["density-2"]) #define the columns to perform calculations on
    if pearson_coef2 > 0.95 and pearson_coef2 < 0.96:
        interval_value = an
        print(interval_value)
        break
if not interval_value:
    raise SystemExit("[ERROR] interval_value is None")
print(f"Density correlated sample size= {len(Exoplanets4)}")
granule_1u = numpy.unique(Exoplanets4["granule_uid-1"])
granule_2u = numpy.unique(Exoplanets4["granule_uid-2"])
granule = list(granule_1u)
granule.extend(x for x in granule_2u if x not in granule)
print(f"density total planets:{len(granule)}")
sys_len =len(numpy.unique(Exoplanets4["star_name"]))
print(f"density total systems:{sys_len}")
Exoplanets4.to_pickle("./data/corellated_density_new.pkl")
z = numpy.arctanh(pearson_coef2)

sigma = (1/((len(Exoplanets4)-3)**0.5))


cint = z + numpy.array([-1, 1]) * sigma * stats.norm.ppf((1+0.95)/2)


confidence_intervals2 =  numpy.tanh(cint)

print("Pearson Correlation Coefficient for density data: ", pearson_coef2, "and P-value of:", p_value2, ", the intervals +/-", interval_value, "confidence_intervals =",confidence_intervals2) # Results


workingTableExoplanets3 = pandas.read_pickle("./data/triples_MR.pkl")

workingTableExoplanets3 = workingTableExoplanets3.query("`period-1`.notna() & `period-2`.notna()& `period-3`.notna()")
workingTableExoplanets3 = workingTableExoplanets3.query("`period_error_min-1`.notna() & `period_error_max-1`.notna() & `period_error_min-2`.notna() & `period_error_max-2` & `period_error_min-3`.notna() & `period_error_max-3`")
workingTableExoplanets3 = workingTableExoplanets3.query("`ima_flag`.isnull()")
workingTableExoplanets3["period_ratio_1"] = workingTableExoplanets3["period-2"]/workingTableExoplanets3["period-1"]
workingTableExoplanets3["period_ratio_2"] = workingTableExoplanets3["period-3"]/workingTableExoplanets3["period-2"]
workingTableExoplanets3["period-ratioratio"] = workingTableExoplanets3["period_ratio_2"] /workingTableExoplanets3["period_ratio_1"]
periodratio_median = numpy.median(workingTableExoplanets3["period-ratioratio"])
print(f"Period ratio median = {periodratio_median}")
periodratio_std = numpy.std(workingTableExoplanets3["period-ratioratio"])
print(f"Period ratio std = {periodratio_std}")

for an in a:
    p1 = periodratio_median - an * periodratio_median
    p2 = periodratio_median + an * periodratio_median


    Exoplanets2 = workingTableExoplanets3.query("`period-ratioratio` < @p2 & `period-ratioratio` > @p1")
    pearson_coef3, p_value3 = stats.pearsonr(Exoplanets2["period_ratio_1"] * radius_coeff, Exoplanets2["period_ratio_2"] * radius_coeff) #define the columns to perform calculations on
    if pearson_coef3 > 0.95 and pearson_coef3 < 0.96:
        interval_value = an
        break
if not interval_value:
    raise SystemExit("[ERROR] interval_value is None")
print(f"Period ratio correlated sample size ={len(Exoplanets2)}")
granule_1u = numpy.unique(Exoplanets2["granule_uid-1"])
granule_2u = numpy.unique(Exoplanets2["granule_uid-2"])
granule_3u = numpy.unique(Exoplanets2["granule_uid-3"])
granule = list(granule_1u)
granule.extend(x for x in granule_2u if x not in granule)
granule.extend(x for x in granule_3u if x not in granule)
print(f"period ratio  total planets:{len(granule)}")
sys_len =len(numpy.unique(Exoplanets2["star_name"]))
print(f"period ratio total systems:{sys_len}")
Exoplanets2.to_pickle("./data/corellated_period_ratio_new.pkl")
z = numpy.arctanh(pearson_coef3)

sigma = (1/((len(Exoplanets2)-3)**0.5))


cint = z + numpy.array([-1, 1]) * sigma * stats.norm.ppf((1+0.95)/2)


confidence_intervals3 =  numpy.tanh(cint)
print("Pearson Correlation Coefficient for period ratio data: ", pearson_coef3, "and a P-value of:", p_value3, ", and the intervals +/-", interval_value, "confidence_intervals =",confidence_intervals3) # Results # Results

