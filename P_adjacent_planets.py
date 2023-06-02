import pyvo
import numpy
from tabulate import tabulate
import pandas
import warnings
import matplotlib
warnings.filterwarnings("ignore", category=UserWarning)
warnings.simplefilter(action="ignore", category=FutureWarning)

pandas.options.mode.use_inf_as_na = True
import math
service = pyvo.dal.TAPService("http://voparis-tap-planeto.obspm.fr/tap")
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
from collections import OrderedDict
from matplotlib.lines import Line2D
import itertools
from collections import OrderedDict
from functools import partial
import matplotlib.ticker as mticker
from cycler import cycler
mass_coeff = constants.M_jup / constants.M_earth
radius_coeff = constants.R_jup / constants.R_earth

quantity_support()
set_matplotlib_formats('svg')
plt.rc('legend', frameon=False)
plt.rc('figure', figsize=(10, 10))
plt.rc('font', size=12)

workingTableExoplanets = pandas.read_pickle("./data/MRP_data_sets_triples.pkl")


workingTableExoplanets = workingTableExoplanets.query("`ima_flag`.isna()")
workingTableExoplanets = workingTableExoplanets[~workingTableExoplanets['star_name'].isin(['KOI-55'])] # marked as controversial in NASA

workingTableExoplanets = workingTableExoplanets.query("`period-1`.notna() & `period-2`.notna() & `period-3`.notna()")
workingTableExoplanets = workingTableExoplanets.query("`period_error_min-1`.notna() & `period_error_max-1`.notna() & `period_error_min-2`.notna() & `period_error_max-2`& `period_error_min-3`.notna() & `period_error_max-3`")

workingTableExoplanets['periodratio-1'] = workingTableExoplanets["period-2"] / workingTableExoplanets["period-1"]
workingTableExoplanets['periodratio-2'] = workingTableExoplanets["period-3"] / workingTableExoplanets["period-2"]

workingTableExoplanets['periodratio_error_min-1'] = numpy.sqrt((workingTableExoplanets["period_error_min-1"]/workingTableExoplanets["period-1"]) ** 2 + (workingTableExoplanets["period_error_min-2"]/workingTableExoplanets["period-2"]) ** 2) * workingTableExoplanets['periodratio-1']
workingTableExoplanets['periodratio_error_max-1'] = numpy.sqrt((workingTableExoplanets["period_error_max-1"]/workingTableExoplanets["period-1"]) ** 2 + (workingTableExoplanets["period_error_max-2"]/workingTableExoplanets["period-2"]) ** 2) * workingTableExoplanets['periodratio-1']
workingTableExoplanets['periodratio_error_min-2'] = numpy.sqrt((workingTableExoplanets["period_error_min-3"]/workingTableExoplanets["period-3"]) ** 2 + (workingTableExoplanets["period_error_min-2"]/workingTableExoplanets["period-2"]) ** 2) * workingTableExoplanets['periodratio-2']
workingTableExoplanets['periodratio_error_max-2'] = numpy.sqrt((workingTableExoplanets["period_error_min-3"]/workingTableExoplanets["period-3"]) ** 2 + (workingTableExoplanets["period_error_min-2"]/workingTableExoplanets["period-2"]) ** 2) * workingTableExoplanets['periodratio-2']

star_age = workingTableExoplanets["star_age"]
star_teff = workingTableExoplanets["star_teff"]
star_mass = workingTableExoplanets["star_mass"]
star_radius = workingTableExoplanets["star_radius"]
star_metallicity = workingTableExoplanets["star_metallicity"]
star_rotp = workingTableExoplanets["st_rotp"]

listForstellarParam = [
    "star_mass",
    "star_radius",
    "star_teff",
    "star_metallicity",
    "star_age"
]

trisuchnosti = [
"periodratio"
]
for suchnost in trisuchnosti:
    currenttable = workingTableExoplanets

    print(suchnost)
    print("[1] Making charts #1...")
    pearson_coef= []
    p_value  = []

    planetstodrop = []

    for index, row in workingTableExoplanets.iterrows():
        if (abs(workingTableExoplanets.at[index,"period_error_min-1"]) + workingTableExoplanets.at[index,"period_error_max-1"])/workingTableExoplanets.at[index,"period-1"] >= 2.0:
            planetstodrop.append(index)
        if (abs(workingTableExoplanets.at[index,"period_error_min-2"]) + workingTableExoplanets.at[index,"period_error_max-2"])/workingTableExoplanets.at[index,"period-2"] >= 2.0:
            planetstodrop.append(index)
        if (abs(workingTableExoplanets.at[index,"period_error_min-3"]) + workingTableExoplanets.at[index,"period_error_max-3"])/workingTableExoplanets.at[index,"period-2"] >= 2.0:
            planetstodrop.append(index)

    currenttable = currenttable.drop(index=planetstodrop)
    granule_1u = numpy.unique(currenttable["granule_uid-1"])
    granule_2u = numpy.unique(currenttable["granule_uid-2"])
    granule_3u = numpy.unique(currenttable["granule_uid-3"])
    granule = list(granule_1u)
    granule.extend(x for x in granule_2u if x not in granule)
    granule.extend(x for x in granule_3u if x not in granule)

    print(f"{suchnost} total planets:{len(granule)}")
    sys_len =len(numpy.unique(currenttable["star_name"]))
    print(f"{suchnost} total systems:{sys_len}")
    print(f"Working with {suchnost}, {len(currenttable)} pairs")

    resultsExoplanetsm2 = currenttable.query(f"`{suchnost}-2` > `{suchnost}-1`")
    resultsExoplanetsm1 = currenttable.query(f"`{suchnost}-2` < `{suchnost}-1`")
    resultsExoplanetsmeq = currenttable.query(f"`{suchnost}-2` == `{suchnost}-1`")


    suchnost_1 = currenttable[f"{suchnost}-1"]
    suchnost_2 = currenttable[f"{suchnost}-2"]
    # errror_frac = (abs(currenttable[f"{suchnost}_error_min-1"]) +currenttable[f"{suchnost}_error_max-1"])/suchnost_1 + (abs(currenttable[f"{suchnost}_error_min-2"]) +currenttable[f"{suchnost}_error_max-2"])/suchnost_2
    xmin= currenttable[[f"{suchnost}-1", f"{suchnost}-2"]].min(axis=1).min(axis=0) * 0.5
    ymax= currenttable[[f"{suchnost}-1", f"{suchnost}-2"]].max(axis=1).max(axis=0) * 1.5
    xerr = numpy.nan_to_num(currenttable[[f"{suchnost}_error_min-1", f"{suchnost}_error_max-1"]].to_numpy().T, posinf=0.)
    yerr = numpy.nan_to_num(currenttable[[f"{suchnost}_error_min-2", f"{suchnost}_error_max-2"]].to_numpy().T, posinf=0.)

    # checking for ordering in pairs
    print(f'Upper side of the line, planets:', len(resultsExoplanetsm2), ', part:', numpy.round(len(resultsExoplanetsm2)/ len(currenttable), 3))
    print(f'Lower side of the line, planets:', len(resultsExoplanetsm1), ', part:', numpy.round(len(resultsExoplanetsm1)/ len(currenttable), 3))
    print(f'On the line, planets:', len(resultsExoplanetsmeq))

    pearson_coef, p_value = stats.pearsonr(suchnost_1, suchnost_2) #define the columns to perform calculations on
    print(f"Pearson Correlation Coefficient for {suchnost} R-value: {pearson_coef} and P-value: {p_value}") # Results
    spearman_coef, s_value = stats.spearmanr(suchnost_1, suchnost_2) #define the columns to perform calculations on
    print("Spearman Correlation Coefficient for raw data: ", spearman_coef, "and a P-value of:", s_value) # Resultsc

    z = numpy.arctanh(pearson_coef)

    sigma = (1/((len(currenttable)-3)**0.5))


    cint = z + numpy.array([-1, 1]) * sigma * stats.norm.ppf((1+0.95)/2)


    confidence_intervals =  numpy.tanh(cint)
    print(f"Pearson Correlation Coefficient confidence intervals: {confidence_intervals}")

    cmap = seaborn.cubehelix_palette(rot=-.3, as_cmap=True)
    seaborn.set_theme(style="white")

    fig, ax = plt.subplots(figsize=(8.5, 8.5))
    ax.grid(True)
    plt.rcParams['figure.figsize']=(6,6)

    seaborn.scatterplot(
        x=suchnost_1, y=suchnost_2,
        palette=cmap, legend=False, ax=ax
    )

    ax.set(xscale="log", yscale="log")
    ax.grid(True)
    ax.set(ylim=(xmin, ymax))
    ax.set(xlim=(xmin, ymax))
    ax.set_ylabel(f'P$_{{i+2}}$/P$_{{i+1}}$', fontsize=14)
    ax.set_xlabel(f'P$_{{i+1}}$/P$_i$', fontsize=14)
    # g._legend.remove()
    x = numpy.linspace(xmin*0.9, ymax*1.1, 2000)
    y = x
    plt.plot(x, y, linewidth=0.8, linestyle='--', color='k')

    plt.errorbar(suchnost_1, suchnost_2, xerr=numpy.abs(xerr), yerr=numpy.abs(yerr), ls='none', fmt='0.8', ecolor='tab:gray', elinewidth=0.8, capsize=None, barsabove=True, zorder=0)
    plt.scatter(suchnost_1, suchnost_2, marker="o", facecolor='tab:blue', zorder=0,)# label=r'R:{{(pearson_coef).round(3)}}, p: {{p_value.round(2)}}x10$^{{{round(p_value)}}}$')

    solarsystemTable = pandas.read_pickle("./data/solarsystemE.pkl")


    for i in range(len(solarsystemTable[f"periodE"])-2):
        if i == 1:
            plt.plot(solarsystemTable[f"periodE"][i+2]/solarsystemTable[f"periodE"][i+1], solarsystemTable[f"periodE"][i+1]/solarsystemTable[f"periodE"][i], 'ro', ms=6, zorder=1, label="Solar system")
        else:
            plt.plot(solarsystemTable[f"periodE"][i+2]/solarsystemTable[f"periodE"][i+1], solarsystemTable[f"periodE"][i+1]/solarsystemTable[f"periodE"][i], 'ro', ms=6, zorder=1)

    f = matplotlib.ticker.ScalarFormatter(useMathText=True)
    f.set_powerlimits((-3,3))
    ax5 = fig.add_axes([0.7, 0.05, 0.25, 0.25])#, sharex=ax)
    dist = seaborn.ecdfplot(x=numpy.log10(suchnost_2/suchnost_1),stat='count', ax=ax5, label=f"{suchnost}")
    ax5.axhline(y=len(numpy.log10(suchnost_2/suchnost_1))/2.,color='r', linestyle="--",linewidth=0.8, )
    ax5.axvline(x=0.,color='k', linestyle="--",linewidth=0.8, )
    ax5.set_xlabel(f'log$_{{10}}$(P$_{{i3/i2}}$/P$_{{i2/i1}}$)',fontsize=14)
    ax.text(23, 83,f'$R$:{(pearson_coef).round(3)}, $p$: {p_value.round(3)}', fontsize=14,weight='bold')
    ax5.set_ylabel(f'CDF', fontsize=14)
    ax5.yaxis.tick_right()
    ax5.yaxis.set_label_position('right')
    ax5.grid(False)
    # ax.legend(frameon=False)
    plt.savefig(f"all_my_systems_{suchnost}.png",bbox_inches="tight")
    plt.savefig(f"all_my_systems_{suchnost}.svg",bbox_inches="tight")

    plt.clf()
    plt.cla()
    plt.close()
    print("[2] Making charts #2...")

    counts = 10000

    suchnost_erD1 = numpy.zeros((len(currenttable), counts))
    suchnost_erD2 = numpy.zeros((len(currenttable), counts))

    j = 0
    for index, row in currenttable.iterrows():

        suchnostExo1 = []
        suchnostExo2 = []
        if suchnost == "mass":
            upperLimitsuchnost1 = (row[f"{suchnost}_error_max-1"] + row[f"{suchnost}-1"])  * mass_coeff
            lowerLimitsuchnost1 = (row[f"{suchnost}-1"] - abs(row[f"{suchnost}_error_min-1"])) * mass_coeff
            upperLimitsuchnost2 = (row[f"{suchnost}_error_max-2"] + row[f"{suchnost}-2"])  * mass_coeff
            lowerLimitsuchnost2 = (row[f"{suchnost}-2"] - abs(row[f"{suchnost}_error_min-2"])) * mass_coeff
        elif suchnost == "radius":
            upperLimitsuchnost1 = (row[f"{suchnost}_error_max-1"] + row[f"{suchnost}-1"])  * radius_coeff
            lowerLimitsuchnost1 = (row[f"{suchnost}-1"] - abs(row[f"{suchnost}_error_min-1"])) * radius_coeff
            upperLimitsuchnost2 = (row[f"{suchnost}_error_max-2"] + row[f"{suchnost}-2"])  * radius_coeff
            lowerLimitsuchnost2 = (row[f"{suchnost}-2"] - abs(row[f"{suchnost}_error_min-2"])) * radius_coeff
        else:
            upperLimitsuchnost1 = (row[f"{suchnost}_error_max-1"] + row[f"{suchnost}-1"])
            lowerLimitsuchnost1 = (row[f"{suchnost}-1"] - abs(row[f"{suchnost}_error_min-1"]))
            upperLimitsuchnost2 = (row[f"{suchnost}_error_max-2"] + row[f"{suchnost}-2"])
            lowerLimitsuchnost2 = (row[f"{suchnost}-2"] - abs(row[f"{suchnost}_error_min-2"]))
        # print(upperLimitmass1, lowerLimitmass1)
        suchnostExo1 = numpy.random.uniform(upperLimitsuchnost1, lowerLimitsuchnost1, counts)
        suchnost_erD1[j] = suchnostExo1
        suchnostExo2 = numpy.random.uniform(upperLimitsuchnost2, lowerLimitsuchnost2, counts)
        suchnost_erD2[j] = suchnostExo2
        j +=1


    pearsonmassive = []
    pvaluemassive = []

    for i in range(counts):
        pearson_coef = 0.
        spearman_coef = 0.
        p_value = 0.
        s_value = 0.
        y1 = []
        x1 = []
        for k in range(len(currenttable)):
            y1.append(suchnost_erD2[k][i])
            x1.append(suchnost_erD1[k][i])

        pearson_coef, p_value = stats.pearsonr(x1, y1) #define the columns to perform calculations on

        pearsonmassive.append(pearson_coef)
        pvaluemassive.append(p_value)

    seaborn.set_theme(style="white")
    fig, ax = plt.subplots(figsize=(7,5))

    seaborn.cubehelix_palette(start=-.2,rot=.6, as_cmap=True)
    data = numpy.array(pearsonmassive)
    df = len(currenttable) - 2
    c=numpy.median(data)
    e = data.mean()
    t = e * numpy.sqrt(df/(1 - e**2))
    p = stats.t.sf(t, df)

    print(f"{suchnost} median={format(c,'.3f')}, mean={format(e,'.3f')}, p-value={'{:.2e}'.format(p)}")
    z1 = numpy.arctanh(e)

    sigma1 = (1/((len(currenttable)-3)**0.5))


    cint1 = z1 + numpy.array([-1, 1]) * sigma1 * stats.norm.ppf((1+0.95)/2)


    confidence_intervals1 =  numpy.tanh(cint1)
    print(f"Confidence intervals for Monte-Carlo sim: {confidence_intervals1}")

    h = seaborn.distplot(data, hist=True, norm_hist=True, kde=False, fit=stats.norm, fit_kws={"color": "b", "lw": 2},bins=100, hist_kws={"weights":None, "density":True, "color": "b"}, ax=ax, label=f"Pearson R")

    ax.axvline(x=c,color='k', linestyle="--", linewidth=2,label=f"Median={format(c,'.3f')}")
    ax.axvline(x=e,color='b', linestyle="-.", linewidth=2,label=f"Mean={format(e,'.3f')}")

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.xaxis.set_major_locator(plt.MaxNLocator(4))
    ax.set_xlabel(r'R(p$_{{i3/i2}}$/p$_{{i2/i1}}$)', fontsize=20)

    ax.grid(False)
    plt.savefig(f"./{suchnost}_{counts}_sim.png", bbox_inches="tight")
    plt.savefig(f"./{suchnost}_{counts}_sim.svg", bbox_inches="tight")

    tblStats = pandas.DataFrame()


    for stellarParam in listForstellarParam:
        x = 40
        n = 2
        PearsonOne = []
        PvalPe = []

        star_median = []
        lenght = range(0, len(currenttable) - x, n)
        mylist =  list(lenght)
        l = len(mylist)

        currenttable1 = currenttable.query(f"`{stellarParam}`.notna()")
        currenttable1 = currenttable1.sort_values(stellarParam)




        for i in range(0, len(currenttable1) - x, n):
            star_t = 0.
            pearson_coef = 0.
            p_value  = 0.

            param_1 = []
            param_2 = []
            newtable = currenttable1.iloc[i:i+x]
            param_1 = newtable[f"{suchnost}-1"]
            param_2 = newtable[f"{suchnost}-2"]

            pearson_coef, p_value = stats.pearsonr(param_1, param_2) #define the columns to perform calculations on

            star_t = newtable[f"{stellarParam}"].sum()/x
            PearsonOne.append(pearson_coef)
            PvalPe.append(p_value)

            star_median.append(star_t)
            if abs(pearson_coef) > 0.3:
                newtable.to_pickle(f"./data/{i}_{pearson_coef}_{stellarParam}_{suchnost}.pkl")


        c = len(star_median)
        for k in range(l-c):
            PearsonOne.append(numpy.nan)
            PvalPe.append(numpy.nan)
            star_median.append(numpy.nan)

        tblStats[f"{stellarParam}_x"] = star_median
        tblStats[f"{stellarParam}_Pearson"] = PearsonOne
        tblStats[f"{stellarParam}_Pvalue"] = PvalPe

    tblStats.to_pickle(f"./data/pvalues_{suchnost}.pkl")
    print("[3] Making charts #3...")

    seaborn.set_theme(style="white")
    fig, ax = plt.subplots(1,5,
                  figsize=(15,3),
                  sharey=True)
    plt.grid(False)
    minaxx = tblStats[["star_mass_Pvalue", "star_radius_Pvalue", "star_metallicity_Pvalue", "star_teff_Pvalue", "star_age_Pvalue"]].min(axis=1)
    maxaxx = tblStats[["star_mass_Pvalue", "star_radius_Pvalue", "star_metallicity_Pvalue", "star_teff_Pvalue", "star_age_Pvalue"]].max(axis=1)

    minmass = tblStats["star_mass_x"]
    maxmass = tblStats["star_mass_x"]

    seaborn.lineplot(data = tblStats, x="star_mass_x",
            y="star_mass_Pvalue",
            ax=ax[0],
    )
    minradius = tblStats["star_radius_x"]
    maxradius = tblStats["star_radius_x"]
    seaborn.lineplot(data = tblStats, x="star_radius_x",
            y="star_radius_Pvalue",
            ax=ax[1]
    )
    minmetallicity = tblStats["star_metallicity_x"]
    maxmetallicity = tblStats["star_metallicity_x"]
    seaborn.lineplot(data = tblStats, x="star_metallicity_x",
            y="star_metallicity_Pvalue",
            ax=ax[2]
    )
    minteff = tblStats["star_teff_x"]
    maxteff = tblStats["star_teff_x"]
    seaborn.lineplot(data = tblStats, x="star_teff_x",
            y="star_teff_Pvalue",
            ax=ax[3]
    )
    minage = tblStats["star_age_x"]
    maxage = tblStats["star_age_x"]
    seaborn.lineplot(data = tblStats, x="star_age_x",
            y="star_age_Pvalue",
            ax=ax[4]
    )
    ax[0].set_ylim((min(minaxx) * 0.8, max(maxaxx) * 1.2))
    ax[0].set_xlim((min(minmass) * 0.8, max(maxmass) * 1.2))

    ax[1].set_ylim((min(minaxx) * 0.8, max(maxaxx)* 1.2))
    ax[1].set_xlim((min(minradius) * 0.8, max(maxradius) * 1.2))

    ax[2].set_ylim((min(minaxx) * 0.8, max(maxaxx)* 1.2))
    ax[2].set_xlim((min(minmetallicity) * 0.8, max(maxmetallicity) * 1.2))

    ax[3].set_ylim((min(minaxx) * 0.8, max(maxaxx)* 1.2))
    ax[3].set_xlim((min(minteff) * 0.8, max(maxteff) * 1.2))

    ax[4].set_ylim((min(minaxx) * 0.8, max(maxaxx)* 1.2))
    ax[4].set_xlim((min(minage) * 0.8, max(maxage) * 1.2))

    ax[0].set_yscale("log")
    ax[1].set_yscale("log")
    ax[2].set_yscale("log")
    ax[3].set_yscale("log")
    ax[4].set_yscale("log")
    ax[0].set_xlabel(f"Mass range [M$_\odot$]", fontsize=16)
    ax[1].set_xlabel(f"Radius range [R$_\odot$]", fontsize=16)
    ax[2].set_xlabel("Metallicity range [Fe/H]", fontsize=16)
    ax[3].set_xlabel("Temperature range [T]", fontsize=16)
    ax[4].set_xlabel("Age range [Gyr]", fontsize=16)
    ax[0].tick_params(axis='both', which='major', labelsize=16)
    ax[1].tick_params(axis='both', which='major', labelsize=16)
    ax[2].tick_params(axis='both', which='major', labelsize=16)
    ax[3].tick_params(axis='both', which='major', labelsize=16)
    ax[4].tick_params(axis='both', which='major', labelsize=16)
    ax[0].set_ylabel("P value", fontsize=16)
    plt.grid(False)
    fig.tight_layout()
    plt.savefig(f"PvalueMovingSample_{suchnost}.png",bbox_inches="tight")
    plt.savefig(f"PvalueMovingSample_{suchnost}.svg",bbox_inches="tight")

    fig, ax = plt.subplots(1,5,
                  figsize=(15,3),
                  sharey=True)
    plt.grid(False)
    seaborn.lineplot(data = tblStats, x="star_mass_x",
            y="star_mass_Pearson",
            ax=ax[0],
    )
    seaborn.lineplot(data = tblStats, x="star_radius_x",
            y="star_radius_Pearson",
            ax=ax[1]
    )
    seaborn.lineplot(data = tblStats, x="star_metallicity_x",
            y="star_metallicity_Pearson",
            ax=ax[2]
    )
    seaborn.lineplot(data = tblStats, x="star_teff_x",
            y="star_teff_Pearson",
            ax=ax[3]
    )
    seaborn.lineplot(data = tblStats, x="star_age_x",
            y="star_age_Pearson",
            ax=ax[4]
    )
    ax[0].set_xlim((min(minmass) * 0.8, max(maxmass) * 1.2))
    ax[1].set_xlim((min(minradius) * 0.8, max(maxradius) * 1.2))
    ax[2].set_xlim((min(minmetallicity) * 0.8, max(maxmetallicity) * 1.2))
    ax[3].set_xlim((min(minteff) * 0.8, max(maxteff) * 1.2))
    ax[4].set_xlim((min(minage) * 0.8, max(maxage) * 1.2))

    ax[0].set_xlabel(f"Mass range [M$_\odot$]", fontsize=16)
    ax[1].set_xlabel(f"Radius range [R$_\odot$]", fontsize=16)
    ax[2].set_xlabel("Metallicity range [Fe/H]", fontsize=16)
    ax[3].set_xlabel("Temperature range [T]", fontsize=16)
    ax[4].set_xlabel("Age range [Gyr]", fontsize=16)
    ax[0].set_ylabel("Pearson R value", fontsize=16)
    ax[0].tick_params(axis='both', which='major', labelsize=16)
    ax[1].tick_params(axis='both', which='major', labelsize=16)
    ax[2].tick_params(axis='both', which='major', labelsize=16)
    ax[3].tick_params(axis='both', which='major', labelsize=16)
    ax[4].tick_params(axis='both', which='major', labelsize=16)
    ax[0].grid(False)
    ax[1].grid(False)
    ax[2].grid(False)
    ax[3].grid(False)
    ax[4].grid(False)
    fig.tight_layout()
    plt.savefig(f"PearsonMovingSample_{suchnost}.png",bbox_inches="tight")
    plt.savefig(f"PearsonMovingSample_{suchnost}.svg",bbox_inches="tight")



