import pyvo
import numpy
from tabulate import tabulate
import pandas
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
import matplotlib
warnings.filterwarnings("ignore", category=UserWarning)
pandas.options.mode.use_inf_as_na = True
import math
service = pyvo.dal.TAPService("http://voparis-tap-planeto.obspm.fr/tap")
import seaborn
from scipy.stats import norm
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
import matplotlib.ticker


quantity_support()
set_matplotlib_formats('svg')
plt.rc('legend', frameon=False)
plt.rc('figure', figsize=(10, 10))
plt.rc('font', size=12)

workingTableExoplanets = pandas.read_pickle("./data/MRP_data_sets.pkl")

workingTableExoplanets = workingTableExoplanets.query("`ima_flag`.isna()")

workingTableExoplanets = workingTableExoplanets[~workingTableExoplanets['star_name'].isin(['KOI-55'])] # marked as controversial in NASA

workingTableExoplanets['density-1'] = (workingTableExoplanets["mass-1"] * mass_coeff * constants.M_earth) * 3 /(4 * numpy.pi*numpy.power(workingTableExoplanets["radius-1"] * radius_coeff * constants.R_earth,3))
workingTableExoplanets['density-2'] = (workingTableExoplanets["mass-2"] * mass_coeff * constants.M_earth) * 3 /(4 * numpy.pi*numpy.power(workingTableExoplanets["radius-2"] * radius_coeff * constants.R_earth,3))

workingTableExoplanets['density_error_min-1'] = numpy.sqrt((workingTableExoplanets["mass_error_min-1"]/workingTableExoplanets["mass-1"]) ** 2 + (3 * workingTableExoplanets["radius_error_min-1"]/workingTableExoplanets["radius-1"]) ** 2) * workingTableExoplanets['density-1']# assigned to a column
workingTableExoplanets['density_error_max-1'] = numpy.sqrt((workingTableExoplanets["mass_error_max-1"]/workingTableExoplanets["mass-1"]) ** 2 + (3 * workingTableExoplanets["radius_error_max-1"]/workingTableExoplanets["radius-1"]) ** 2) * workingTableExoplanets['density-1']
workingTableExoplanets['density_error_min-2'] = numpy.sqrt((workingTableExoplanets["mass_error_min-2"]/workingTableExoplanets["mass-2"]) ** 2 + (3 * workingTableExoplanets["radius_error_min-2"]/workingTableExoplanets["radius-2"]) ** 2) * workingTableExoplanets['density-2']# assigned to a column
workingTableExoplanets['density_error_max-2'] = numpy.sqrt((workingTableExoplanets["mass_error_max-2"]/workingTableExoplanets["mass-2"]) ** 2 + (3 * workingTableExoplanets["radius_error_max-2"]/workingTableExoplanets["radius-2"]) ** 2) * workingTableExoplanets['density-2']

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
    "radius",
    "mass",
    "density"
]
seaborn.set_theme(style="white")
for suchnost in trisuchnosti:
    if suchnost != "radius":
        currenttable = workingTableExoplanets[~workingTableExoplanets['mass_detection_type'].isin(['Theoretical'])]
    else:
        currenttable = workingTableExoplanets

    print("[1] Making charts #1...")
    currenttable = currenttable.query(f"`{suchnost}-1`.notna() & `{suchnost}-2`.notna()")
    currenttable = currenttable.query(f"`{suchnost}_error_min-1`.notna() & `{suchnost}_error_max-1`.notna() & `{suchnost}_error_min-2`.notna() & `{suchnost}_error_max-2`")


    pearson_coef= []
    p_value  = []
    spearman_coef = []
    s_value = []

    planetstodrop = []

    for index, row in currenttable.iterrows():
        if (abs(currenttable.at[index, f"{suchnost}_error_min-1"]) + currenttable.at[index, f"{suchnost}_error_max-1"])/currenttable.at[index, f"{suchnost}-1"] >= 2.0:
            planetstodrop.append(index)
        if (abs(currenttable.at[index, f"{suchnost}_error_min-2"]) + currenttable.at[index, f"{suchnost}_error_max-2"])/currenttable.at[index, f"{suchnost}-2"] >= 2.0:
            planetstodrop.append(index)
    print(len(planetstodrop))

    currenttable = currenttable.drop(index=planetstodrop)

    resultsExoplanetsm2 = currenttable.query(f"`{suchnost}-2` > `{suchnost}-1`")
    resultsExoplanetsm1 = currenttable.query(f"`{suchnost}-2` < `{suchnost}-1`")
    resultsExoplanetsmeq = currenttable.query(f"`{suchnost}-2` == `{suchnost}-1`")
    errror_frac = (abs(currenttable[f"{suchnost}_error_min-1"]) +currenttable[f"{suchnost}_error_max-1"])/currenttable[f"{suchnost}-1"] + (abs(currenttable[f"{suchnost}_error_min-2"]) +currenttable[f"{suchnost}_error_max-2"])/currenttable[f"{suchnost}-2"]

    if suchnost == "mass":

        suchnost_1 = currenttable[f"{suchnost}-1"] * mass_coeff
        suchnost_2 = currenttable[f"{suchnost}-2"] * mass_coeff
        # errror_frac = (abs(currenttable[f"{suchnost}_error_min-1"]) +currenttable[f"{suchnost}_error_max-1"]) * mass_coeff/suchnost_1 + (abs(currenttable[f"{suchnost}_error_min-2"]) +currenttable[f"{suchnost}_error_max-2"]) * mass_coeff/suchnost_2
        xmin= currenttable[[f"{suchnost}-1", f"{suchnost}-2"]].min(axis=1).min(axis=0) * mass_coeff * 0.4
        ymax= currenttable[[f"{suchnost}-1", f"{suchnost}-2"]].max(axis=1).max(axis=0) * mass_coeff * 1.6
        xerr = numpy.nan_to_num(currenttable[[f"{suchnost}_error_min-1", f"{suchnost}_error_max-1"]].to_numpy().T, posinf=0.) * mass_coeff
        yerr = numpy.nan_to_num(currenttable[[f"{suchnost}_error_min-2", f"{suchnost}_error_max-2"]].to_numpy().T, posinf=0.) * mass_coeff
        size = 0.34280936454849503
    elif suchnost == "radius":
        suchnost_1 = currenttable[f"{suchnost}-1"] * radius_coeff
        suchnost_2 = currenttable[f"{suchnost}-2"] * radius_coeff
        # errror_frac = (abs(currenttable[f"{suchnost}_error_min-1"]) +currenttable[f"{suchnost}_error_max-1"]) * radius_coeff/suchnost_1 + (abs(workingTableExoplanets[f"{suchnost}_error_min-2"]) +workingTableExoplanets[f"{suchnost}_error_max-2"]) * radius_coeff/suchnost_2
        xmin= currenttable[[f"{suchnost}-1", f"{suchnost}-2"]].min(axis=1).min(axis=0) * radius_coeff * 0.5
        ymax= currenttable[[f"{suchnost}-1", f"{suchnost}-2"]].max(axis=1).max(axis=0) * radius_coeff * 1.5
        xerr = numpy.nan_to_num(currenttable[[f"{suchnost}_error_min-1", f"{suchnost}_error_max-1"]].to_numpy().T, posinf=0.) * radius_coeff
        yerr = numpy.nan_to_num(currenttable[[f"{suchnost}_error_min-2", f"{suchnost}_error_max-2"]].to_numpy().T, posinf=0.) * radius_coeff
        size = 0.2665551839464883
    else:
        suchnost_1 = currenttable[f"{suchnost}-1"]
        suchnost_2 = currenttable[f"{suchnost}-2"]
        xmin= currenttable[[f"{suchnost}-1", f"{suchnost}-2"]].min(axis=1).min(axis=0) * 0.5
        ymax= currenttable[[f"{suchnost}-1", f"{suchnost}-2"]].max(axis=1).max(axis=0) * 1.5
        xerr = numpy.nan_to_num(currenttable[[f"{suchnost}_error_min-1", f"{suchnost}_error_max-1"]].to_numpy().T, posinf=0.)
        yerr = numpy.nan_to_num(currenttable[[f"{suchnost}_error_min-2", f"{suchnost}_error_max-2"]].to_numpy().T, posinf=0.)
        size = 0.
    print(f"Working with {suchnost}, {len(currenttable)} pairs")
    currenttable.to_pickle(f"./data/{suchnost}_only.pkl")
    # checking for ordering in pairs
    print(f'Upper side of the line, planets:', len(resultsExoplanetsm2), ', part:', numpy.round(len(resultsExoplanetsm2)/ len(currenttable), 3))
    print(f'Lower side of the line, planets:', len(resultsExoplanetsm1), ', part:', numpy.round(len(resultsExoplanetsm1)/ len(currenttable), 3))
    print(f'On the line, planets:', len(resultsExoplanetsmeq))

    pearson_coef, p_value = stats.pearsonr(suchnost_1, suchnost_2) #define the columns to perform calculations on
    print(f"Pearson Correlation Coefficient for {suchnost}: {pearson_coef}, and a P-value of: {p_value}") # Results
    # spearman_coef, s_value = stats.spearmanr(suchnost_1, suchnost_2) #define the columns to perform calculations on
    # print("Spearman Correlation Coefficient for raw data: ", spearman_coef, "and a P-value of:", s_value) # Resultsc
    res = stats.pearsonr(suchnost_1, suchnost_2)
    print(res.confidence_interval())
    z = numpy.arctanh(pearson_coef)

    sigma = (1/((len(currenttable)-3)**0.5))


    cint = z + numpy.array([-1, 1]) * sigma * stats.norm.ppf((1+0.95)/2)


    confidence_intervals =  numpy.tanh(cint)
    print(f"Pearson Correlation Coefficient confidence intervals: {confidence_intervals}")

    cmap = seaborn.cubehelix_palette(rot=-.3, as_cmap=True)
    # fig, ax = plt.subplots(figsize=(6,6))
    seaborn.set_theme(style="white")
    fig, ax = plt.subplots(
                        figsize=(8.5, 8.5),
                        )


    plt.rcParams['figure.figsize']=(6,6)
    seaborn.scatterplot(
        x=suchnost_1, y=suchnost_2,
        palette=cmap, legend=False, ax=ax, zorder=2
    )
    ax.set(xscale="log", yscale="log")
    ax.grid(True)
    ax.set(ylim=(xmin, ymax))
    ax.set(xlim=(xmin, ymax))
    # g._legend.remove()
    x = numpy.linspace(xmin*0.9, ymax*1.1, 2000)
    y = x
    ax.plot(x, y, linewidth=0.8, linestyle='--', color='k',zorder=0)
    ax.errorbar(suchnost_1, suchnost_2, xerr=numpy.abs(xerr), yerr=numpy.abs(yerr), ls='none', fmt='0.8', ecolor='tab:gray', elinewidth=0.8, capsize=None, barsabove=True, zorder=1)
    ax.scatter(suchnost_1, suchnost_2, marker="o", facecolor='tab:blue', zorder=0, label="Data set")
    solarsystemTable = pandas.read_pickle("./data/solarsystemE.pkl")

    if suchnost != "mass_radius_ratio":
        for i in range(len(solarsystemTable[f"{suchnost}E"])-1):
            if i == 1:
                ax.plot(solarsystemTable[f"{suchnost}E"][i], solarsystemTable[f"{suchnost}E"][i+1], 'ro', ms=6, zorder=3, label="Solar system")
            else:
                ax.plot(solarsystemTable[f"{suchnost}E"][i], solarsystemTable[f"{suchnost}E"][i+1], 'ro', ms=6, zorder=3)

    f = matplotlib.ticker.ScalarFormatter(useMathText=True)
    f.set_powerlimits((-3,3))
    ax5 = fig.add_axes([0.7, 0.05, 0.25, 0.25])#, sharex=ax)
    dist = seaborn.ecdfplot(x=numpy.log10(suchnost_2/suchnost_1),stat='count', ax=ax5, label=f"{suchnost}")
    ax5.axhline(y=len(numpy.log10(suchnost_2/suchnost_1))/2.,color='r', linestyle="--",linewidth=0.8, )
    ax5.axvline(x=0.,color='k', linestyle="--",linewidth=0.8, )
    ax5.grid(False)
    mean = numpy.median(suchnost_2/suchnost_1)
    if suchnost == "density":
        ax.set_ylabel(f"D$_{{p(i+1)}}$ [kg m$^{{-3}}$]", fontsize=14)
        ax.set_xlabel(f"D$_{{p(i)}}$ [kg m$^{{-3}}$]", fontsize=14)
        ax5.set_xlabel(f'log$_{{10}}$(D$_{{i+1}}$/D$_{{i}}$)',fontsize=14)

    elif suchnost == "radius":
        ax.set_xlabel(f"R$_{{p(i)}}$ [R$_\\oplus$]", fontsize=14)
        ax.set_ylabel(f"R$_{{p(i+1)}}$  [R$_\\oplus$] ", fontsize=14)
        ax5.set_xlabel(f'log$_{{10}}$(R$_{{i+1}}$/R$_{{i}}$)', fontsize=14)


    elif suchnost == "mass":
        ax.set_xlabel(f"M$_{{p(i)}}$ [M$_\\oplus$]", fontsize=14)
        ax.set_ylabel(f"M$_{{p(i+1)}}$  [M$_\\oplus$] ", fontsize=14)
        ax.set_xlim(8e-2,1e3)
        ax.set_ylim(8e-2,1e3)
        ax5.set_xlabel(f'log$_{{10}}$(M$_{{i+1}}$/M$_{{i}}$)', fontsize=14)

    ax5.set_ylabel(f'CDF', fontsize=14)
    ax5.yaxis.tick_right()
    ax5.yaxis.set_label_position('right')
    plt.grid(False)
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
        suchnostExo1 = numpy.random.uniform(upperLimitsuchnost1, lowerLimitsuchnost1, counts)
        suchnost_erD1[j] = suchnostExo1
        suchnostExo2 = numpy.random.uniform(upperLimitsuchnost2, lowerLimitsuchnost2, counts)
        suchnost_erD2[j] = suchnostExo2
        j +=1



    pearsonmassive = []
    pvaluemassive = []

    spearmanmassive = []
    svaluemassive = []
    mean_radius_massive = []
    mean_mass_massive = []
    median_radius_massive = []
    median_mass_massive = []
    for i in range(counts):
        pearson_coef = 0.
        spearman_coef = 0.
        p_value = 0.
        s_value = 0.
        mean_radius = 0.
        mean_mass = 0.
        median_radius = 0.
        median_mass = 0.
        y1 = []
        x1 = []
        for k in range(len(currenttable)):
            y1.append(suchnost_erD2[k][i])
            x1.append(suchnost_erD1[k][i])
            # p_zero1 = [1.03,0.29] p0=p_zero1,
        mean_radius = numpy.mean(y1)
        mean_mass = numpy.mean(x1)
        median_radius = numpy.median(y1)
        median_mass = numpy.median(x1)

        pearson_coef, p_value = stats.pearsonr(x1, y1) #define the columns to perform calculations on
        # print("Pearson Correlation Coefficient for raw data: ", pearson_coef, "and a P-value of:", p_value) # Results

        pearsonmassive.append(pearson_coef)
        pvaluemassive.append(p_value)

        mean_radius_massive.append(mean_radius)
        mean_mass_massive.append(mean_mass)
        median_radius_massive.append(median_radius)
        median_mass_massive.append(median_mass)
    seaborn.set_theme(style="white")
    fig, ax = plt.subplots(1,2,figsize=(14,5))

    seaborn.cubehelix_palette(start=-.2,rot=.6, as_cmap=True)
    data = numpy.array(pearsonmassive)

    seaborn.set(rc={'figure.figsize':(14,5)})


    h = seaborn.distplot(data, hist=True, norm_hist=True, kde=False, fit=stats.norm, fit_kws={"color": "b", "lw": 2},bins=100, hist_kws={"weights":None, "density":True, "color": "b"}, ax=ax[0], label=f"R-value")

    data1 = numpy.array(pvaluemassive)
    h1 = seaborn.distplot(numpy.log10(data1), hist=True, norm_hist=True, kde=False,
        fit=stats.norm, fit_kws={"color": "b", "lw": 2},bins=100, hist_kws={"weights":None, "density":True, "color": "b"}, ax=ax[1], label=f"P-value")
    e = data.mean()
    d = numpy.log10(data1.mean())
    c=numpy.median(data)
    b=numpy.log10(numpy.median(data1))
    f = matplotlib.ticker.ScalarFormatter(useMathText=True)
    f.set_powerlimits((-3,3))

    ax[0].axvline(x=c,color='g', label=f"Median={format(c,'.3f')}")
    ax[0].axvline(x=e,color='r', label=f"Mean={format(e,'.3f')}")

    frac1, whole1 = math.modf(b)
    b2 = whole1 - 1
    b1 = 10 ** b * 10 ** abs(b2)
    ax[1].axvline(x=b,color='g', label=f"Median={b1.round(1)}x10$^{{{round(b2)}}}$")
    frac2, whole2 = math.modf(d)
    d2 = whole2 - 1
    d1 = 10 ** d * 10 ** abs(d2)
    ax[1].axvline(x=d,color='r', label=f"Mean={d1.round(1)}x10$^{{{round(d2)}}}$")

    ax[0].legend(loc="upper right", prop={"size":13.5})
    ax[1].legend(loc="upper left", prop={"size":13.5})
    if suchnost == "radius":
        ax[0].legend(loc="upper left", prop={"size":13.5})
        ax[1].legend(loc="upper left", prop={"size":13.5})
    ax[0].tick_params(axis='both', which='major', labelsize=16)

    ax[1].tick_params(axis='both', which='major', labelsize=16)
    # ax[0,0].legend()
    # ax[0,1].legend()
    ax[0].grid(False)
    ax[1].grid(False)

    plt.savefig(f"./{suchnost}_{counts}_sim.png", bbox_inches="tight")
    plt.savefig(f"./{suchnost}_{counts}_sim.svg", bbox_inches="tight")

    tblStats = pandas.DataFrame()

    for stellarParam in listForstellarParam:
        x = 40
        n = 1
        PearsonOne = []
        PvalPe = []
        SpearmanOne = []
        PvalSp = []
        star_median = []
        lenght = range(0, len(currenttable) - x, n)
        mylist =  list(lenght)
        l = len(mylist)

        currenttable1 = currenttable.query(f"`{stellarParam}`.notna()")

        currenttable1 = currenttable1.sort_values(stellarParam)




        for i in range(0, len(currenttable1) - x, n): #len(workingTableExoplanets) - x):
            star_t = 0.
            pearson_coef = 0.
            p_value  = 0.
            spearman_coef = 0.
            s_value = 0.
            param_1 = []
            param_2 = []
            newtable = currenttable1.iloc[i:i+x]
            param_1 = newtable[f"{suchnost}-1"]
            param_2 = newtable[f"{suchnost}-2"]

            pearson_coef, p_value = stats.pearsonr(param_1, param_2) #define the columns to perform calculations on
            spearman_coef, s_value = stats.spearmanr(param_1, param_2) #define the columns to perform calculations on

            star_t = newtable[f"{stellarParam}"].sum()/x
            PearsonOne.append(pearson_coef)
            PvalPe.append(p_value)
            SpearmanOne.append(spearman_coef)
            PvalSp.append(s_value)
            star_median.append(star_t)
            if abs(pearson_coef) > 0.84:
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

tblStats1 = pandas.read_pickle(f"./data/pvalues_density.pkl")
tblStats2 = pandas.read_pickle(f"./data/pvalues_radius.pkl")
tblStats0 = pandas.read_pickle(f"./data/pvalues_mass.pkl")

print("[3] Making charts #3...")
seaborn.set_theme(style="white")
fig, ax = plt.subplots(1,5,
                  figsize=(20,4.2),#15,3
                  sharey=True)


seaborn.lineplot(data = tblStats0, x="star_mass_x",
        y="star_mass_Pvalue",
        ax=ax[0], color="tab:green",

)
seaborn.lineplot(data = tblStats1, x="star_mass_x",
        y="star_mass_Pvalue",
        ax=ax[0], color="k",linestyle="--",

)
seaborn.lineplot(data = tblStats2, x="star_mass_x",
        y="star_mass_Pvalue",
        ax=ax[0], color="tab:blue",

)

seaborn.lineplot(data = tblStats0, x="star_radius_x",
        y="star_radius_Pvalue",
        ax=ax[1], color="tab:green",

)
seaborn.lineplot(data = tblStats1, x="star_radius_x",
        y="star_radius_Pvalue",
        ax=ax[1], color="k",linestyle="--",

)
seaborn.lineplot(data = tblStats2, x="star_radius_x",
        y="star_radius_Pvalue",
        ax=ax[1], color="tab:blue",

)

seaborn.lineplot(data = tblStats0, x="star_metallicity_x",
        y="star_metallicity_Pvalue",
        ax=ax[2], color="tab:green",

)
seaborn.lineplot(data = tblStats1, x="star_metallicity_x",
        y="star_metallicity_Pvalue",
        ax=ax[2], color="k",linestyle="--",

)
seaborn.lineplot(data = tblStats2, x="star_metallicity_x",
        y="star_metallicity_Pvalue",
        ax=ax[2], color="tab:blue",

)


seaborn.lineplot(data = tblStats0, x="star_teff_x",
        y="star_teff_Pvalue",
        ax=ax[3], color="tab:green",

)
seaborn.lineplot(data = tblStats1, x="star_teff_x",
        y="star_teff_Pvalue",
        ax=ax[3], color="k", linestyle="--",

)
seaborn.lineplot(data = tblStats2, x="star_teff_x",
        y="star_teff_Pvalue",
        ax=ax[3], color="tab:blue",

)

seaborn.lineplot(data = tblStats0, x="star_age_x",
        y="star_age_Pvalue",
        ax=ax[4], color="tab:green",

)
seaborn.lineplot(data = tblStats1, x="star_age_x",
        y="star_age_Pvalue",
        ax=ax[4], color="k", linestyle="--",

)
seaborn.lineplot(data = tblStats2, x="star_age_x",
        y="star_age_Pvalue",
        ax=ax[4], color="tab:blue",

)


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
fig.tight_layout()
ax[0].grid(False)
ax[1].grid(False)
ax[2].grid(False)
ax[3].grid(False)
ax[4].grid(False)

plt.savefig(f"PvalueMovingSample_all.png",bbox_inches="tight")
plt.savefig(f"PvalueMovingSample_all.svg",bbox_inches="tight")
seaborn.set_theme(style="white")

fig, ax = plt.subplots(1,5,
              figsize=(20,4.2),#15,3
              sharey=True)
seaborn.lineplot(data = tblStats0, x="star_mass_x",
        y="star_mass_Pearson",
        ax=ax[0], color="tab:green", zorder=2
    #     estimator=None,
    # linewidth=0.8
)
seaborn.lineplot(data = tblStats1, x="star_mass_x",
        y="star_mass_Pearson",
        ax=ax[0], color="tab:gray",linestyle="--", zorder=3

)
seaborn.lineplot(data = tblStats2, x="star_mass_x",
        y="star_mass_Pearson",
        ax=ax[0], color="tab:blue", zorder=1

)

seaborn.lineplot(data = tblStats0, x="star_radius_x",
        y="star_radius_Pearson",
        ax=ax[1], color="tab:green",

)
seaborn.lineplot(data = tblStats1, x="star_radius_x",
        y="star_radius_Pearson",
        ax=ax[1], color="k",linestyle="--",

)
seaborn.lineplot(data = tblStats2, x="star_radius_x",
        y="star_radius_Pearson",
        ax=ax[1], color="tab:blue",

)

seaborn.lineplot(data = tblStats0, x="star_metallicity_x",
        y="star_metallicity_Pearson",
        ax=ax[2], color="tab:green",

)
seaborn.lineplot(data = tblStats1, x="star_metallicity_x",
        y="star_metallicity_Pearson",
        ax=ax[2], color="k",linestyle="--",

)
seaborn.lineplot(data = tblStats2, x="star_metallicity_x",
        y="star_metallicity_Pearson",
        ax=ax[2], color="tab:blue",

)

seaborn.lineplot(data = tblStats0, x="star_teff_x",
        y="star_teff_Pearson",
        ax=ax[3], color="tab:green",

)
seaborn.lineplot(data = tblStats1, x="star_teff_x",
        y="star_teff_Pearson",
        ax=ax[3], color="k",linestyle="--",

)
seaborn.lineplot(data = tblStats2, x="star_teff_x",
        y="star_teff_Pearson",
        ax=ax[3], color="tab:blue",

)

seaborn.lineplot(data = tblStats0, x="star_age_x",
        y="star_age_Pearson",
        ax=ax[4], color="tab:green",

)
seaborn.lineplot(data = tblStats1, x="star_age_x",
        y="star_age_Pearson",
        ax=ax[4], color="k",linestyle="--",

)
seaborn.lineplot(data = tblStats2, x="star_age_x",
        y="star_age_Pearson",
        ax=ax[4], color="tab:blue",

)



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
plt.savefig(f"PearsonMovingSample_all.png",bbox_inches="tight")
plt.savefig(f"PearsonMovingSample_all.svg",bbox_inches="tight")



