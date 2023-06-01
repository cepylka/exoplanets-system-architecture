import sys
import numpy
import pandas
import pyvo
import warnings
from tabulate import tabulate
from astropy import units
from astropy import constants
from astropy.visualization import quantity_support
from collections import OrderedDict
import itertools
from functools import partial
from cycler import cycler
from astropy import units
from astropy import constants
import pathlib

from uio.utility.files import pickle
from uio.utility.databases import tap

if len(sys.argv) != 2:
    raise SystemExit("[ERROR] You need to provide a path to a pickle")

pkl = pathlib.Path(sys.argv[1])
workingTable = pickle.openPickleAsPandasTable(pkl)

print(f"Processing pickle: {pkl.name}")

serviceNameNASA = "NASA"

dictForCheckingPlanet = {
    "semi_major_axis": "pl_orbsmax",
    "semi_major_axis_error_min": "pl_orbsmaxerr2",
    "semi_major_axis_error_max": "pl_orbsmaxerr1",
    "pl_orbsmaxlim": "pl_orbsmaxlim",
    "period": "pl_orbper",
    "period_error_min": "pl_orbpererr2",
    "period_error_max": "pl_orbpererr1",
    "pl_orbperlim": "pl_orbperlim",
    "eccentricity": "pl_orbeccen",
    "inclination": "pl_orbincl"
}

dictForCheckingStar = {
    "star_teff": "st_teff",
    "star_radius": "st_rad",
    "star_mass": "st_mass",
    "star_age": "st_age",
    "star_metallicity": "st_met",
    "st_metratio": "st_metratio",
    "st_lum": "st_lum",
    "sy_snum": "sy_snum",
    "st_rotp": "st_rotp",
    "sy_pnum": "sy_pnum",
    "cb_flag": "cb_flag",
    "ra": "ra",
    "sy_dist": "sy_dist"
}

dictForChecking = {}
dictForChecking.update(dictForCheckingPlanet.items())
dictForChecking.update(dictForCheckingStar.items())

for valueToAdd in dictForChecking.keys():
    workingTable[valueToAdd] = (
        numpy.array(numpy.NaN, dtype=float)
        if valueToAdd not in tap.services[serviceNameNASA][
            "parameters-that-are-strings"
        ]
        else numpy.array(numpy.NaN, dtype=str)
    )

print(workingTable)

print("\n[1/2] Getting parameters from NASA...")

# print(workingTable)

for index, row in workingTable.iterrows():
    systemName = row["star_name"]
    print(f"Iterating {index}, star name: {systemName}")
    for valueToAdd in dictForChecking.keys():
        if (
            dictForChecking[valueToAdd].endswith("err1")
            or dictForChecking[valueToAdd].endswith("err2")
        ):
            continue
        else:
            valueFromNASA = tap.getParameterFromNASA(
                systemName,
                index.replace("'", ""),
                dictForChecking[valueToAdd]
            )
            if valueFromNASA:
                workingTable.at[index, valueToAdd] = valueFromNASA
                # if it's a parameter with min/max errors, get them
                if valueToAdd in tap.services[serviceNameNASA][
                    "parameters-that-have-errors"
                ]:
                    print(f"---\nTrying to get min/max errors of {valueToAdd}")
                    valueFromNASAerr1, valueFromNASAerr2 = tap.getParameterErrorsFromNASA(
                        systemName,
                        index.replace("'", ""),
                        f"{dictForChecking[valueToAdd]}"
                    )
                    if valueFromNASAerr1:
                        workingTable.at[index, f"{valueToAdd}_error_max"] = valueFromNASAerr1
                    if valueFromNASAerr2:
                        workingTable.at[index, f"{valueToAdd}_error_min"] = valueFromNASAerr2
                    print("---")

print(workingTable)

print("\n[2/2] Habitable zone calculation...")

# first we create 6 rows for HZ coefficients
workingTable["st_recentVenus"] = numpy.array(numpy.NaN, dtype=float)
workingTable["st_runawayGreenhouse"] = numpy.array(numpy.NaN, dtype=float)
workingTable["st_maxGreenhouse"] = numpy.array(numpy.NaN, dtype=float)
workingTable["st_earlyMars"] = numpy.array(numpy.NaN, dtype=float)
workingTable["st_half_Earth"] = numpy.array(numpy.NaN, dtype=float)
workingTable["st_five_Earth"] = numpy.array(numpy.NaN, dtype=float)

# coeffcients to be used in the analytical expression to calculate
# habitable zone flux boundaries

seffsun = [1.776, 1.107, 0.356, 0.320, 1.188, 0.99]
a = [2.136e-4, 1.332e-4, 6.171e-5, 5.547e-5, 1.433e-4, 1.209e-4]
b = [2.533e-8, 1.580e-8, 1.698e-9, 1.526e-9, 1.707e-8, 1.404e-8]
c = [-1.332e-11, -8.308e-12, -3.198e-12, -2.874e-12, -8.968e-12, -7.418e-12]
d = [-3.097e-15, -1.931e-15, -5.575e-16, -5.011e-16, -2.084e-15, -1.713e-15]

for index, row in workingTable.iterrows():
    workingTable.at[index, "st_recentVenus"] = seffsun[0] + a[0] * (workingTable.at[index, "star_teff"] - 5780.0) + b[0] * (workingTable.at[index,"star_teff"] - 5780.0) ** 2 + c[0] * (workingTable.at[index,"star_teff"] - 5780.0) ** 3 + d[0] * (workingTable.at[index,"star_teff"] - 5780.0) ** 4
    workingTable.at[index, "st_runawayGreenhouse"] = seffsun[1] + a[1] * (workingTable.at[index, "star_teff"] - 5780.0) + b[1] * (workingTable.at[index,"star_teff"] - 5780.0) ** 2 + c[1] * (workingTable.at[index,"star_teff"] - 5780.0) ** 3 + d[1] * (workingTable.at[index,"star_teff"] - 5780.0) ** 4
    workingTable.at[index, "st_maxGreenhouse"] = seffsun[2] + a[2] * (workingTable.at[index, "star_teff"] - 5780.0) + b[2] * (workingTable.at[index,"star_teff"] - 5780.0) ** 2 + c[2] * (workingTable.at[index,"star_teff"] - 5780.0) ** 3 + d[2] * (workingTable.at[index,"star_teff"] - 5780.0) ** 4
    workingTable.at[index, "st_earlyMars"] = seffsun[3] + a[3] * (workingTable.at[index, "star_teff"] - 5780.0) + b[3] * (workingTable.at[index,"star_teff"] - 5780.0) ** 2 + c[3] * (workingTable.at[index,"star_teff"] - 5780.0) ** 3 + d[3] * (workingTable.at[index,"star_teff"] - 5780.0) ** 4
    workingTable.at[index, "st_half_Earth"] = seffsun[4] + a[4] * (workingTable.at[index, "star_teff"] - 5780.0) + b[4] * (workingTable.at[index,"star_teff"] - 5780.0) ** 2 + c[4] * (workingTable.at[index,"star_teff"] - 5780.0) ** 3 + d[4] * (workingTable.at[index,"star_teff"] - 5780.0) ** 4
    workingTable.at[index, "st_five_Earth"] = seffsun[5] + a[5] * (workingTable.at[index, "star_teff"] - 5780.0) + b[5] * (workingTable.at[index,"star_teff"] - 5780.0) ** 2 + c[5] * (workingTable.at[index,"star_teff"] - 5780.0) ** 3 + d[5] * (workingTable.at[index,"star_teff"] - 5780.0) ** 4

outputDirectory = pathlib.Path("./data/star_pickles_enriched")
outputDirectory.mkdir(exist_ok=True)
workingTable.to_pickle(outputDirectory / pkl.name)
