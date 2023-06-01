import sys
import numpy
import pandas
import pyvo

from astropy import units
from astropy import constants
from astropy.visualization import quantity_support

from uio.utility.databases import tap


if len(sys.argv) < 2:
    raise SystemExit("No star name provided")

if len(sys.argv) > 2:
    raise SystemExit("Too many arguments")

queryStar = sys.argv[1]
print(f"Processing star name: {queryStar}")
queryStarNoSpace = queryStar.strip().replace(" ", "-")

quantity_support()

mass_coeff = constants.M_jup / constants.M_earth
radius_coeff = constants.R_jup / constants.R_earth

serviceNameNASA = "NASA"
if not tap.services.get(serviceNameNASA):
    raise SystemExit(f"The list of services doesn't have {serviceNameNASA}")
serviceEndpointNASA = tap.getServiceEndpoint(serviceNameNASA)

dictForCheckingPlanet = {
    "pl_massj": "mass",
    "pl_massjerr2": "mass_error_min",
    "pl_massjerr1": "mass_error_max",
    "pl_massjlim": "pl_massjlim",
    "pl_radj": "radius",
    "pl_radjerr2": "radius_error_min",
    "pl_radjerr1": "radius_error_max",
    "pl_radjlim": "pl_radjlim",
    "rv_flag": "rv_flag",
    "tran_flag": "tran_flag",
    "ttv_flag": "ttv_flag",
    "ima_flag": "ima_flag",
    "pl_name": "granule_uid"
}

dictForCheckingStar = {
    "st_spectype": "star_spec_type",
    "hostname": "star_name"
}

dictForChecking = {}
dictForChecking.update(dictForCheckingPlanet.items())
dictForChecking.update(dictForCheckingStar.items())

# NASA

print("\n[1/3] Querying NASA...")

resultsNASA = tap.queryService(
    serviceEndpointNASA,
    " ".join((
        "SELECT DISTINCT hostname, pl_name",
        f"FROM ps",
        f"WHERE hostname = '{queryStar}'"
    ))
)
planetsFound = len(resultsNASA)
print(f"Total planets found in {queryStar}: {planetsFound}")
if planetsFound == 0:
    print("No planets found!")
    raise SystemExit(0)

workingTableExoplanets = resultsNASA.to_table().to_pandas(index="pl_name")

for valueToAdd in dictForChecking.keys():
    if valueToAdd in ["hostname", "pl_name"]:
        continue
    workingTableExoplanets[valueToAdd] = (
        numpy.array(numpy.NaN, dtype=float)
        if valueToAdd not in tap.services[serviceNameNASA][
            "parameters-that-are-strings"
        ]
        else numpy.array(numpy.NaN, dtype=str)
    )

print("\n[2/3] Getting parameters from NASA...")

# print(workingTableExoplanets)

for index, row in workingTableExoplanets.iterrows():
    systemName = row["hostname"]
    print(f"Iterating {index}, star name: {queryStar}")
    for valueToAdd in dictForChecking.keys():
        if (
            valueToAdd in ["hostname", "pl_name"]
            or valueToAdd.endswith("err1")
            or valueToAdd.endswith("err2")
        ):
            continue
        else:
            valueFromNASA = tap.getParameterFromNASA(
                systemName,
                index.replace("'", ""),
                valueToAdd
            )
            if valueFromNASA:
                workingTableExoplanets.at[index, valueToAdd] = valueFromNASA
                # if it's a parameter with min/max errors, get them
                if valueToAdd in (
                    tap.services[serviceNameNASA][
                        "parameters-that-have-errors"
                    ]
                ):
                    print(f"---\nTrying to get min/max errors of {valueToAdd}")
                    valueFromNASAerr1, valueFromNASAerr2 = tap.getParameterErrorsFromNASA(
                        systemName,
                        index.replace("'", ""),
                        valueToAdd
                    )
                    if valueFromNASAerr1:
                        workingTableExoplanets.at[index, f"{valueToAdd}err1"] = valueFromNASAerr1
                    if valueFromNASAerr2:
                        workingTableExoplanets.at[index, f"{valueToAdd}err2"] = valueFromNASAerr2
                    print("---")

print("\n[3/3] Cross-checking with PADC, checking...")

workingTableExoplanets["flag_mE"] = numpy.array(False, dtype=bool)
workingTableExoplanets["flag_rE"] = numpy.array(False, dtype=bool)
workingTableExoplanets["mass_detection_type"] = numpy.array("", dtype=str)
workingTableExoplanets["radius_detection_type"] = numpy.array("", dtype=str)

# checking mass and radius

enrichedCount = 0
for index, row in workingTableExoplanets.iterrows():
    enriched = False
    # mass
    if pandas.isna(row["pl_massj"]):
        valueFromPADC = tap.getParameterFromPADC(
            index.replace("'", ""),
            tap.mappings["NASA-to-PADC"]["planets"]["pl_massj"]
        )
        if valueFromPADC:
            enriched = True
            workingTableExoplanets.at[index, "pl_massj"] = valueFromPADC

            valErrorMin, valErrorMax = tap.getParameterErrorsFromPADC(
                index.replace("'", ""),
                tap.mappings["NASA-to-PADC"]["planets"]["pl_massj"]
            )
            workingTableExoplanets.at[index, "pl_massjerr1"] = valErrorMin
            workingTableExoplanets.at[index, "pl_massjerr2"] = valErrorMax

            workingTableExoplanets.at[index, "flag_mE"] = True
            valueFromPADC1 = tap.getParameterFromPADC(
                index.replace("'", ""),
                "mass_detection_type"
            )
            workingTableExoplanets.at[index, "mass_detection_type"] = valueFromPADC1

    # radius
    if pandas.isna(row["pl_radj"]):
        valueFromPADC = tap.getParameterFromPADC(
            index.replace("'", ""),
            tap.mappings["NASA-to-PADC"]["planets"]["pl_radj"]
        )
        if valueFromPADC:
            enriched = True
            workingTableExoplanets.at[index, "pl_radj"] = valueFromPADC

            valErrorMin, valErrorMax = tap.getParameterErrorsFromPADC(
                index.replace("'", ""),
                tap.mappings["NASA-to-PADC"]["planets"]["pl_radj"]
            )
            workingTableExoplanets.at[index, "pl_radjerr1"] = valErrorMin
            workingTableExoplanets.at[index, "pl_radjerr2"] = valErrorMax

            workingTableExoplanets.at[index, "flag_rE"] = True
            valueFromPADC1 = tap.getParameterFromPADC(
                index.replace("'", ""),
                "radius_detection_type"
            )
            workingTableExoplanets.at[index, "radius_detection_type"] = valueFromPADC1

    # spectral type
    if pandas.isna(row["st_spectype"]):
        valueFromPADC = tap.getParameterFromPADC(
            index.replace("'", ""),
            tap.mappings["NASA-to-PADC"]["stars"]["st_spectype"]
        )
        if valueFromPADC:
            enriched = True
            workingTableExoplanets.at[index, "st_spectype"] = valueFromPADC

    if enriched:
        enrichedCount += 1

print(f"Enriched planets count: {enrichedCount}")

# now rename the table to its historically usual columns names

workingTableExoplanets.index.rename("granule_uid", inplace=True)
workingTableExoplanets.rename(
    columns=(dictForCheckingPlanet | dictForCheckingStar),
    inplace=True
)

print()
print(workingTableExoplanets)

workingTableExoplanets.to_pickle(f"./data/star_pickles/{queryStar}.pkl")
