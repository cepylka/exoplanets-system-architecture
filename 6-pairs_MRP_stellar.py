import pathlib
import pandas
import numpy

from uio.utility.files import pickle

# original table
tbl = pickle.openPickleAsPandasTable("./data/all_my_systems.pkl")

# sorted table
tblSorted = tbl.sort_values(["star_name", "period"], ascending=[True, True])
# display(tblSorted)

tblPairs = pandas.DataFrame(
    {
        "granule_uid-1": pandas.Series(dtype="string"),
        "radius-1": pandas.Series(dtype="float"),
        "radius_error_min-1": pandas.Series(dtype="float"),
        "radius_error_max-1": pandas.Series(dtype="float"),
        "mass-1": pandas.Series(dtype="float"),
        "mass_error_min-1": pandas.Series(dtype="float"),
        "mass_error_max-1": pandas.Series(dtype="float"),
        "semi_major_axis-1": pandas.Series(dtype="float"),
        "semi_major_axis_error_min-1": pandas.Series(dtype="float"),
        "semi_major_axis_error_max-1": pandas.Series(dtype="float"),
        "period-1": pandas.Series(dtype="float"),
        "period_error_min-1": pandas.Series(dtype="float"),
        "period_error_max-1": pandas.Series(dtype="float"),
        "granule_uid-2": pandas.Series(dtype="string"),
        "radius-2": pandas.Series(dtype="float"),
        "radius_error_min-2": pandas.Series(dtype="float"),
        "radius_error_max-2": pandas.Series(dtype="float"),
        "mass-2": pandas.Series(dtype="float"),
        "mass_error_min-2": pandas.Series(dtype="float"),
        "mass_error_max-2": pandas.Series(dtype="float"),
        "semi_major_axis-2": pandas.Series(dtype="float"),
        "semi_major_axis_error_min-2": pandas.Series(dtype="float"),
        "semi_major_axis_error_max-2": pandas.Series(dtype="float"),
        "period-2": pandas.Series(dtype="float"),
        "period_error_min-2": pandas.Series(dtype="float"),
        "period_error_max-2": pandas.Series(dtype="float"),
        "star_teff": pandas.Series(dtype="float"),
        "star_radius": pandas.Series(dtype="float"),
        "star_mass": pandas.Series(dtype="float"),
        "star_metallicity": pandas.Series(dtype="float"),
        "star_age": pandas.Series(dtype="float"),
        "st_rotp": pandas.Series(dtype="float"),
        "sy_snum": pandas.Series(dtype="float"),
        "cb_flag": pandas.Series(dtype="float"),
        "sy_pnum": pandas.Series(dtype="float"),
        "pl_radjlim": pandas.Series(dtype="float"),
        "pl_massjlim": pandas.Series(dtype="float"),
        "mass_detection_type": pandas.Series(dtype="string")

    }
)

frames = []

grouped = tblSorted.groupby("star_name")
for name, group in grouped:
    currentSystemSize = len(group)
    # print(f"Planets: {currentSystemSize}\n")
    for i in range(currentSystemSize):
        # check if next iteration exceeds the available range
        if (i + 1) == currentSystemSize:
            break
        # make a pair
        row = pandas.DataFrame(
            {
                "granule_uid-1": group.index[i],
                "radius-1": group.iloc[i]["radius"],
                "radius_error_min-1": group.iloc[i]["radius_error_min"],
                "radius_error_max-1": group.iloc[i]["radius_error_max"],
                "mass-1": group.iloc[i]["mass"],
                "mass_error_min-1": group.iloc[i]["mass_error_min"],
                "mass_error_max-1": group.iloc[i]["mass_error_max"],
                "semi_major_axis-1": group.iloc[i]["semi_major_axis"],
                "semi_major_axis_error_min-1": group.iloc[i]["semi_major_axis_error_min"],
                "semi_major_axis_error_max-1": group.iloc[i]["semi_major_axis_error_max"],
                "period-1": group.iloc[i]["period"],
                "period_error_min-1": group.iloc[i]["period_error_min"],
                "period_error_max-1": group.iloc[i]["period_error_max"],
                "granule_uid-2": group.index[i+1],
                "radius-2": group.iloc[i+1]["radius"],
                "radius_error_min-2": group.iloc[i+1]["radius_error_min"],
                "radius_error_max-2": group.iloc[i+1]["radius_error_max"],
                "mass-2": group.iloc[i+1]["mass"],
                "mass_error_min-2": group.iloc[i+1]["mass_error_min"],
                "mass_error_max-2": group.iloc[i+1]["mass_error_max"],
                "semi_major_axis-2": group.iloc[i+1]["semi_major_axis"],
                "period-2": group.iloc[i+1]["period"],
                "period_error_min-2": group.iloc[i+1]["period_error_min"],
                "period_error_max-2": group.iloc[i+1]["period_error_max"],
                "semi_major_axis_error_min-2": group.iloc[i+1]["semi_major_axis_error_min"],
                "semi_major_axis_error_max-2": group.iloc[i+1]["semi_major_axis_error_max"],
                "star_name": group.iloc[i]["star_name"],
                "star_teff": group.iloc[i]["star_teff"],
                "star_radius": group.iloc[i]["star_radius"],
                "star_mass": group.iloc[i]["star_mass"],
                "star_metallicity": group.iloc[i]["star_metallicity"],
                "star_age": group.iloc[i]["star_age"],
                "st_rotp": group.iloc[i]["st_rotp"],
                "sy_snum": group.iloc[i]["sy_snum"],
                "sy_snum": group.iloc[i]["sy_snum"],
                "sy_pnum": group.iloc[i]["sy_pnum"],
                "pl_radjlim": group.iloc[i]["pl_radjlim"],
                "pl_massjlim":group.iloc[i]["pl_massjlim"],
                "mass_detection_type":group.iloc[i]["mass_detection_type"],
                "ttv_flag":group.iloc[i]["ttv_flag"],
                "rv_flag":group.iloc[i]["rv_flag"],
                "tran_flag":group.iloc[i]["tran_flag"],
                "ima_flag":group.iloc[i]["ima_flag"]

            },
            index=[f"{group.iloc[i]['star_name']}-{i+1}"]
        )
        # add that pair to the resulting table
        if (
            pandas.isnull(group.iloc[i]["pl_radjlim"])
            and
            pandas.isnull(group.iloc[i+1]["pl_radjlim"])
            and
            pandas.isnull(group.iloc[i]["pl_massjlim"])
            and
            pandas.isnull(group.iloc[i+1]["pl_massjlim"])
        ):
            frames.append(row)

tblPairs = pandas.concat(frames)  # .sort_index()
# display(tblPairs)
tblPairs.to_pickle("./data/MRP_data_sets.pkl")
