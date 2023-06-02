from pathlib import Path
import pandas


# original table
tbl = pandas.read_pickle("./data/all_my_systems.pkl")
# tbl = tbl[~tbl['mass_detection_type'].isin(['Theoretical'])]
print(f"Planets before filtering: {len(tbl.index)}")
# keep only systems with more that 2 planets
tbl = tbl.groupby(["star_name"]).filter(lambda x: len(x) > 2)
print(f"Planets after filtering: {len(tbl.index)}")
# sorted table
tblSorted = tbl.sort_values(["star_name", "period"], ascending=[True, True])
# display(tblSorted)

tblTriples = pandas.DataFrame(
    {
        "granule_uid-1": pandas.Series(dtype="float"),
        "granule_uid-2": pandas.Series(dtype="float"),
        "granule_uid-3": pandas.Series(dtype="float"),
        "semi_major_axis-1": pandas.Series(dtype="float"),
        "semi_major_axis-2": pandas.Series(dtype="float"),
        "semi_major_axis-3": pandas.Series(dtype="float"),
        "semi_major_axis_error_min-1": pandas.Series(dtype="float"),
        "semi_major_axis_error_max-1": pandas.Series(dtype="float"),
        "semi_major_axis_error_min-2": pandas.Series(dtype="float"),
        "semi_major_axis_error_max-2": pandas.Series(dtype="float"),
        "semi_major_axis_error_min-3": pandas.Series(dtype="float"),
        "semi_major_axis_error_max-3": pandas.Series(dtype="float"),
        "period-1": pandas.Series(dtype="float"),
        "period-2": pandas.Series(dtype="float"),
        "period-3": pandas.Series(dtype="float"),
        "period_error_min-1": pandas.Series(dtype="float"),
        "period_error_max-1": pandas.Series(dtype="float"),
        "period_error_min-2": pandas.Series(dtype="float"),
        "period_error_max-2": pandas.Series(dtype="float"),
        "period_error_min-3": pandas.Series(dtype="float"),
        "period_error_max-3": pandas.Series(dtype="float"),
        "radius-1": pandas.Series(dtype="float"),
        "radius_error_min-1": pandas.Series(dtype="float"),
        "radius_error_max-1": pandas.Series(dtype="float"),
        "mass-1": pandas.Series(dtype="float"),
        "mass_error_min-1": pandas.Series(dtype="float"),
        "mass_error_max-1": pandas.Series(dtype="float"),
        "semi_major_axis-1": pandas.Series(dtype="float"),
        "semi_major_axis_error_min-1": pandas.Series(dtype="float"),
        "semi_major_axis_error_max-1": pandas.Series(dtype="float"),
        "radius-2": pandas.Series(dtype="float"),
        "radius_error_min-2": pandas.Series(dtype="float"),
        "radius_error_max-2": pandas.Series(dtype="float"),
        "mass-2": pandas.Series(dtype="float"),
        "mass_error_min-2": pandas.Series(dtype="float"),
        "mass_error_max-2": pandas.Series(dtype="float"),
        "radius-3": pandas.Series(dtype="float"),
        "radius_error_min-3": pandas.Series(dtype="float"),
        "radius_error_max-3": pandas.Series(dtype="float"),
        "mass-3": pandas.Series(dtype="float"),
        "mass_error_min-3": pandas.Series(dtype="float"),
        "mass_error_max-3": pandas.Series(dtype="float"),
    }
)

frames = []

for name, group in tblSorted.groupby("star_name"):
    group = group.sort_values(["period"], ascending=[True])
    currentSystemSize = len(group)
    # print(f"Number of planets in {name}: {currentSystemSize}")
    for i in range(currentSystemSize):
        # check if next triple exceeds the available range
        if (i + 2) == currentSystemSize:
            break
        # make a triple
        row = pandas.DataFrame(
            {
                "granule_uid-1": group.index[i],
                "semi_major_axis-1": group.iloc[i]["semi_major_axis"],
                "semi_major_axis_error_min-1": group.iloc[i]["semi_major_axis_error_min"],
                "semi_major_axis_error_max-1": group.iloc[i]["semi_major_axis_error_max"],
                "period-1": group.iloc[i]["period"],
                "period_error_min-1": group.iloc[i]["period_error_min"],
                "period_error_max-1": group.iloc[i]["period_error_max"],
                "granule_uid-2": group.index[i+1],
                "semi_major_axis-2": group.iloc[i+1]["semi_major_axis"],
                "semi_major_axis_error_min-2": group.iloc[i+1]["semi_major_axis_error_min"],
                "semi_major_axis_error_max-2": group.iloc[i+1]["semi_major_axis_error_max"],
                "radius-1": group.iloc[i]["radius"],
                "radius_error_min-1": group.iloc[i]["radius_error_min"],
                "radius_error_max-1": group.iloc[i]["radius_error_max"],
                "mass-1": group.iloc[i]["mass"],
                "mass_error_min-1": group.iloc[i]["mass_error_min"],
                "mass_error_max-1": group.iloc[i]["mass_error_max"],
                "period-2": group.iloc[i+1]["period"],
                "period_error_min-2": group.iloc[i+1]["period_error_min"],
                "period_error_max-2": group.iloc[i+1]["period_error_max"],
                "radius-2": group.iloc[i+1]["radius"],
                "radius_error_min-2": group.iloc[i+1]["radius_error_min"],
                "radius_error_max-2": group.iloc[i+1]["radius_error_max"],
                "mass-2": group.iloc[i+1]["mass"],
                "mass_error_min-2": group.iloc[i+1]["mass_error_min"],
                "mass_error_max-2": group.iloc[i+1]["mass_error_max"],
                "granule_uid-3": group.index[i+2],
                "semi_major_axis-3": group.iloc[i+2]["semi_major_axis"],
                "semi_major_axis_error_min-3": group.iloc[i+2]["semi_major_axis_error_min"],
                "semi_major_axis_error_max-3": group.iloc[i+2]["semi_major_axis_error_max"],
                "period-3": group.iloc[i+2]["period"],
                "period_error_min-3": group.iloc[i+2]["period_error_min"],
                "period_error_max-3": group.iloc[i+2]["period_error_max"],
                "radius-3": group.iloc[i+2]["radius"],
                "radius_error_min-3": group.iloc[i+2]["radius_error_min"],
                "radius_error_max-3": group.iloc[i+2]["radius_error_max"],
                "mass-3": group.iloc[i+2]["mass"],
                "mass_error_min-3": group.iloc[i+2]["mass_error_min"],
                "mass_error_max-3": group.iloc[i+2]["mass_error_max"],
                "star_name": group.iloc[i]["star_name"],
                "star_teff": group.iloc[i]["star_teff"],
                "star_radius": group.iloc[i]["star_radius"],
                "star_mass": group.iloc[i]["star_mass"],
                "star_metallicity": group.iloc[i]["star_metallicity"],
                "star_age": group.iloc[i]["star_age"],
                "st_rotp": group.iloc[i]["st_rotp"],
                "sy_snum": group.iloc[i]["sy_snum"],
                "sy_pnum": group.iloc[i]["sy_pnum"],
                "ima_flag":group.iloc[i]["ima_flag"]
            },

            index=[f"{group.iloc[i]['star_name']}-{i+1}"]
        )
        frames.append(row)

tblTriples = pandas.concat(frames)

tblTriples.to_pickle("./data/MRP_data_sets_triples.pkl")
