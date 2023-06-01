import os
import sys
import pandas
import shutil
import pathlib
from glob import glob

from uio.utility.files import pickle

if len(sys.argv) != 2:
    raise SystemExit(
        "[ERROR] You need to provide path to the folder with pickles"
    )
folder = pathlib.Path(sys.argv[1])
print(f"Folder: {folder.absolute()}")
if not folder or not folder.is_dir():
    raise SystemExit(
        "[ERROR] Provided folder does not exist"
    )

pickles = glob(os.path.join(folder, "*.pkl"))
if len(pickles) == 0:
    raise SystemExit("[ERROR] There are no pickles in this folder")

pathlib.Path.mkdir(folder / "minimum2withmassandradius", exist_ok=True)
for pkl in pickles:
    pth = pathlib.Path(pkl)
    workingTable = pickle.openPickleAsPandasTable(pth)
    numberOfPlanetsWithMassAndRadius = len(
        workingTable.query("mass.notna() & radius.notna()")
    )
    print(
        " ".join((
            f"{pth.stem}: {numberOfPlanetsWithMassAndRadius}",
            "planets with mass and radius"
        ))
    )
    if numberOfPlanetsWithMassAndRadius > 1:
        shutil.copy(
            folder / pth.name,
            folder / "minimum2withmassandradius" / pth.name
        )
