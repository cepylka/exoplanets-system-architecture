import pathlib
import pandas

from uio.utility.files import pickle

pickle.mergePickles(
    "./merging/",
    "./data/all_my_systems.pkl"
)
