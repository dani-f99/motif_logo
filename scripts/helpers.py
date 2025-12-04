import pandas as pd
import numpy as np
import json
import os

# Reading information from json file. Used to extract the parameters from the config.json.
def read_json(path:str = "config.json") -> dict:
    """
    path : str -> path of the json file
    """

    with open('config.json') as config:
        config_f = json.load(config)

    return config_f

# Creating folder according to the `congif.json` paths and program scheme
def create_folders():
    config_dict = read_json()
    main_folders = [config_dict["input_folder"], config_dict["output_folder"]]
    subf_output = [config_dict["output_folder"] + i for i in ["//motif_figure", "//motif_data"]]

    for folder in main_folders + subf_output:
        if os.path.exists(folder) is False:
            os.mkdir(folder)
            print(f"> folder `{folder}` was created.")
        else:
            print(f"> folder `{folder}` exists, continuing.")


# Importing