"""
This script can be used to create test data for define wind. To use it, do the
following:

- Run the model you want to create test data using
- Run windsave2table
- Run the script to create a {model}.grid.txt file for all the .master and .grad
  files found in the working directory
- Move and rename files if appropriate

Example:
    $ py --grid-only my_model
    $ windsave2table my_model
    $ python create_grid_files.py

Notes:
    - May want to update script to take in the model as a command line
      argument instead of using glob
"""

import pandas
import glob

master_files = sorted(glob.glob("*.master.txt"))
model_names = [file.replace(".master.txt", "") for file in master_files]

# There are the parameters we'll mush into one file. Note that some parameters
# may not make sense for all models. There is a catch later in the script to
# only use the parameters it finds in the windsave2table output
master_parameters = [
    "i",
    "j",
    "x",
    "z",
    "xcen",
    "zcen",
    "r",
    "rcen",   # 1d
    "r_cen",  # 2d :)
    "theta",
    "theta_cen",
    "inwind",
    "v_x",
    "v_y",
    "v_z",
    "vol",
    "rho",
    "ne",
    "t_e",
    "t_r",
    "h1",
    "c4",
]
gradient_parameters = [
    "dv_x_dx",
    "dv_y_dx",
    "dv_z_dx",
    "dv_x_dy",
    "dv_y_dy",
    "dv_z_dy",
    "dv_x_dz",
    "dv_y_dz",
    "dv_z_dz",
    "div_v",
    "dvds_max",
    # "dvds_ave",
    "gamma",
]

# do it for all models found in working directory
for model, file in zip(model_names, master_files):
    master_df = pandas.read_csv(file, delim_whitespace=True)
    gradient_df = pandas.read_csv(f"{model}.gradient.txt", delim_whitespace=True)

    if model == "star":
        print(model, master_df.columns, gradient_df.columns)

    with open(f"{model}.grid.txt", "w", encoding="utf-8") as file_out:
        # Create header row
        line = "# "
        for parameter in master_parameters + gradient_parameters:
            if parameter in master_df.columns or parameter in gradient_df.columns:
                line += f"{parameter} "
        line += "\n"
        file_out.write(line)

        # Create contents of file. Write out the parameters defined above for
        # each cell
        for (_1, master_row), (_2, gradient_row) in zip(
            master_df.iterrows(), gradient_df.iterrows()
        ):
            line = ""
            for parameter in master_parameters:
                if parameter not in master_df:
                    continue
                if parameter in ["i", "j", "inwind"]:
                    line += f"{int(master_row[parameter]):d} "
                else:
                    line += f"{master_row[parameter]:9.8e} "
            for parameter in gradient_parameters:
                if parameter not in gradient_df:
                    continue
                line += f"{gradient_row[parameter]:9.8e} "
            line += "\n"
            file_out.write(line)
