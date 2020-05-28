from tkinter import filedialog, Tk
import sys
sys.path.insert(-1, 'C:\\Users\\MTH\\Documents\\Gradwork\\Spring2020\\Software_development\\FUEL-package')
from FUEL import household
from FUEL import olivier_file_convert as convert
import pandas as pd
#
# def file_paths():
#     '''This will create a list of file paths.'''
#
#
#     root = Tk()
#     root.filename = filedialog.askopenfilenames(initialdir="C:/Users\MTH\Documents\Gradwork\FUEL\olivier_files",
#                                                  title="Select file",
#                                                  filetypes=(("csv files","*.csv"),("all files","*.*")))
#     return root.filename


filepaths = ['HH_38_2018-08-26_15-01-40_processed_v3.csv',
             'HH_44_2018-08-17_13-49-22_processed_v2.csv',
             'HH_141_2018-08-17_17-50-31_processed_v2.csv',
             'HH_318_2018-08-25_18-35-07_processed_v2.csv',
             'HH_319_2018-08-25_19-27-32_processed_v2.csv',
             'HH_326_2018-08-25_17-52-16_processed_v2.csv',
             'HH_345_2018-08-25_15-52-57_processed_v2.csv',
             'HH_371_2018-08-17_15-31-52_processed_v2.csv'
             ]


def find_household_number():
    '''This will find the household numbers associated with each file.'''

    households = []

    for file in filepaths:
        number = file.split("_")[1]
        households.append(number)
    return households


def info():
    '''List of all stoves and fuels in study.'''

    all_stoves = []
    all_fuels = []

    for file in filepaths:
        df_stoves, stoves, fuels = convert.reformat_olivier_files(file)
        data = household.Household(df_stoves, stoves, fuels)
        for s in data.stoves:
            if s == '3pierres' or s == '3':
                s = '3stone'
            if s not in all_stoves:
                all_stoves.append(s)
        for f in data.fuels:
            if f not in all_fuels:
                all_fuels.append(f)

    return all_stoves, all_fuels


if __name__ == "__main__":

    households = find_household_number()
    study_stoves, study_fuels = info()

    df_usage = pd.DataFrame(columns=['Day', 'Household', 'telia usage (mins)', 'om30 usage (mins)',
                                     '3stone usage (mins)', 'malgchch usage (mins)',
                                     'firewood usage (kg)', 'charcoal usage (kg)', 'lpg usage (kg)'])
    # look through each of the files
    for i, house in enumerate(households):
        df_stoves, stoves, fuels = convert.reformat_olivier_files(filepaths[i])
        data = household.Household(df_stoves, stoves, fuels)

        # getting household cooking and fuel use info
        stove_data = data.cooking_duration()
        fuel_data = data.fuel_usage()

        # study is 3 days long so we will iterate through 3 times (could be more elegant)
        for j in range(3):
            info = {'Day': j+1, 'Household': house}
            # go through each of the stoves in the study and then put record how much that stove was
            # used on that day in that household
            for stv in study_stoves:
                if stv in stove_data.columns:
                    time = stove_data[stv].values[j]
                    info.update({stv + ' usage (mins)': time})
                else:
                    info.update({stv + ' usage (mins)': 0})
            # same process as with the stoves
            for fls in study_fuels:
                if fls in fuel_data.columns:
                    weight = fuel_data[fls].values[j]
                    info.update({fls + ' usage (kg)': weight})
                else:
                    info.update({fls + ' usage (kg)': 0})

            # append the dataframe with the complete information for that household on that day
            df_usage = df_usage.append(info, ignore_index=True)

    # writing info to csv file
    df_usage.to_csv(r'usage_data.csv')



