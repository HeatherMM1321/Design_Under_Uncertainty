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

if __name__ == "__main__":

    study_stoves, study_fuels = info()

    df_usage = pd.DataFrame(columns=['Day', 'Household', 'telia usage (mins)', 'om30 usage (mins)',
                                     '3stone usage (mins)', 'malgchch usage (mins)',
                                     'firewood usage (kg)', 'charcoal usage (kg)', 'lpg usage (kg)'])
    # look through each of the files
    for i, file in enumerate(filepaths):
        df_stoves, stoves, fuels, hh_id = convert.reformat_olivier_files(file)
        data = household.Household(df_stoves, stoves, fuels, hh_id)

        # getting household cooking and fuel use info
        stove_data = data.cooking_duration()
        fuel_data = data.fuel_usage()

        days = round(data.study_duration.total_seconds()/86400)

        for j in range(days):
            info = {'Day': j+1, 'Household': hh_id}
            for s in stoves:
                time = stove_data[s].values[j]
                if s == '3pierres' or s == '3':
                    s = '3stone'
                info.update({s + ' usage (mins)': time})
                # else:
                #     info.update({stv + ' usage (mins)': 'NA'})
            for f in fuels:
                time = fuel_data[f].values[j]
                info.update({f + ' usage (kg)': time})
            # append the dataframe with the complete information for that household on that day
            df_usage = df_usage.append(info, ignore_index=True)

    # writing info to csv file
    df_usage.to_csv(r'usage_data.csv')



