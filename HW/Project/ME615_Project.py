from tkinter import filedialog, Tk
from household import Household

def file_paths():
    '''This will create a list of file paths.'''

    root = Tk()
    root.filename = filedialog.askopenfilenames(initialdir="C:/Users\MTH\Documents\Gradwork\FUEL\olivier_files",
                                                 title="Select file",
                                                 filetypes=(("csv files","*.csv"),("all files","*.*")))
    return root.filename


def find_household_number(filepaths):
    '''This will find the household numbers associated with each file.'''
    paths = [filepaths]
    for file in paths:
        filename = file.split("/")
    return filename[-1].split("_")[1]


def stoves_in_study():
    '''This function will return a dictionary with the stove name and the households it is used in.'''

    stove_types = {}

    paths = file_paths()
    for f in paths:
        household = find_household_number(f)
        data = Household(f)
        stoves = data.stove_types()
        for s in stoves:
            if s.lower() not in stove_types:
                info = {s.lower(): [household]}
                stove_types.update(info)
            else:
                stove_types[s.lower()].append(household)
    return stove_types

def stove_prevelance():
    '''This function will return the number of time a stove was used in a study.'''

    stoves = stoves_in_study()
    names = [key for key in stoves]
    households = [stoves[key] for key in stoves]
    usage = [[names[i], len(num)] for i, num in enumerate(households)]
    return usage


if __name__ == "__main__":
    #print(stoves_in_study())
    print(stove_prevelance())

