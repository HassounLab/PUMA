import xlsxwriter
import numpy as np
import pandas as pd


def write_data(result, filename="file_name.xlsx", sheetname='sheet_name', header=None):
    """
    write out data to excel sheet;  better change to csv.
    :param result: data to write out
    :param filename
    :param sheetname
    :return:
    None
    """
    workbook = xlsxwriter.Workbook(filename)
    worksheet = workbook.add_worksheet(sheetname)
    col = 0
    worksheet.write_row(0, col, header)
    for row, data in enumerate(result):
        worksheet.write_row(row+1, col, data)
    workbook.close()



def load_samples_from_csv(filename_samples):
    """
    Load samples after pathway activity prediction
    :param filename_samples:
    :return:
    """
    samples = pd.read_csv(filename_samples)
    samples.drop(['Unnamed: 0'], axis=1, inplace=True)
    samples.head()
    return np.array(samples)


def load_samples_from_excel(filename_samples):
    """
    Load samples after pathway activity prediction
    :param filename_samples:
    :return:

    """
    samples = pd.read_excel(filename_samples, header=None)
    samples.head()
    return np.array(samples)


