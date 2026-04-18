
# Python STD
# import enum
import math
import re

# Data
import numpy as np
import pandas as pd
import scipy as sp

from numba import njit

# Machine Learning

# Image Processing / Computer Vision

# Optimization

# Auxiliary

# Visualization
import matplotlib.pyplot as plt

# Miscellaneous
from enum import auto, Enum, unique

# Typing
from typing import Any, Callable, Dict, List, Literal, Optional, Self, Set, Tuple, Union

# Local Packages


# Auxiliary Functions

def ParseYahooFinanceDivTable( filePath: str ) -> pd.DataFrame:
    """
    Parses a Yahoo Finance dividends table export and returns a DataFrame.
    
    Args:
        filePath (str): The path to a text file containing Yahoo Finance dividends table HTML.
    Returns:
        pd.DataFrame: A DataFrame containing the parsed dividends data with columns 'Date', 'Dividend'.
    """
    with open(filePath, 'r', encoding = 'utf-8') as hFile:
        fileText = hFile.read()

    # Match each table row date and dividend amount from the Yahoo Finance table markup.
    lMatches = re.findall(
        r"<tr[^>]*>\s*"
        r"<td[^>]*>\s*([^<]+?)\s*</td>"
        r".*?<span[^>]*>\s*([0-9]*\.?[0-9]+)\s*</span>\s*Dividend",
        fileText,
        flags = re.IGNORECASE | re.DOTALL,
    )

    if len(lMatches) == 0:
        raise ValueError(f'No dividend rows were found in file: {filePath}')

    dfDividends = pd.DataFrame(lMatches, columns = ['Date', 'Dividend'])
    dfDividends['Date'] = pd.to_datetime(dfDividends['Date'], format = '%b %d, %Y', errors = 'coerce')
    dfDividends['Dividend'] = pd.to_numeric(dfDividends['Dividend'], errors = 'coerce')
    dfDividends = dfDividends.dropna(subset = ['Date', 'Dividend']).reset_index(drop = True)

    if len(dfDividends) == 0:
        raise ValueError(f'Dividend rows were found but could not be parsed in file: {filePath}')

    return dfDividends

def ParseYahooFinanceCloseTable( filePath: str ) -> pd.DataFrame:
    """
    Parses a Yahoo Finance historical prices table export and returns a DataFrame.
    
    Args:
        filePath (str): The path to a text file containing Yahoo Finance historical prices table HTML.
    Returns:
        pd.DataFrame: A DataFrame containing the parsed table with columns
                      'Date', 'Open', 'High', 'Low', 'Close', 'AdjClose', and 'Volume'.
    """
    with open(filePath, 'r', encoding = 'utf-8') as hFile:
        fileText = hFile.read()

    lRows = re.findall(r'<tr[^>]*>(.*?)</tr>', fileText, flags = re.IGNORECASE | re.DOTALL)
    lData = []

    for rowText in lRows:
        lCells = re.findall(r'<td[^>]*>(.*?)</td>', rowText, flags = re.IGNORECASE | re.DOTALL)

        # Regular historical-price rows contain Date, Open, High, Low, Close, Adj Close, Volume.
        if len(lCells) < 7:
            continue

        lCellText = [re.sub(r'<[^>]+>', '', cellText).strip() for cellText in lCells]

        tsDate = pd.to_datetime(lCellText[0], format = '%b %d, %Y', errors = 'coerce')
        valOpen = pd.to_numeric(lCellText[1].replace(',', ''), errors = 'coerce')
        valHigh = pd.to_numeric(lCellText[2].replace(',', ''), errors = 'coerce')
        valLow = pd.to_numeric(lCellText[3].replace(',', ''), errors = 'coerce')
        valClose = pd.to_numeric(lCellText[4].replace(',', ''), errors = 'coerce')
        valAdjClose = pd.to_numeric(lCellText[5].replace(',', ''), errors = 'coerce')
        valVolume = pd.to_numeric(lCellText[6].replace(',', ''), errors = 'coerce')

        if pd.isna(tsDate) or pd.isna(valOpen) or pd.isna(valHigh) or pd.isna(valLow) or pd.isna(valClose) or pd.isna(valAdjClose) or pd.isna(valVolume):
            continue

        lData.append((tsDate, valOpen, valHigh, valLow, valClose, valAdjClose, valVolume))

    if len(lData) == 0:
        raise ValueError(f'No historical price rows were found in file: {filePath}')

    dfClose = pd.DataFrame(lData, columns = ['Date', 'Open', 'High', 'Low', 'Close', 'AdjClose', 'Volume'])

    return dfClose.reset_index(drop = True)


