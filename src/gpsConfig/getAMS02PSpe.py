import pandas as pd

originalDataFrame = pd.read_csv("2017AMS02HeilumSpe.csv", sep=",", header=None)
sortedDataFrame = originalDataFrame.sort_values(by=0)
print(sortedDataFrame)

gpsColumn = []

for i in range(0, len(sortedDataFrame)):
    gpsColumn.append("/gps/hist/point")
    MeVenergy  = sortedDataFrame.iloc[i, 0] * 1000 # Convert from MeV to GeV
    # flux = sortedDataFrame.iloc[i, 1] / (sortedDataFrame.iloc[i, 0] ** 2.7) 
    flux = sortedDataFrame.iloc[i, 1]

    sortedDataFrame.iloc[i, 0] = MeVenergy
    sortedDataFrame.iloc[i, 1] = flux

sortedDataFrame.insert(0, 'header', gpsColumn)
print(sortedDataFrame)

sortedDataFrame.to_csv("2017AMS02HeilumEnergySpectrum.csv", header=None, index=None, sep=" ")