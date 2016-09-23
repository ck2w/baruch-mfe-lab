#!/usr/bin/anaconda3/bin/python

import csv
import math
import numpy as np

#RETURN_TYPE="PCT"
RETURN_TYPE="LOG"

with open('financials2016-8.csv') as csvfile:
    dataReader = csv.reader(csvfile, delimiter=',')
    rowCount = 0
    colCount = 0
    timeSeries = []
    for row in dataReader:
        if rowCount == 0:
            colCount = len(row)-1 # minus the data column
            for i in range(colCount):
                timeSeries.append([])
        else:
            for i in range(colCount):
                timeSeries[i].append(float(row[i+1]))
        rowCount += 1

    numSeries = len(timeSeries)
    seriesLength = len(timeSeries[0])

    returnSeriesLength = seriesLength-1
    returnSeries = []
    for i in range(numSeries):
        returnSeries.append([None]*returnSeriesLength)

    if RETURN_TYPE == "PCT":
        for i in range(numSeries):
            for j in range(returnSeriesLength):
                returnSeries[i][j] = (timeSeries[i][j]-timeSeries[i][j+1])/timeSeries[i][j+1]
    elif RETURN_TYPE == "LOG":
        for i in range(numSeries):
            for j in range(returnSeriesLength):
                returnSeries[i][j] = math.log(timeSeries[i][j]/timeSeries[i][j+1])

    mean = np.zeros(numSeries)
    for i in range(numSeries):
        s = 0
        for j in range(returnSeriesLength):
            s += returnSeries[i][j]
        s /= returnSeriesLength
        mean[i] = s

    cov = np.zeros(shape=(numSeries,numSeries))
    for i in range(numSeries):
        for j in range(numSeries):
            v = 0
            for k in range(returnSeriesLength):
                v += (returnSeries[i][k]-mean[i])*(returnSeries[j][k]-mean[j])
            v /= (returnSeriesLength-1)
            cov[i][j] = v

    print("Covariance matrix:")
    for i in range(numSeries):
        for j in range(numSeries):
            print("{0:.9f}".format(cov[i][j]), end=",")
        print("")

    cov_matrix = np.array(cov)
    L = np.linalg.cholesky(cov_matrix)
    U = np.transpose(L)
    print("U-componet in Cholesky:")
    for i in range(numSeries):
        for j in range(numSeries):
            print("{0:.9f}".format(U[i][j]), end=",")
        print("")
    
    print("L-componet in Cholesky:")
    for i in range(numSeries):
        for j in range(numSeries):
            print("{0:.9f}".format(U[j][i]), end=",")
        print("")

    returnSeriesTranspose = np.transpose(returnSeries)

    yIndex = 5 
    ySeries = np.transpose(returnSeries[yIndex])
    xSeries = np.delete(returnSeriesTranspose, yIndex, 1)
    # if include constant term, uncomment the following two lines
    oneSeries = np.array([1]*returnSeriesLength)
    xSeries = np.insert(xSeries, 0, oneSeries, axis=1)
    xSeriesTranspose = np.transpose(xSeries)
    A = np.dot(xSeriesTranspose, xSeries)
    b = np.dot(xSeriesTranspose, ySeries)
    x = np.linalg.solve(A,b)
    print("Linear regression coefficient:")
    for i in range(numSeries):
        print("{0:.9f}".format(x[i]))

    print("Linear regression residual error:")
    print(np.linalg.norm(np.dot(x,xSeriesTranspose)-ySeries))



    




