#!/usr/bin/python3

import sys
from osgeo import gdal

pszFilename = sys.argv[1]
if pszFilename:
    hDataset = gdal.Open( pszFilename, gdal.GA_ReadOnly )

    for i in range(hDataset.RasterCount):
        hBand = hDataset.GetRasterBand(i+1)
        (lo, hi, avg, std) = hBand.ComputeStatistics(False)
#        print(f"Band {i} min={lo} max={hi} avg={avg} dev={std}")

        hist = hBand.GetHistogram(min=-0.5, max=hi+0.5, buckets=int(hi+1), approx_ok = False)
#        print(hist)