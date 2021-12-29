"""

@author: Adrien Wehrlé, Jason E. Box, GEUS (Geological Survey of Denmark and Greenland), 2019

"""

import os
import glob
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import Counter
import pickle


def load_data(file, year):

    """

    Loading PROMICE data for a given station and all or given year(s)
    into a DataFrame.


    INTPUTS:
        file: Path to the desired file containing PROMICE data [string]
        year: Year to import. If 'all', all the years are imported [int,string]

    OUTPUTS:
        promice_data: Dataframe containing PROMICE data for the desired settings [DataFrame]

    """

    # extract site name
    global site
    site = file.split("/")[-1].split("_")[0] + "_" + file.split("/")[-1].split("_")[1]

    if site[:3] == "MIT" or site[:3] == "EGP" or site[:3] == "CEN":
        site = site.split("_")[0]

    # load data
    promice_data = pd.read_csv(file, delim_whitespace=True)

    # set invalid values (-999) to nan
    promice_data[promice_data == -999.0] = np.nan

    # only keep selected year(s) if needed
    global yr
    yr = year

    if year != "all":
        promice_data = promice_data[promice_data.Year == year]
    elif isinstance(year, list):
        promice_data[promice_data.Year.isin(year)]

    if promice_data.empty:
        print("ERROR: Selected year not available")
        return

    return promice_data


def BIC_processing(promice_data, fig_save=False, fig_path=None, visualisation=False):

    """

    Processing, filtering and exclusion of ice ablation and albedo timeseries
    spanning the onset of bare ice conditions.


    INPUTS:
        promice_data: Dataset imported using load_data() [DataFrame]
        visualisation: If True, displays raw and processed DPT as well as air temperature,
                       boom height and albedo time series (default: False) [boolean]
        fig_save: If True, saves figure generated with visualisation option
                  (default: False) [boolean]
        fig_path: Path where to save figure generated with visualisation option
                  (default: None) [boolean]

    OUTPUTS:
        promice_data_proc: Processed ice ablation and albedo time
                           series [DataFrame]

    """

    def DPT_processing(df, year):

        """

        Processing of Depth Pressure Sensor (DPT) time series using PROMICE 2019-08-02
        datasets.

        INPUTS:
            promice_data: Dataframe imported using load_data() [DataFrame]
            year: Year of promice_data to process [int]

        OUTPUTS:
            DPT_proc: Processed ice ablation time series [pandas.series]
            albedo: Albedo time series associated with ice ablation [pandas.series]
            DPT_flag: Assessement of the ease to determine bare ice appearance
                      from ice ablation (0=no data, 1=, 2=low confidence,
                      3= high confidence) [int]
            albedo_flag: Determines if an albedo time series will be excluded (0)
                         or not (1) [int]
            BID: Bare Ice Day, day of bare ice appearance based on ice ablation [float]

        """

        # load ice ablation and associated time
        global z, doy
        z = df["DepthPressureTransducer_Cor(m)"].copy()
        doy = df["DayOfYear"]

        def assign_nans(days):
            if type(days) == int:
                z[days == doy] = np.nan
            else:
                for d in days:
                    z[d == doy] = np.nan

        # initialize albedo_flag
        albedo_flag = 1

        # manually process DPT measurements for a given site and year and
        # identify theoretical ice ablation onset (IAO)
        global no_proc
        no_proc = 0

        if site == "NUK_U":
            if year in [2007, 2014, 2015]:
                DPT_flag = 1
            elif year in [2008, 2009]:
                DPT_flag = 0
            elif year == 2010:
                z[doy > 203] -= 2.6
                z[doy > 204] -= 8.8
                z[doy >= 256] = np.nan
                DPT_flag = 3
                IAO = 121
            elif year == 2011:
                assign_nans(235)
                z[doy >= 236] -= 2
                z[doy >= 276] = np.nan
                DPT_flag = 3
                IAO = 156
                albedo_flag = 0
            elif year == 2012:
                z[(doy >= 260) & (z >= -3)] = np.nan
                DPT_flag = 3
                IAO = 148
            elif year == 2013:
                z[doy > 202] -= 2
                z[doy > 203] -= 5.8
                DPT_flag = 3
                IAO = 152
            elif year == 2016:
                z[doy >= 263] = np.nan
                z[(doy <= 70) & (z <= -0.5)] = np.nan
                z[(doy > 38) & (doy < 43)] = np.nan
                assign_nans(np.arange(71, 76))
                assign_nans([32, 92, 93])
                DPT_flag = 3
                IAO = 160
            elif year == 2017:
                z[(doy >= 250) & (z >= -0.3)] = np.nan
                DPT_flag = 3
                IAO = 204
            elif year == 2018:
                z[doy <= 103] = np.nan
                z[doy >= 212] -= 8
                assign_nans([123, 124, 210, 211])
                z[doy >= 210] -= 1.25
                DPT_flag = 3
                IAO = 161
            elif year == 2019:
                DPT_flag = 3
                IAO = 128
                z[doy < 97] = np.nan
            else:
                no_proc = 1

        elif site == "NUK_L":
            if year == 2007:
                DPT_flag = 1
            elif year == 2008:
                z[doy > 212] -= 7.1
                z[doy > 211] -= 5
                assign_nans(156)
                DPT_flag = 3
                IAO = 121
            elif year == 2009:
                assign_nans([90, 92, 95, 102, 105, 142, 154, 159, 169, 174])
                z[doy > 236] -= 0.5
                DPT_flag = 3
                IAO = 136
            elif year == 2010:
                assign_nans(207)
                z[doy > 206] -= 11.1
                z[doy > 264] = np.nan
                DPT_flag = 3
                IAO = 107
            elif year == 2011:
                assign_nans([210, 238, 239, 240])
                z[doy > 209] -= 1.15
                z[(doy >= 237) & (doy <= 239)] -= 0.4
                z[(doy >= 241) & (doy <= 242)] -= 0.2
                z[doy > 280] = np.nan
                DPT_flag = 3
                IAO = 156
                albedo_flag = 0
            elif year == 2012:
                assign_nans(241)
                z[doy > 241] -= 12
                DPT_flag = 3
                IAO = 150
            elif year == 2013:
                assign_nans([23, 24, 204])
                z[doy > 24] += 2.6
                z[doy > 204] -= 5.25
                DPT_flag = 3
                IAO = 155
                albedo_flag = 0
            elif year == 2014:
                z[doy > 206] -= 2.5
                z[doy > 207] -= 3.35
                DPT_flag = 3
                IAO = 121
                assign_nans([54, 55])
            elif year in [2015, 2016]:
                DPT_flag = 3
                if year == 2015:
                    IAO = 159
                else:
                    IAO = 100
            elif year == 2017:
                assign_nans([117, 118, 219, 220, 233, 234, 235])
                z[doy > 117] -= 14.92
                z[doy >= 142] -= 0.5
                DPT_flag = 3
                IAO = 142
            elif year == 2018:
                assign_nans([194, 196, 211])
                z[doy > 59] -= 14.75
                z[(doy >= 197) & (doy <= 210)] += 7.35
                z[(doy >= 194) & (doy <= 211)] += 0.27
                z[(doy >= 50) & (doy <= 70)] = np.nan
                DPT_flag = 2
            elif year == 2019:
                DPT_flag = 3
                IAO = 97
                albedo_flag = 0
            else:
                no_proc = 1

        elif site == "NUK_N":
            if year in [2010, 2013]:
                DPT_flag = 1
            elif year == 2011:
                z[doy >= 280] = np.nan
                DPT_flag = 3
                IAO = 164
            elif year == 2012:
                z[doy >= 240] = np.nan
                DPT_flag = 3
                IAO = 158
                z[(doy >= 90) & (doy <= 156) & (z < 0)] = np.nan
                assign_nans(96)
                z[doy == 168] -= 0.19
                assign_nans(np.arange(145, 152))
            elif year == 2014:
                DPT_flag = 3
                IAO = 165
            else:
                no_proc = 1

        elif site == "NUK_K":
            if year == 2014:
                z[(doy >= 225) & (doy <= 238)] += 1.19
                assign_nans(239)
                DPT_flag = 2
            elif year == 2015:
                DPT_flag = 3
                IAO = 200
            elif year == 2016:
                assign_nans([23, 24, 124])
                DPT_flag = 3
                IAO = 166
            elif year == 2017:
                DPT_flag = 3
                IAO = 204
            elif year in [2018, 2019]:
                DPT_flag = 3
                IAO = 176
                if year == 2019:
                    assign_nans([94, 127, 128, 129])
                elif year == 2018:
                    assign_nans(177)
                    z[doy > 177] += 0.07
            else:
                no_proc = 1

        elif site == "UPE_L":
            if year == 2009:
                z[(doy >= 120)] = np.nan
                DPT_flag = 1
            elif year == 2010:
                assign_nans(range(271, 275))
                DPT_flag = 3
                IAO = 129
            elif year == 2011:
                z[(doy >= 262) & (doy <= 311)] = np.nan
                DPT_flag = 3
                IAO = 157
            elif year == 2012:
                z[doy > 150] -= 0.3
                z[doy > 151] += 0.2
                z[doy > 152] += 0.1
                z[doy >= 225] -= 0.1
                z[doy >= 226] -= 0.15
                z[(doy >= 258) & (doy <= 308)] = np.nan
                DPT_flag = 3
                IAO = 148
                assign_nans([149, 150, 152, 153])
            elif year == 2013:
                z[doy <= 70] = np.nan
                assign_nans([74, 214, 215])
                z[(doy >= 73) & (doy <= 75)] = np.nan
                z[doy >= 214] -= 15.3
                DPT_flag = 3
                IAO = 160
            elif year == 2014:
                z[(doy >= 266) & (doy <= 270)] = np.nan
                DPT_flag = 3
                IAO = 158
            elif year == 2015:
                assign_nans(213)
                z[doy >= 214] -= 0.32
                z[(doy >= 253) & (doy <= 261)] = np.nan
                DPT_flag = 3
                IAO = 151
            elif year == 2016:
                DPT_flag = 0
            elif year == 2017:
                DPT_flag = 1
            elif year == 2018:
                DPT_flag = 3
                IAO = 167
            elif year == 2019:
                DPT_flag = 3
                IAO = 124
                z[doy < 99] = np.nan
            else:
                no_proc = 1

        elif site == "UPE_U":
            if year == 2009:
                DPT_flag = 1
            elif (year >= 2010) & (year <= 2012):
                DPT_flag = 3
                if year == 2010:
                    IAO = 143
                    assign_nans(114)
                elif year == 2011:
                    IAO = 159
                elif year == 2012:
                    IAO = 153
            elif year in [2013, 2018]:
                DPT_flag = 2
            elif year == 2014:
                DPT_flag = 3
                IAO = 171
            elif year == 2015:
                z[doy >= 216] -= 5.48
                assign_nans(215)
                DPT_flag = 2
            elif year == 2016:
                DPT_flag = 3
                IAO = 183
            elif year == 2017:
                if year == 2017:
                    IAO = 190
                    DPT_flag = 3
            elif year == 2019:
                DPT_flag = 3
                IAO = 161
            else:
                no_proc = 1

        elif site == "KPC_L":
            if year in [2008, 2010, 2012]:
                DPT_flag = 1
            elif year == 2009:
                DPT_flag = 0
            elif year == 2011:
                DPT_flag = 0
            elif year == 2013:
                DPT_flag = 3
                IAO = 163
            elif year == 2014:
                DPT_flag = 3
                IAO = 165
            elif year == 2015:
                DPT_flag = 3
                IAO = 172
            elif year == 2016:
                z[doy > 209] -= 6.297
                assign_nans([210, 211])
                DPT_flag = 3
                IAO = 166
                z[doy > 211] -= 0.1
            elif year == 2017:
                DPT_flag = 3
                IAO = 160
            elif year == 2018:
                DPT_flag = 2
            elif year == 2019:
                z[doy > 192] -= 4.478
                assign_nans(193)
                DPT_flag = 3
                IAO = 163
            else:
                no_proc = 1

        elif site == "KPC_U":
            if (year >= 2008) & (year <= 2019):
                DPT_flag = 1
            else:
                no_proc = 1

        elif site == "KAN_L":
            if year == 2008:
                DPT_flag = 1
            elif year == 2009:
                z[(doy >= 250) & (z >= -3.2)] = np.nan
                z[(doy >= 260) & (np.isnan(z))] = -3.5
                assign_nans(276)
                DPT_flag = 3
                IAO = 136
            elif year == 2010:
                z[(doy >= 260) & (z >= -5.1)] = np.nan
                z[doy >= 280] -= 0.1
                z[(doy >= 260) & (np.isnan(z))] = -5.45
                assign_nans([136, 284])
                DPT_flag = 3
                IAO = 121
            elif year == 2011:
                z[doy > 155] -= 4
                assign_nans(155)
                DPT_flag = 3
                IAO = 156
            elif year == 2012:
                z[doy >= 234] -= 11.56
                z[(doy >= 282) & (doy <= 298)] = np.nan
                assign_nans([234, 235])
                DPT_flag = 3
                IAO = 149
            elif year == 2013:
                DPT_flag = 3
                IAO = 141
            elif year == 2014:
                z[(doy >= 290) & (doy <= 319)] = np.nan
                DPT_flag = 3
                IAO = 146
                z[doy > 133] += 0.07
            elif year == 2015:
                z[doy >= 119] -= 7.08
                assign_nans([118, 278])
                z[doy >= 188] -= 0.29
                DPT_flag = 3
                IAO = 162
            elif year == 2016:
                z[(doy > 198)] -= 1.6
                z[(doy >= 260) & (z >= -4.5)] = np.nan
                z[(doy >= 260) & (z <= -5.5)] = np.nan
                z[(doy >= 273) & (doy <= 280)] = np.nan
                assign_nans(198)
                DPT_flag = 3
                IAO = 100
            elif year == 2017:
                z[doy > 244] -= 3.6
                z[doy > 58] -= z[doy == 59]
                assign_nans(244)
                DPT_flag = 3
                IAO = 123
            elif year == 2018:
                assign_nans(240)
                z[doy >= 241] -= 10.65
                DPT_flag = 3
                IAO = 137
            elif year == 2019:
                z[doy >= 248] -= 1.3
                DPT_flag = 3
                IAO = 115
            else:
                no_proc = 1

        elif site == "KAN_M":
            if (year >= 2008) & (year <= 2011):
                DPT_flag = 0
            elif year == 2012:
                DPT_flag = 3
                IAO = 169
            elif year == 2013:
                DPT_flag = 3
                IAO = 173
                assign_nans(183)
            elif year in [2014, 2015]:
                DPT_flag = 1
            elif year == 2016:
                DPT_flag = 3
                IAO = 156
            elif year == 2017:
                z[(doy > 100) & (doy < 135) & (z > 0.15)] = np.nan
                DPT_flag = 3
                IAO = 202
            elif year == 2018:
                z[doy >= 238] = np.nan
                DPT_flag = 2
            elif year == 2019:
                z[doy >= 247] = np.nan
                DPT_flag = 3
                IAO = 172
            else:
                no_proc = 1

        elif site == "QAS_M":
            if year == 2016:
                DPT_flag = 2
            elif year == 2017:
                assign_nans(range(141, 146))
                assign_nans([151, 152, 187, 188])
                z[(doy >= 196) & (doy <= 236)] = np.nan
                DPT_flag = 1
            elif year == 2018:
                z[doy > 242] -= 5.0
                assign_nans(range(241, 244))
                DPT_flag = 3
                IAO = 192
                z[doy >= 244] -= 0.2
            elif year == 2019:
                z[doy >= 242] -= 4.92
                DPT_flag = 3
                IAO = 178
            else:
                no_proc = 1

        elif site == "QAS_A":
            if (year >= 2012) & (year <= 2015):
                DPT_flag = 0
            else:
                no_proc = 1

        elif site == "QAS_L":
            if year in [2007, 2009, 2013]:
                DPT_flag = 1
            elif year == 2008:
                z[doy >= 219] -= 7
                z[doy >= 277] += 5.1
                assign_nans([77, 78, 92, 103, 110, 220, 277])
                z[doy >= 220] -= 1.9
                DPT_flag = 3
                IAO = 121
            elif year == 2010:
                z[(doy >= 127) & (doy <= 128)] = np.nan
                assign_nans(147)
                z[(doy >= 129) & (doy <= 146) & (doy <= 148)] -= 3.5
                z[(doy >= 144) & (doy <= 150)] = np.nan
                z[doy >= 129] -= 1.3
                DPT_flag = 2
                IAO = 129
            elif year == 2011:
                z[doy >= 222] = np.nan
                DPT_flag = 3
                IAO = 155
            elif year == 2012:
                assign_nans([133, 134])
                z[doy > 132] -= 0.55
                z[doy > 230] -= 0.8
                DPT_flag = 3
                IAO = 151
            elif year == 2014:
                assign_nans([125, 235])
                z[doy > 124] -= 1.48
                z[doy >= 236] -= 7.9
                z[doy > 159] -= 0.08
                DPT_flag = 3
                IAO = 149
            elif year == 2015:
                DPT_flag = 3
                IAO = 171
            elif year == 2016:
                assign_nans(224)
                z[doy >= 225] -= 13.3
                DPT_flag = 3
                IAO = 133
            elif year == 2017:
                DPT_flag = 3
                IAO = 152
                z[doy >= 152] -= 0.05
            elif year == 2018:
                assign_nans([123, 238])
                z[doy > 122] -= 1
                z[doy > 132] -= 0.9
                z[doy > 130] += 0.2
                z[doy > 126] -= 0.2
                z[doy >= 239] -= 8.2
                DPT_flag = 3
                IAO = 155
                assign_nans(np.arange(123, 133))
                assign_nans(np.arange(133, 137))
            elif year == 2019:
                z[doy > 142] = np.nan
                z[doy >= 177] -= 1
                z[doy >= 140] -= 0.4
                DPT_flag = 1
            else:
                no_proc = 1

        elif site == "QAS_U":
            if year in [2008, 2009, 2013]:
                DPT_flag = 1
            elif year == 2010:
                DPT_flag = 3
                IAO = 175
            elif year in [2011, 2015, 2018]:
                DPT_flag = 2
            elif year == 2012:
                DPT_flag = 3
                IAO = 196
                albedo_flag = 0
            elif year == 2014:
                DPT_flag = 3
                IAO = 206
            elif year == 2016:
                z[doy > 224] = np.nan
                DPT_flag = 3
                IAO = 202
            elif year == 2017:
                z[doy <= 143] = np.nan
                DPT_flag = 3
                IAO = 220
            elif year == 2019:
                z[doy >= 243] -= 2.24
                z[(doy >= 241) & (doy <= 242)] = np.nan
                z[doy <= 100] = np.nan
                DPT_flag = 3
                IAO = 200
            else:
                no_proc = 1

        elif site == "THU_L":
            if year in [2010, 2011, 2012, 2013, 2018]:
                DPT_flag = 1
            elif year == 2014:
                z[doy > 255] = np.nan
                DPT_flag = 3
                IAO = 177
            elif year == 2015:
                z[doy >= 20] += 0.1
                z[doy >= 182] -= 0.3
                z[doy >= 260] -= 0.2
                z[doy > 250] = np.nan
                DPT_flag = 3
                IAO = 163
            elif year == 2016:
                z[doy >= 204] -= 0.2
                DPT_flag = 3
                IAO = 180
            elif year == 2017:
                DPT_flag = 3
                IAO = 194
            elif year == 2019:
                z[(doy > 43) & (doy < 59)] = np.nan
                z[(doy < 100) & (z < -0.2)] = np.nan
                DPT_flag = 3
                IAO = 161
                albedo_flag = 0
            else:
                no_proc = 1

        elif site == "THU_U":
            if year in [2010, 2018]:
                DPT_flag = 1
            elif year in [2011, 2013, 2015, 2019]:
                DPT_flag = 2
            elif year == 2012:
                DPT_flag = 3
                IAO = 191
            elif year == 2014:
                z[doy >= 258] += 0.85
                assign_nans(258)
                DPT_flag = 2
            elif year == 2016:
                z[doy >= 204] += 2.5
                assign_nans(203)
                z[z < -0.5] = np.nan
                DPT_flag = 2
            elif year == 2017:
                z[z < -0.06] = np.nan
                DPT_flag = 1
            else:
                no_proc = 1

        elif site == "THU_U2":
            if year == 2018:
                assign_nans(216)
                DPT_flag = 2
            elif year == 2019:
                DPT_flag = 1  # instrument broken

        elif site == "SCO_L":
            if year == 2008:
                DPT_flag = 2
            elif (year >= 2009) & (year <= 2013):
                DPT_flag = 3
                if year == 2009:
                    IAO = 143
                elif year == 2010:
                    IAO = 144
                    albedo_flag = 0
                elif year == 2011:
                    IAO = 158
                elif year == 2012:
                    IAO = 150
                elif year == 2013:
                    IAO = 151
            elif year == 2014:
                z[doy >= 222] -= 14.79
                assign_nans(221)
                DPT_flag = 2
            elif year in [2015, 2016]:
                if year == 2015:
                    DPT_flag = 1
                else:
                    assign_nans(np.arange(115, 120))
                    DPT_flag = 3
                    IAO = 158
            elif year == 2017:
                assign_nans(np.arange(110, 116))
                z[doy >= 219] -= 9.88
                assign_nans(218)
                DPT_flag = 3
                IAO = 146
            elif year == 2018:
                DPT_flag = 3
                IAO = 148
            elif year == 2019:
                DPT_flag = 1
            else:
                no_proc = 1

        elif site == "SCO_U":
            if year in [2008, 2018]:
                DPT_flag = 2
            elif year == 2009:
                z[doy >= 316] = np.nan
                DPT_flag = 3
                IAO = 184
                albedo_flag = 0
            elif year == 2010:
                DPT_flag = 1
            elif year == 2011:
                DPT_flag = 3
                IAO = 161
            elif year == 2012:
                z[doy >= 242] -= 1.47
                assign_nans(241)
                DPT_flag = 3
                IAO = 150
            elif year in [2013, 2014]:
                DPT_flag = 3
                if year == 2013:
                    z[(doy > 22) & (doy < 26)] = np.nan
                    IAO = 155
                else:
                    IAO = 162
            elif year == 2015:
                assign_nans([101, 102, 103, 113])
                DPT_flag = 3
                IAO = 170
            elif year == 2016:
                assign_nans([47, 60, 61])
                DPT_flag = 3
                IAO = 160
            elif year == 2017:
                z[doy >= 217] -= 12.08
                assign_nans(216)
                DPT_flag = 3
                IAO = 153
                z[(doy > 49) & (doy < 54)] = np.nan
                assign_nans(np.arange(82, 93))
            elif year == 2019:
                DPT_flag = 3
                IAO = 160
            else:
                no_proc = 1

        elif site == "TAS_A":
            if year in [2013, 2015, 2018]:
                DPT_flag = 1
            elif year == 2014:
                z[(doy > 41) & (doy < 75)] = np.nan
                DPT_flag = 3
                IAO = 188
            elif year == 2016:
                DPT_flag = 3
                IAO = 176
            elif year == 2017:
                DPT_flag = 3
                IAO = 210
            elif year == 2019:
                DPT_flag = 3
                IAO = 195
                albedo_flag = 0
            else:
                no_proc = 1

        elif site == "TAS_L":
            if year in range(2007, 2010):
                DPT_flag = 0
            elif year in [2010, 2011, 2014, 2015, 2016]:
                DPT_flag = 1
            elif year == 2012:
                z[doy >= 252] -= 2.86
                assign_nans(251)
                z[doy <= 83] = np.nan
                DPT_flag = 3
                IAO = 153
            elif year == 2013:
                z[doy >= 235] -= np.nan
                DPT_flag = 2
                IAO = 140
            elif year == 2017:
                z[doy >= 209] -= 10.36
                assign_nans([207, 208])
                DPT_flag = 2
            elif year == 2018:
                z[doy >= 272] += 1.29
                assign_nans(271)
                DPT_flag = 3
                IAO = 156
            elif year == 2019:
                z[doy >= 121] -= 1.135
                DPT_flag = 3
                IAO = 143
                albedo_flag = 0
            else:
                no_proc = 1

        elif site == "TAS_U":
            if year == 2008:
                DPT_flag = 3
                IAO = 167
                assign_nans(143)
            elif year in [2009, 2012]:
                DPT_flag = 1
            elif year == 2010:
                z[doy >= 219] -= 8.07
                assign_nans(218)
                DPT_flag = 2
            elif year == 2011:
                z[doy >= 203] -= 1.558
                assign_nans([202, 222])
                DPT_flag = 3
                IAO = 161
                albedo_flag = 0
            elif year == 2013:
                DPT_flag = 3
                IAO = 167
            elif year == 2014:
                z[doy >= 216] -= 1.29
                assign_nans([138, 215])
                DPT_flag = 3
                IAO = 136
                z[doy >= 138] += 1.27
                assign_nans(137)
            elif year == 2015:
                DPT_flag = 2
            else:
                no_proc = 1

        elif site == "MIT":
            if year == 2009:
                z[doy >= 224] -= 6.6
                assign_nans([124, 223])
                DPT_flag = 3
                IAO = 222
                albedo_flag = 0
            elif year == 2010:
                DPT_flag = 3
                IAO = 190
            elif year == 2011:
                assign_nans(224)
                DPT_flag = 3
                IAO = 203
                z[doy > 223] -= 0.08
            elif year == 2012:
                z[doy >= 250] -= 13.70
                assign_nans(249)
                DPT_flag = 2
            elif year in [2013, 2015, 2016]:
                DPT_flag = 1
            elif year == 2014:
                DPT_flag = 3
                IAO = 188
            elif year == 2017:
                assign_nans(207)
                DPT_flag = 3
                IAO = 208
                albedo_flag = 0
            elif year == 2018:
                DPT_flag = 2
            elif year == 2019:
                z[doy < 12] += 4.4
                DPT_flag = 3
                IAO = 186
                albedo_flag = 0
            else:
                no_proc = 1

        # stations in the accumulation area are not used for this application
        elif site == "EGP" or site == "CEN" or site == "KAN_U":
            print("WARNING: No processing available for %s" % site)
            DPT_flag = 0

        if no_proc == 1:
            print("WARNING: No processing available for %s %s" % (site, year))
            DPT_flag = 0

        if DPT_flag != 3:
            IAO = np.nan

        albedo = df["Albedo_theta<70d"]

        # correction for measurement platform obstruction of the radiometer field
        # of view after Aoki et al (2011) that increases average PROMICE BBA by 0.034
        albedo += 0.034

        if albedo_flag == 0:
            albedo[:] = np.nan
        if np.sum(np.isnan(albedo)) == len(albedo):
            albedo_flag = 0

        # adjust DPT to have a null pre-melt season ice ablation
        z -= np.nanmean(z[(doy < IAO) & (doy > IAO - 45)])

        # determine the start of significant ice ablation (>6cm after IAO)
        if DPT_flag == 3 and albedo_flag == 1:
            try:
                indx = np.where(z[doy > IAO] < -0.06)[0][0]
                BID = doy[doy > IAO].iloc[indx]
            except IndexError:
                BID = np.nan
                albedo_flag = 0
        else:
            BID = np.nan

        return z, albedo, DPT_flag, BID, albedo_flag

    def plot_variables():

        fs = 13
        mpl.rc("xtick", labelsize=fs)
        mpl.rc("ytick", labelsize=fs)
        mpl.rc("lines", markersize=3)

        plt.figure(figsize=(10, 15))
        ax1 = plt.subplot(411)
        ln1 = ax1.plot(
            df_y_init["DayOfYear"],
            df_y_init["DepthPressureTransducer_Cor(m)"],
            "ro-",
            label="Raw",
            zorder=3,
        )
        ax1.set_ylabel("Ice ablation (meters)", fontsize=fs, color="b")
        ax2 = ax1.twinx()
        ln2 = ax2.plot(df_y["DayOfYear"], DPT_proc, "go-", label="Processed", zorder=2)
        ax2.set_ylabel("Ice ablation (meters)", fontsize=fs, color="g")
        ax2.axhline(0, color="gray", LineStyle="--", zorder=1)
        ax2.legend(loc="upper left")
        lns = ln1 + ln2
        labs = [l.get_label() for l in lns]
        plt.legend(lns, labs)
        ax1.get_xaxis().set_visible(False)

        ax2 = plt.subplot(412, sharex=ax1)
        ax2.plot(
            df_y_init["DayOfYear"], df_y_init["Albedo_theta<70d"], "o-", color="purple"
        )
        ax2.set_ylabel("Albedo (unitless)", fontsize=fs, color="purple")
        ax2.get_xaxis().set_visible(False)

        ax3 = plt.subplot(413, sharex=ax1)
        ax3.plot(df_y_init["DayOfYear"], HeightSensorBoom_m, "o-", color="orange")
        ax3.set_ylabel("Boom height (meters)", fontsize=fs, color="orange")
        ax3.get_xaxis().set_visible(False)

        ax4 = plt.subplot(414, sharex=ax1)
        ax4.plot(df_y_init["DayOfYear"], df_y_init["AirTemperature(C)"], "ro-")
        ax4.set_ylabel("Air temperature (°C)", fontsize=fs, color="red")
        ax4.set_xlabel("Day of year (DOY)", fontsize=fs)
        ax4.axhline(0, color="gray", LineStyle="--")
        plt.suptitle("%s %s" % (site, y), fontsize=fs + 6)

        if fig_save:
            plt.savefig(
                fig_path + site + "_" + str(y) + "_DPT_proc.png",
                dpi=300,
                bbox_inches="tight",
            )

    if yr == "all" or isinstance(yr, list):

        # search for available years
        if yr == "all":
            years = list(Counter(promice_data.Year))
        else:
            years = yr

        # run BIC_processing() for each year
        for i, y in enumerate(years):

            df_y = promice_data[promice_data.Year == int(y)]
            df_y_init = df_y.copy()

            DPT_proc, albedo, DPT_flag, BID, albedo_flag = DPT_processing(
                df=df_y, year=int(y)
            )
            DPT_flags = np.zeros(len(DPT_proc))
            DPT_flags[:] = DPT_flag
            BIDs = np.zeros(len(DPT_proc))
            BIDs[:] = BID
            albedo_flags = np.zeros(len(DPT_proc))
            albedo_flags[:] = albedo_flag

            HeightSensorBoom_m = df_y["HeightSensorBoom(m)"]

            # reallocate processed variables to outputs
            if i == 0:
                promice_data_proc = pd.DataFrame(
                    {
                        "DPT_proc": DPT_proc,
                        "DPT_flag": DPT_flags,
                        "BID": BIDs,
                        "DOY": df_y["DayOfYear"],
                        "Year": df_y["Year"],
                        "Air_temperature_C": df_y["AirTemperature(C)"],
                        "Albedo_theta_inf_70d": albedo,
                        "Albedo_flag": albedo_flags,
                        "HeightSensorBoom_m": HeightSensorBoom_m,
                    }
                )

            else:
                promice_data_proc_sc = pd.DataFrame(
                    {
                        "DPT_proc": DPT_proc,
                        "DPT_flag": DPT_flags,
                        "BID": BIDs,
                        "DOY": df_y["DayOfYear"],
                        "Year": df_y["Year"],
                        "Air_temperature_C": df_y["AirTemperature(C)"],
                        "Albedo_theta_inf_70d": albedo,
                        "Albedo_flag": albedo_flags,
                        "Height_sensor_boom_m": HeightSensorBoom_m,
                    }
                )

                promice_data_proc = promice_data_proc.append(promice_data_proc_sc)

            if visualisation:
                plot_variables()

    else:

        df_y = promice_data.copy()

        # run BIC_processing() for the selected year
        DPT_proc, albedo, DPT_flag, BID, albedo_flag = DPT_processing(
            df=df_y, year=int(y)
        )
        DPT_flags = np.zeros(len(DPT_proc))
        DPT_flags[:] = DPT_flag
        BIDs = np.zeros(len(DPT_proc))
        BIDs[:] = BID
        albedo_flags = np.zeros(len(DPT_proc))
        albedo_flags[:] = albedo_flag

        HeightSensorBoom_m = df_y["HeightSensorBoom(m)"]

        if y == 2015 and site == "THU_L":
            HeightSensorBoom_m[doy == 184] = np.nan
        elif y == 2012 and site == "THU_U":
            HeightSensorBoom_m[doy == 172] = np.nan
            HeightSensorBoom_m[doy == 211] = np.nan

        # reallocate processed variables to outputs
        promice_data_proc = pd.DataFrame(
            {
                "DPT_proc": DPT_proc,
                "DPT_flag": DPT_flags,
                "BID": BIDs,
                "DOY": df_y["DayOfYear"],
                "Year": df_y["Year"],
                "Air_temperature_C": df_y["AirTemperature(C)"],
                "Albedo_theta_inf_70d": albedo,
                "Albedo_flag": albedo_flags,
                "HeightSensorBoom_m": HeightSensorBoom_m,
            }
        )

        if visualisation:
            plot_variables()

    return promice_data_proc


def BIC_composite(inpath, dt=45, verbose=True, visualisation=True, save_pkl=False):

    """

    Composites of air temperature, snow height, ice ablation and albedo time series
    centered on bare ice appearance and spanning ± dt days.


    INPUTS:
        inpath: Path to the folder containing PROMICE .txt files [string]
        dt: Number of days to consider around bare ice appearance
            (default: 45) [int]
        verbose: If True, prints and counts station years (default: True) [boolean]
        visualisation: If True, displays composite (default: True) [boolean]
        save_pkl: If True, saves composite dictionnary into a .pkl file
                  (default: False) [boolean]

    OUTPUTS:
        composite: Dictionnary containing air_temperature_degrees, ice_ablation_meters,
                   snow_height_meters, albedo_unitless keys and associated dataframes
                   of station years (per row) spanning BID±dt (along columns). [dictionnary]

    """

    files = glob.glob(inpath + "*day_v03_upd.txt")
    sites = []

    for f in files:
        site = f.split(os.sep)[-1].split("_")[0] + "_" + f.split(os.sep)[-1].split("_")[1]
        if site[:3] == "MIT" or site[:3] == "CEN" or site[:3] == "EGP":
            site = site.split("_")[0]
        sites.append(site)

    sites = list(Counter(sites))

    variables = {
        "air_temperature_degrees": "Air_temperature_C",
        "ice_ablation_meters": "DPT_proc",
        "snow_height_meters": "Height_sensor_boom_m",
        "albedo_unitless": "Albedo_theta_inf_70d",
    }

    # compute timespan associated with ±dt
    time_span = 2 * dt + 1

    composite = {var: np.zeros(time_span) for var in variables}

    station_years = []
    stations = []

    for s in sites:

        # stations in the accumulation area are not used for this application
        if s == "CEN" or s == "EGP" or s == "KAN_U":
            continue
        else:
            df = load_data(file=inpath + s + "_day_v03_upd.txt", year="all")
            df_corr = BIC_processing(df)

            years = list(Counter(df_corr.Year))

            for year in years:
                df_yr = df_corr[df_corr.Year == int(year)]

                if df_yr.DPT_flag.iloc[0] == 3 and df_yr.Albedo_flag.iloc[0] == 1:
                    station_years.append(s)
                    stations.append(s + " " + str(year))
                    BID = int(df_yr.BID.iloc[0])
                    df_yr_comp = df_yr[
                        (df_yr.DOY >= BID - dt) & (df_yr.DOY <= BID + dt)
                    ]

                    if verbose:
                        print(s, year)

                    # append nan values if time_span higher than data set coverage
                    if len(df_yr_comp) != time_span:
                        to_append = np.zeros(
                            (time_span - len(df_yr_comp), np.shape(df_yr_comp)[1])
                        )
                        to_append[:] = np.nan
                        df_yr_comp = df_yr_comp.append(pd.DataFrame(to_append))

                    mean_boom_height = np.nanmean(
                        df_yr_comp.Height_sensor_boom_m[dt + 10 : dt + 30]
                    )
                    df_yr_comp.Height_sensor_boom_m = (
                        mean_boom_height - df_yr_comp.Height_sensor_boom_m
                    )

                    composite = {
                        key: np.vstack((composite[key], df_yr_comp[value]))
                        for key, value in variables.items()
                    }

    composite = {var: composite[var][1:, :] for var in variables}

    if verbose:
        nb_stations = len(list(Counter(station_years)))
        print("%s station years from %s stations" % (len(station_years), nb_stations))

    if visualisation:

        var_names = []
        [var_names.append(key) for key in composite.keys()]
        plt.style.use("default")
        plt.figure(figsize=(15, 15))
        op = 0

        for i, var in enumerate(composite.values()):

            if var_names[i] in [
                "ice_ablation_meters",
                "albedo_unitless",
                "air_temperature_degrees",
                "snow_height_meters",
            ]:
                op += 1
                fs = 20
                ax = plt.subplot(2, 2, op)
                plt.plot(np.nanmean(var, axis=0), "k-", zorder=1)
                plt.plot(np.nanmean(var, axis=0), "ko", alpha=0.3, markersize=5)
                plt.fill_between(
                    np.arange(0, len(var.T)),
                    np.nanmean(var, axis=0) - np.nanstd(var, axis=0),
                    np.nanmean(var, axis=0) + np.nanstd(var, axis=0),
                    color="gray",
                    alpha=0.2,
                )
                plt.xticks(
                    np.arange(0, time_span, 10),
                    np.arange(-dt, time_span - dt, 10),
                    fontsize=fs - 2,
                )
                plt.yticks(fontsize=fs - 2)
                ylims = ax.get_ylim()
                xlims = ax.get_xlim()
                plt.xlabel("bare ice day", fontsize=fs)

                if var_names[i] == "albedo_unitless":

                    plt.hlines(
                        np.nanmean(var, axis=0)[dt],
                        LineStyle="--",
                        color="red",
                        xmin=xlims[0],
                        xmax=dt,
                        label="%.3f±%.3f"
                        % (np.nanmean(var, axis=0)[dt], np.nanstd(var, axis=0)[dt]),
                        zorder=3,
                    )
                    plt.ylabel(
                        "%s, %s"
                        % (var_names[i].split("_")[0], var_names[i].split("_")[1]),
                        fontsize=fs,
                    )
                    plt.legend(fontsize=fs - 5, loc="upper right")
                    plt.vlines(
                        dt,
                        ymin=ylims[0],
                        ymax=np.nanmean(var, axis=0)[dt],
                        LineStyle="--",
                        color="black",
                        zorder=3,
                    )

                else:

                    plt.hlines(
                        np.nanmean(var, axis=0)[dt],
                        LineStyle="--",
                        color="black",
                        xmin=xlims[0],
                        xmax=dt,
                        zorder=3,
                    )
                    plt.axhline(0, color="gray", zorder=2, LineStyle="--")
                    plt.ylabel(
                        "%s, %s"
                        % (
                            var_names[i].split("_")[0]
                            + " "
                            + var_names[i].split("_")[1],
                            var_names[i].split("_")[2],
                        ),
                        fontsize=fs,
                    )
                    plt.vlines(
                        dt,
                        ymin=ylims[0],
                        ymax=np.nanmean(var, axis=0)[dt],
                        LineStyle="--",
                        color="black",
                        zorder=3,
                    )

                plt.ylim([ylims[0], ylims[1]])
                plt.xlim([0, 2 * dt])
                plt.subplots_adjust(hspace=0.25, wspace=0.25)

    # convert station years to index-friendly format
    stations = [st.replace(" ", "_") for st in stations]

    # load Dataframes into output dictionnary
    composite = {
        var: pd.DataFrame(
            composite[var], index=stations, columns=np.arange(-dt, time_span - dt)
        )
        for var in composite
    }

    if save_pkl:

        file_name = "PROMICE_composite_%dBID" % dt

        f = open(inpath + file_name + ".pkl", "wb")
        pickle.dump(composite, f)
        f.close()

    return composite
