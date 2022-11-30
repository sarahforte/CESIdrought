# CESIdrought
#######################################################################
#
#    -###- CANADIAN ENVIRONMENTAL SUSTAINABILITY INDICATORS -###-
#              -- Water Quantity Indicator Calculator --
#
#  This project contains scripts used to automate the calculation
#          of CESI drought/dry year indicator.
#
#  Environment and Climate Change Canada
#  Created: November 2022
#
#######################################################################

# These scripts calculate the presence/absence of hydrological drought
# during the summer at hydrometric station based on a 30 year reference
# period.
# Different methodologies are used for perennial and intermittent
# watercourses, but both involve defining a station-specific threshold
# and comparing the yearly summer data to the threshold.
# Trends over 50 years are calculated with the results and a kriged
# pan-Canadian map is created for the trends using regional kriging. Five
# grouped ecozones are used as regions because they can be expected to have
# similar drought responses.
# The specific steps are:
#    1) Define variables
#    2) Create tables to compile results
#    3) Retrieve hydrometric data from hydat database
#    4) Check if there is sufficient data for calculations
#    5) Identify summer period
#    6A) For perennial watercourses:
#       i. calculate rolling average
#      ii. define threshold by fitting curve to 30-year reference data
#     iii. evaluate each year against threshold
#    6B) For intermittent watercourses:
#       i. calculate dry spell durations
#      ii. define threshold from cumulative distribution of dry spell duration
#     iii. evaluate each year against threshold
#    7) Identify 50-year trends using 4 tests: hurdle, negative binomial, Wald-
#         Wolfowitz and all zeros
#    8) Save output as csv files
#    9) Krige results and create pan-Canadian map

# Assuming there are three sub-folders in working directory: 
# 1. "Dependencies" - with following files
#             station list [RHBN_U.csv]
#             country boundaries [CanadaBound.shp]
#             kriging zones [5 merged ecozones.shp]
#             outlines of big lakes [CANwaterbodies.shp]
#             CESI north arrow [esri6north.png]
# 2. "Output" - in which to save the output
# 3. "R" - in which 2 script files are located:
#             [CESI hydrological drought.R-THIS ONE]
#             [CESI hydrological drought-functions.R]
#######################################################################
