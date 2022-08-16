#!/usr/bkin/env python3.6

# helper functions for genomics analysis
# Ben Kinnersley (b.kinnersley@ucl.ac.uk)

import datetime
import re

def convert_date_to_datetime(x):
  x_split = re.findall(r'[^\/]+', x)
  
  # check for cases where year is incorrect (e.g. 2018 rather than 18)
  if len(str(x_split[2])) == 4:
    year = int(''.join(list(x_split[2])[2:]))
  else:
    year = int(x_split[2])
    month = int(x_split[1])
    day = int(x_split[0])
   
  x_datetime = datetime.date(year, month, day)
   
  return x_datetime
   
def convert_timedelta_to_months(x):
  x_days = x.days
  x_months = x_days / 30
  return x_months
