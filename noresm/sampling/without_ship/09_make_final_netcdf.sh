#!/bin/bash

R --slave --no-save track="$1" < make_final_netcdf.R
