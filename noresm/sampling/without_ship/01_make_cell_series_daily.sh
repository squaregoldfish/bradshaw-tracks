#!/bin/bash

R --slave --no-save track="$1" < make_cell_series_daily.R
