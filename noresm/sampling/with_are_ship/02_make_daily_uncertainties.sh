#!/bin/bash

R --slave --no-save track="$1" < make_daily_uncertainties.R
