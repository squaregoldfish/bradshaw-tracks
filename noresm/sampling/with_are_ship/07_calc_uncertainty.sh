#!/bin/bash

R --slave --no-save track="$1" < calc_uncertainty.R
