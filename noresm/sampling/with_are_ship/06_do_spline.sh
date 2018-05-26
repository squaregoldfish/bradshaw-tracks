#!/bin/bash

R --slave --no-save track="$1" < do_spline.R
