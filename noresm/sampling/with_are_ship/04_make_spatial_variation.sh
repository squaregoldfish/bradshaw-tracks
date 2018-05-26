#!/bin/bash

R --slave --no-save track="$1" < make_spatial_variation.R
