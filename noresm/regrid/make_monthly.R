# Make a monthly data set from a daily one
library(ncdf4)

MONTH_START_DAYS <- c(0, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366)
MONTH_START_PORTION <- c(0.0, 0.088, 0.164, 0.249, 0.332, 0.416, 0.499, 0.584, 0.668, 0.751, 0.836, 0.918) 

makeMonthly <- function(series) {

    series_month <- vector(mode="numeric", length=0)
    monthly <- vector(mode="numeric", length=(length(series) / 365 * 12))

    index <- 0
    month_index <- 0
    current_month <- 1

    while (index < length(series)) {
        month_finished <- FALSE
        year_finished <- FALSE

        index <- index + 1

        day_of_year <- index %% 365
        if (day_of_year == 0) {
            day_of_year <- 365
        }

        if (day_of_year == 1 && current_month != 1) {
            month_finished <- TRUE
            year_finished <- TRUE
        } else if (day_of_year >= MONTH_START_DAYS[current_month + 1]) {
            month_finished <- TRUE
        }

        if (month_finished) {
            month_index <- month_index + 1

            series_mean <- NA
            if (length(series_month) > 0) {
                series_mean <- mean(series_month)
            }
            monthly[month_index] <- series_mean

            current_month <- current_month + 1
            if (year_finished) {
                current_month <- 1
                year_finished <- FALSE
            }

            month_finished <- FALSE

            series_month <- vector(mode="numeric", length=0)
        }

        if (!is.na(series[index])) {
            series_month[length(series_month) + 1] <- series[index]
        }

    }

    month_index <- month_index + 1
    series_mean <- NA
    if (length(series_month) > 0) {
        series_mean <- mean(series_month)
    }
    monthly[month_index] <- series_mean

    return (monthly)
}

#============================================================

# Command line parameters
for (arg in commandArgs()) {
    argset <- strsplit(arg, "=", fixed=TRUE)
    if (!is.na(argset[[1]][2])) {
        if (argset[[1]][1] == "infile") {
            assign("infile",argset[[1]][2])
        } else if (argset[[1]][1] == "outfile") {
            assign("outfile",argset[[1]][2])
        } else if (argset[[1]][1] == "start_year") {
            assign("start_year",argset[[1]][2])
        }
    }
}

start_year <- as.integer(start_year)

nc <- nc_open(infile)
lons <- ncvar_get(nc, "lon")
lats <- ncvar_get(nc, "lat")
pco2 <- ncvar_get(nc, "pco2")
nc_close(nc)

years <- 31

out_times = vector(mode="numeric", length = years * 12)
current_year <- start_year
current_month <- 1
for (y in 1:length(out_times)) {
	out_times[y] <- current_year + MONTH_START_PORTION[current_month]
	current_month <- current_month + 1
	if (current_month > 12) {
		current_year <- current_year + 1
		current_month <- 1
	}
}

out_pco2 <- vector(mode="numeric", length=(144 * 72 * length(out_times)))
dim(out_pco2) <- c(144, 72, length(out_times))

for (i in 1:144) {
	cat("\r",i)
	for (j in 1:72) {
		# Truncate to match output length
		out_pco2[i, j, ] <- makeMonthly(pco2[i, j, ])[1:length(out_times)]
	}
}
cat("\n")
lon_dim <- ncdim_def("lon", "degrees_east", lons)
lat_dim <- ncdim_def("lat", "degrees_north", lats)
time_dim <- ncdim_def("time", "year", out_times, unlim=TRUE)

pco2_var <- ncvar_def("pco2", "ppm", list(lon_dim, lat_dim, time_dim), -1e35, prec="double")

nc <- nc_create(outfile, list(pco2_var))
ncvar_put(nc, pco2_var, out_pco2)
nc_close(nc)
