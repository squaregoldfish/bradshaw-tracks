# Add uncertainty time series files alongside the actual time series

# Command line parameters
for (arg in commandArgs()) {
    argset <- strsplit(arg, "=", fixed=TRUE)
    if (!is.na(argset[[1]][2])) {
        if (argset[[1]][1] == "track") {
            assign("track",argset[[1]][2])
        }
    }
}

OUTPUT_ROOT <- paste(as.vector(read.table("output_root.txt")[[1]]), "/", track, "/without_ship", sep="")
SERIES_DIR <- paste(OUTPUT_ROOT, "/cell_series_daily", sep="")
SERIES_LENGTH <- 11680

for (lon in 1:144) {
    for (lat in 1:72) {
        cat("\r",lon,lat,"   ")

        series_file <- paste(SERIES_DIR, "/cell_series_", lon, "_", lat, ".csv", sep="")
        series <- as.double(unlist(read.csv(series_file,header=F)[2]))

        uncertainties <- array(NA, c(SERIES_LENGTH))
        uncertainties[!is.na(series)] <- 2.5
   
        weights <- array(NA, c(SERIES_LENGTH))
        weights[!is.na(series)] <- 1.0

        uncertainty_file <- paste(SERIES_DIR, "/cell_uncertainties_", lon, "_", lat, ".csv", sep="")
        sink(uncertainty_file)
        for (i in 1:SERIES_LENGTH) {
            cat(i, ",", uncertainties[i], "\n", sep="")
        }
        sink()

        weights_file <- paste(SERIES_DIR, "/cell_weights_", lon, "_", lat, ".csv", sep="")
        sink(weights_file)
        for (i in 1:SERIES_LENGTH) {
            cat(i, ",", weights[i], "\n", sep="")
        }
        sink()
    }
}

cat("\n")
