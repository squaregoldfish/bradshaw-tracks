library(ncdf4)

MODEL_FILE <- as.vector(read.table("model_file.txt")[[1]])
SAMPLE_FILE <- as.vector(read.table("sample_file.txt")[[1]])
OUTPUT_ROOT <- as.vector(read.table("output_root.txt")[[1]])

cat("Reading model file...")
nc <- nc_open(MODEL_FILE)
model_pco2 <- ncvar_get(nc, "pco2")
nc_close(nc)

cat("\r\033[KReading sample file...")
nc <- nc_open(SAMPLE_FILE)
lons <- ncvar_get(nc, "LON")
lats <- ncvar_get(nc, "LAT")
times <- ncvar_get(nc, "TIME")
sample_pco2 <- ncvar_get(nc, "pco2")
nc_close(nc)

cat("\r\033[KSampling model output...")
model_pco2[is.na(sample_pco2)] <- NA

cat("\r\033[KWriting sampled model output...")

lon_dim <- ncdim_def("LON", "degrees_east", lons)
lat_dim <- ncdim_def("LAT", "degrees_north", lats)
time_dim <- ncdim_def("TIME", "year", times)

pco2_var <- ncvar_def("pco2", "uatm", list(lon_dim, lat_dim, time_dim), -999, prec="double")

nc <- nc_create(paste(OUTPUT_ROOT, "/daily.nc", sep=""), list(pco2_var))
ncvar_put(nc, "pco2", model_pco2)
nc_close(nc)

for (lon_loop in 1:144) {
    for (lat_loop in 1:72) {
        cat("\r\033[KWriting series file ", lon_loop, " ", lat_loop, sep="")
        sink(paste(OUTPUT_ROOT,"/cell_series_daily/cell_series_", lon_loop, "_", lat_loop, ".csv", sep=""))
        for (time_loop in 1:length(times)) {
            csv_value <- model_pco2[lon_loop, lat_loop, time_loop]
            if (csv_value == -999) {
                csv_value <- NA
            }
            cat(time_loop, ",", csv_value, "\n", sep="")
        }
        sink()
    }
}

cat("\n")
