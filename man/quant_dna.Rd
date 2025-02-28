% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quant_dna.R
\name{quant_dna}
\alias{quant_dna}
\title{Analyze quant fluorescence data}
\usage{
quant_dna(quant_data, trained_model)
}
\arguments{
\item{quant_data}{list. Output of \code{\link[=load_quant_files]{load_quant_files()}}, containing both the raw data and plate map information.}

\item{trained_model}{list. Output of \code{\link[=train_quant_standards]{train_quant_standards()}}, containing the trained regression model, background RFU, and a base ggplot object.}
}
\value{
Returns a list containing:
\itemize{
\item \code{quant_output}: data frame. Contains the sample data with predicted DNA concentrations (ng/well and ng/uL).
\item \code{plot}: ggplot object. Updated visualization of the regression model with sample predictions overlaid as jittered points.
}
}
\description{
Use a trained regression model to predict DNA concentrations for sample data in a fluorescent quant assay. This function processes the sample data from the quant assay, calculates zeroed RFUs, predicts DNA concentration (ng/well) based on the regression model, and normalizes it based on the sample volume (ng/uL). This takes the object generated via \code{\link[=train_quant_standards]{train_quant_standards()}} as input.

This is Step 3 of 4 in a series of functions for analyzing fluorescent data.
}
\details{
The function uses the background RFU from the trained model to calculate zeroed RFUs for the sample data. Predictions are based on the regression model from \code{trained_model}. For power models, predictions are transformed back from log space.
}
\examples{
require (tamuccGCL)

# Import data files
raw_data <- system.file("extdata", "raw_data.csv", package = "tamuccGCL")
plate_map <- system.file("extdata", "plate_map.csv", package = "tamuccGCL")

quant_data <- load_quant_files(raw_data, plate_map)
trained_model <- train_quant_standards (quant_data)
quant_report <- quant_dna(quant_data = quant_data, trained_model = trained_model)

}
\seealso{
\code{\link[=load_quant_files]{load_quant_files()}}, \code{\link[=train_quant_standards]{train_quant_standards()}}
}
