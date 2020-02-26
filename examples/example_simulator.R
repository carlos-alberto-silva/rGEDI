library(rGEDI)

las.sample=system.file("extdata/Amazon.las", package="rGEDI")
lasfile = lidR::readLAS(las.sample)

output = tempfile(fileext=".h5")

out.WF=gediWFSimulator(
  input=las.sample,
  output=output
)

output2 = tempfile(fileext="_noised.h5")

noised = gediWFNoise(
  input = output,
  output = output2,
)


outputMetric = tempfile(fileext="")
out.Metric2 = gediWFMetric(
  input = output2,
  outRoot = "E:/Documentos/sample_metric",
  deconTol = 1e-14
)
