library(rGEDI)

las.sample=system.file("extdata/sample.las", package="rGEDI")

output = tempfile(fileext=".h5")

out.WF=gediWFSimulator(
  input=las.sample,
  output=output,
  coords =  c(278215, 602215)
)

output2 = tempfile(fileext="_noised.h5")

noised = gediWFNoise(
  input = output,
  output = output2,

)
noised$close()

outputMetric = tempfile(fileext="")
out.Metric2 = gediWFMetric(
  input = output2,
  outRoot = "E:/Documentos/sample_metric",
  deconTol = 1e-14
)
out.WF$close()
