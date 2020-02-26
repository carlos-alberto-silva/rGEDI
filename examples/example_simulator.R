library(rGEDI)

las.sample=system.file("extdata/Cerrado.las", package="rGEDI")
lasfile = lidR::readLAS(las.sample)

x = mean(lasfile@bbox[1,])
y = mean(lasfile@bbox[2,])
output = tempfile(fileext=".h5")

out.WF=gediWFSimulator(
  input=las.sample,
  output=output,
  coords = c(x, y)
)
out.WF$ls()[,1]
z = seq(out.WF[["Z0"]][], out.WF[["ZN"]][], length.out = out.WF[["NBINS"]][])
plot(z, out.WF[["RXWAVEINT"]][,])
summary(lasfile@data$Z)

output2 = tempfile(fileext="_noised.h5")

noised = gediWFNoise(
  input = output,
  output = output2
)
z = seq(noised[["Z0"]][], noised[["ZN"]][], length.out = noised[["NBINS"]][])
plot(z, noised[["RXWAVEINT"]][,])

outputMetric = tempfile(fileext="")
out.Metric = gediWFMetrics(
  input = output,
  outRoot = outputMetric
)

out.WF$ls()[,1]
noised$ls()[,1]
noised[["ZGDEM"]] = out.WF[["ZGDEM"]][]
noised$flush()

out.WF[["RXWAVEINT"]][,] = noised[["RXWAVEINT"]][,]
out.WF[["RXWAVEINT"]]$flush()


i = 0
i = i + 1
dt = out.WF$ls()[i,1]

dt = "RXWAVECOUNT"
print(dt)
noised[[dt]][]
out.WF[[dt]][]
noised[[dt]][406:486,]-out.WF[[dt]][406:486,]
plot(z, noised[[dt]][,]/), col="red")
points(z, out.WF[[dt]][,])


noised[["RXWAVECOUNT"]][,] = (noised[["RXWAVECOUNT"]][,]-min(noised[["RXWAVECOUNT"]][,]))/abs(diff(range(noised[["RXWAVECOUNT"]][,])))
noised[["RXWAVECOUNT"]][,] = noised[["RXWAVECOUNT"]][,]*0.25
noised[["RXWAVECOUNT"]]
noised$flush()

out.Metric2 = gediWFMetrics(
  input = output,
  outRoot = outputMetric
)
names(out.Metric2)
out.Metric2 - out.Metric
