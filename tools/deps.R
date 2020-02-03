downloadDep = function(name, file, url) {
  if(!file.exists(paste0("./",name,"/",file))) {
    print(paste0("Downloading ",name,"..."))
    download.file(url, "lib.tar.gz", quiet=FALSE)
    untar("lib.tar.gz", exdir=".")
    unlink("lib.tar.gz")
  }
}


downloadDepBitBucket = function(name, file, origName) {
  if(!file.exists(paste0("./",name,"/",file))) {
    print(paste0("Downloading ",name,"..."))
    url = paste0("https://bitbucket.org/caiohamamura/",name,"/get/v0.1.2.zip")
    download.file(url, "lib.zip", quiet=FALSE)
    unzip("lib.zip", exdir=".")
    unlink("lib.zip")
    file.rename(origName, name)
  }
}

downloadDep("cmpfit-1.2",
            "mpfit.h",
            "https://www.physics.wisc.edu/~craigm/idl/down/cmpfit-1.2.tar.gz")
downloadDepBitBucket("gedisimulator",
                     "gediRat.c",
                     "caiohamamura-gedisimulator-daf127f577c4")
downloadDepBitBucket("tools",
                     "tools.c",
                     "caiohamamura-tools-51f4756ba29f")
downloadDepBitBucket("libclidar",
                     "libLasProcess.h",
                     "caiohamamura-libclidar-bf4cf5fe087a")
