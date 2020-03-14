downloadDep = function(name, file, url, origName = "") {
  if(!file.exists(file.path(".",name,file))) {
    print(paste0("Downloading ",name,"..."))
    download.file(url, "lib.tar.gz", quiet=FALSE)
    untar("lib.tar.gz", exdir=".")
    unlink("lib.tar.gz")
    if (origName != "") {
      unlink(name, recursive=TRUE)
    	file.rename(origName, name)
    }
  }
}


downloadDepBitBucket = function(name, file, origName) {
  fileCheck = file.path(".",name,file)
  if(!file.exists(fileCheck)) {
    print(paste0("Downloading ",name,"..."))
    url = paste0("https://bitbucket.org/caiohamamura/",name,"/get/v0.4.0.zip")
    download.file(url, "lib.zip", quiet=FALSE)
    unzip("lib.zip", exdir=".")
    unlink("lib.zip")
    unlink(name, recursive=TRUE)
    file.rename(origName, name)
  }
}

downloadDep("cmpfit-1.2",
            "mpfit.h",
            "https://www.physics.wisc.edu/~craigm/idl/down/cmpfit-1.2.tar.gz")
downloadDep("libclidar",
	    "libLasProcess.h",
	    "https://github.com/caiohamamura/libclidar/archive/v0.4.0.tar.gz",
	    "libclidar-0.4.0")
downloadDepBitBucket("gedisimulator",
                     "gediRat.c",
                     "caiohamamura-gedisimulator-63539ca0598c")
downloadDep("tools",
            "tools.c",
            "https://github.com/caiohamamura/tools/archive/v0.4.0.tar.gz",
	    "tools-0.4.0")
