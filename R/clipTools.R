createAttributesWithinGroup = function(h5, newFile, group="/") {
  for (attr in hdf5r::list.attributes(h5[[group]])) {
    hdf5r::h5attr(newFile[[group]], attr) = hdf5r::h5attr(h5[[group]], attr)
  }
}
