#include <Rcpp.h>
#include <gdal_priv.h>

using namespace Rcpp;

GDALDataset *create_dataset(const char *output, int nbands, int datatype, const char *projection, double lat_min, double lat_max, double lon_min, double lon_max, std::vector<double> res, double nodata, CharacterVector co)
{
  CPLErr err = CE_None;
  int width = (int)ceil((lon_max - lon_min) / res[0]);
  int height = (int)ceil((lat_min - lat_max) / res[1]);

  GDALAllRegister();
  GDALDriverManager *driverMan = GetGDALDriverManager();
  GDALDriver *tiffDriver = driverMan->GetDriverByName("GTiff");
  if (tiffDriver == NULL)
    Rcpp::stop("Could not retrieve the driver!");

  std::vector<char *> charVec{};

  for (auto &option : co)
  {
    charVec.push_back(option);
  }
  charVec.push_back(nullptr);

  // Rcout << width << " " << height << "\n";
  GDALDataset *ds = tiffDriver->Create(output, width, height, nbands, (GDALDataType)datatype, charVec.data());

  if (ds == NULL)
    Rcpp::stop("Could not create file!");

  double transform[6] = {lon_min, res[0], 0, lat_max, 0, res[1]};
  ds->SetGeoTransform(transform);
  for (int i = 1; i <= nbands; i++)
  {
    GDALRasterBand *band = ds->GetRasterBand(i);
    err = band->SetNoDataValue(nodata);
    if (err == CE_Failure)
      Rcpp::stop(CPLGetLastErrorMsg());
  }
  err = ds->SetProjection(projection);

  if (err == CE_Failure)
    Rcpp::stop(CPLGetLastErrorMsg());

  return ds;
}

template <typename T, typename S>
S ReadBlock(GDALRasterBand *band, int iXBlock, int iYBlock)
{
  CPLErr res = CE_None;
  S output = NULL;

  int nXBlockSize, nYBlockSize;
  band->GetBlockSize(&nXBlockSize, &nYBlockSize);
  if (std::is_same<T, GInt32>::value || std::is_same<T, double>::value)
  {
    S vec(nXBlockSize * nYBlockSize);
    res = band->ReadBlock(iXBlock, iYBlock, vec.begin());
    output = vec;
  }
  else
  {
    std::vector<T> buffer(nXBlockSize * nYBlockSize);
    res = band->ReadBlock(iXBlock, iYBlock, buffer.data());
    output = Rcpp::wrap(buffer.begin(), buffer.end());
  }

  if (res == CE_Failure)
    Rcpp::stop(CPLGetLastErrorMsg());

  return output;
}

template <typename T>
void WriteBlock(GDALRasterBand *band, int iXBlock, int iYBlock, T buffer)
{
  GDALDataType dtype = band->GetRasterDataType();
  CPLErr res = CE_None;

  switch (dtype)
  {
  case GDALDataType::GDT_UInt16:
  {
    std::vector<GUInt16> vec(buffer.begin(), buffer.end());
    res = band->WriteBlock(iXBlock, iYBlock, vec.data());
    break;
  }

  case GDALDataType::GDT_Int16:
  {
    std::vector<GInt16> vec(buffer.begin(), buffer.end());
    res = band->WriteBlock(iXBlock, iYBlock, vec.data());
    break;
  }

  case GDALDataType::GDT_UInt32:
  {
    std::vector<GUInt32> vec(buffer.begin(), buffer.end());
    res = band->WriteBlock(iXBlock, iYBlock, vec.data());
    break;
  }

  case GDALDataType::GDT_Float32:
  {
    std::vector<float> vec(buffer.begin(), buffer.end());
    res = band->WriteBlock(iXBlock, iYBlock, vec.data());
    break;
  }

  default:
  {
    res = band->WriteBlock(iXBlock, iYBlock, buffer.begin());
    break;
  }
  }

  if (res == CE_Failure)
    Rcpp::stop(CPLGetLastErrorMsg());

  res = band->FlushBlock(iXBlock, iYBlock, 1);

  if (res == CE_Failure)
    Rcpp::stop(CPLGetLastErrorMsg());
}

IntegerVector GetBlockXSize(GDALRasterBand *band)
{
  IntegerVector xsize(1);
  band->GetBlockSize(xsize.begin(), NULL);
  return xsize;
}

IntegerVector GetBlockYSize(GDALRasterBand *band)
{
  IntegerVector ysize(1);
  band->GetBlockSize(NULL, ysize.begin());
  return ysize;
}

void RGDALClose(GDALDataset *ds)
{
  Rcout << "Closing ds\n";
  GDALDatasetH handle = GDALDataset::ToHandle(ds);
  GDALClose(handle);
  Rcout << "Closed!\n";
}

RCPP_MODULE(gdal_module)
{

  class_<GDALDataset>("CPP_GDALDataset")
      .method("GetRasterBand", &GDALDataset::GetRasterBand)
      .method("GetRasterXSize", &GDALDataset::GetRasterXSize)
      .method("GetRasterYSize", &GDALDataset::GetRasterYSize)
      .method("Close", &RGDALClose);

  class_<GDALRasterBand>("CPP_GDALRasterBand")
      .method("GetXSize", &GDALRasterBand::GetXSize)
      .method("ReadBlock1", &ReadBlock<GByte, RawVector>)
      .method("ReadBlock2", &ReadBlock<GUInt16, IntegerVector>)
      .method("ReadBlock3", &ReadBlock<GInt16, IntegerVector>)
      .method("ReadBlock4", &ReadBlock<GUInt32, IntegerVector>)
      .method("ReadBlock5", &ReadBlock<GInt32, IntegerVector>)
      .method("ReadBlock6", &ReadBlock<float, NumericVector>)
      .method("ReadBlock7", &ReadBlock<double, NumericVector>)
      .method("WriteBlock1", &WriteBlock<RawVector>)
      .method("WriteBlock2", &WriteBlock<IntegerVector>)
      .method("WriteBlock3", &WriteBlock<IntegerVector>)
      .method("WriteBlock4", &WriteBlock<IntegerVector>)
      .method("WriteBlock5", &WriteBlock<IntegerVector>)
      .method("WriteBlock6", &WriteBlock<NumericVector>)
      .method("WriteBlock7", &WriteBlock<NumericVector>)
      .method("GetBlockXSize", &GetBlockXSize)
      .method("GetBlockYSize", &GetBlockYSize);

  function("create_dataset", &create_dataset);
}
