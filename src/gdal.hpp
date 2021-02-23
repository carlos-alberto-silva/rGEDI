#include <Rcpp.h>
#include <gdal_priv.h>

using namespace Rcpp;

class GDALRasterBandR
{
private:
  GDALRasterBand *band;

public:
  GDALRasterBandR(GDALRasterBand *the_band)
  {
    band = the_band;
  }

  IntegerVector GetBlockXSize()
  {
    IntegerVector result(1);
    band->GetBlockSize(result.begin(), NULL);
    return result;
  }

  IntegerVector GetBlockYSize()
  {
    IntegerVector result(1);
    band->GetBlockSize(NULL, result.begin());
    return result;
  }

  template <typename T, typename S>
  S ReadBlock(int iXBlock, int iYBlock)
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
  void WriteBlock(int iXBlock, int iYBlock, T buffer)
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
};

class GDALDatasetR
{
private:
  GDALDataset *ds = NULL;

public:
  GDALDatasetR(GDALDataset* _ds) {
    ds = _ds;
  }
  
  GDALRasterBandR* GetRasterBand(int nband)
  {
    GDALRasterBandR* band = new GDALRasterBandR(ds->GetRasterBand(nband));
    return band;
  }

  int GetRasterXSize()
  {
    return ds->GetRasterXSize();
  }

  int GetRasterYSize()
  {
    return ds->GetRasterYSize();
  }

  void Close()
  {
    GDALClose(GDALDataset::ToHandle(ds));
  }
};

GDALDatasetR* create_dataset(const char *output, int nbands, int datatype, const char *projection, double lat_min, double lat_max, double lon_min, double lon_max, std::vector<double> res, double nodata, CharacterVector co)
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

  GDALDatasetR* outDs = new GDALDatasetR(ds);

  return outDs;
}

RCPP_MODULE(gdal_module)
{

  class_<GDALDatasetR>("CPP_GDALDataset")
      .method("GetRasterBand", &GDALDatasetR::GetRasterBand)
      .method("GetRasterXSize", &GDALDatasetR::GetRasterXSize)
      .method("GetRasterYSize", &GDALDatasetR::GetRasterYSize)
      .method("Close", &GDALDatasetR::Close);

  class_<GDALRasterBandR>("CPP_GDALRasterBand")
      .method("ReadBlock1", &GDALRasterBandR::ReadBlock<GByte, RawVector>)
      .method("ReadBlock2", &GDALRasterBandR::ReadBlock<GUInt16, IntegerVector>)
      .method("ReadBlock3", &GDALRasterBandR::ReadBlock<GInt16, IntegerVector>)
      .method("ReadBlock4", &GDALRasterBandR::ReadBlock<GUInt32, IntegerVector>)
      .method("ReadBlock5", &GDALRasterBandR::ReadBlock<GInt32, IntegerVector>)
      .method("ReadBlock6", &GDALRasterBandR::ReadBlock<float, NumericVector>)
      .method("ReadBlock7", &GDALRasterBandR::ReadBlock<double, NumericVector>)
      .method("WriteBlock1", &GDALRasterBandR::WriteBlock<RawVector>)
      .method("WriteBlock2", &GDALRasterBandR::WriteBlock<IntegerVector>)
      .method("WriteBlock3", &GDALRasterBandR::WriteBlock<IntegerVector>)
      .method("WriteBlock4", &GDALRasterBandR::WriteBlock<IntegerVector>)
      .method("WriteBlock5", &GDALRasterBandR::WriteBlock<IntegerVector>)
      .method("WriteBlock6", &GDALRasterBandR::WriteBlock<NumericVector>)
      .method("WriteBlock7", &GDALRasterBandR::WriteBlock<NumericVector>)
      .method("GetBlockXSize", &GDALRasterBandR::GetBlockXSize)
      .method("GetBlockYSize", &GDALRasterBandR::GetBlockYSize);

  function("create_dataset", &create_dataset);
}
