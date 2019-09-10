#include "itkImage.h"
#include "itkSimpleFuzzyConnectednessScalarImageFilter.h"
#include "itkConfidenceConnectedImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkCurvatureFlowImageFilter.h"
#include "itkImageFileReader.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkIndex.h"
#include "itkImageFileWriter.h"
#include "mat.h"
#include <iostream>
#include "itkImageToVTKImageFilter.h"
#include "vtkAutoInit.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkBinaryMask3DMeshSource.h"
#include "itkDeformableMesh3DFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkMesh.h"
#include "itkCovariantVector.h"
#include "itkPointSetToImageFilter.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkImageReader.h"
#include "vtkVolume16Reader.h"
#include "vtkImageCast.h"
#include "vtkMarchingCubes.h"
#include "vtkDecimatePro.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkPolyDataNormals.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h" 
#include "vtkOutlineFilter.h" 
#include "vtkCamera.h" 
#include "vtkProperty.h" 
#include "vtkContourFilter.h"
#include "vtkStripper.h"
#include "vtkAxesActor.h"
#include "vtkOrientationMarkerWidget.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkComposeImageFilter.h"
#include "itkScalarToArrayCastImageFilter.h"
#include "itkMRFImageFilter.h"
#include "itkDistanceToCentroidMembershipFunction.h"
#include "itkMinimumDecisionRule.h"
#include "itkScalarImageKmeansImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkExtractImageFilter.h"
#include "vtkImageViewer.h"
#include "vtkRenderWindowInteractor.h"
#include "itkMetaImageIO.h"
#include "itkMetaImageIOFactory.h"
#include "itkAddImageFilter.h"

#include "itkPermuteAxesImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkCurvatureFlowImageFilter.h"
#include "vtkAutoInit.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"


#include "itkVoronoiSegmentationImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkTIFFImageIOFactory.h"
#include "itkImageSeriesWriter.h"
#include "itkNumericSeriesFileNames.h"

using namespace std;


int main(int argc, char * argv[])
{
	
	itk::TIFFImageIOFactory::RegisterOneFactory();

	typedef itk::Image<double, 3> ImageType;

	MATFile *pmatFile = NULL;
	mxArray *pMxArray = NULL;
	double *S1;
	_int64 M, N, P;
	pmatFile = matOpen("WWZResult3.mat", "r");//Hessian matrix from matlab
	if (pmatFile == NULL) {
		cout << "MatOpen error!!!" << endl;
	}
	pMxArray = matGetVariable(pmatFile, "t");


	S1 = (double*)mxGetData(pMxArray);

	const mwSize *distrSize = mxGetDimensions(pMxArray);
	M = distrSize[0];
	N = distrSize[1];
	P = distrSize[2];

	
	//create 3D image
	ImageType::Pointer image = ImageType::New();
	ImageType::IndexType start;
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	ImageType::SizeType size;
	size[0] = M;
	size[1] = N;
	size[2] = P;
	ImageType::RegionType region;


	region.SetSize(size);
	region.SetIndex(start);
	image->SetRegions(region);
	image->Allocate();
	//Accessing voxels
	for (int k = 0;k < size[2];k++)
	{
		for (int j = 0;j < size[1];j++)
		{
			for (int i = 0;i < size[0];i++)
			{
				ImageType::IndexType point_temp;
				point_temp[0] = i;//x
				point_temp[1] = j;//y
				point_temp[2] = k;//z
				image->SetPixel(point_temp, S1[k*M*N + j*N + i]);
			}
		}
	}
	image->Update();



	/*FuzzyConnected,rough Segmentation*/
	
	/*typedef unsigned char LabelPixelType;
	typedef itk::Image<LabelPixelType, 3>  LabelImageType;*/

	typedef itk::Image<unsigned char, 3>  BinaryImageType;
	typedef itk::ConfidenceConnectedImageFilter<ImageType, BinaryImageType> ConfidenceConnectedFilterType;
	ConfidenceConnectedFilterType::Pointer confidenceConnectedFilter_1 = ConfidenceConnectedFilterType::New();
	
	ImageType::IndexType index_1;
	index_1[0] = 275;
	index_1[1] = 250;
	index_1[2] = 159;

	// 参数
	const double varianceMultiplier_1 = 1;
	confidenceConnectedFilter_1->SetInput(image);
	confidenceConnectedFilter_1->SetMultiplier(varianceMultiplier_1);
	confidenceConnectedFilter_1->SetNumberOfIterations(2);
	confidenceConnectedFilter_1->AddSeed(index_1);
	confidenceConnectedFilter_1->SetInitialNeighborhoodRadius(20);
	confidenceConnectedFilter_1->Update();
	confidenceConnectedFilter_1->SetReplaceValue(1);

	double meanEstimation_1 = confidenceConnectedFilter_1->GetMean();
	double varianceEstimation_1 = confidenceConnectedFilter_1->GetVariance();
	std::cout << "Mean estimation_1=" << meanEstimation_1 << std::endl;
	std::cout << "Variance eatimation_1=" << varianceEstimation_1 << std::endl;

	//meanEstimation_1 = 0.5;
	//varianceEstimation_1 = 0.002;


	typedef itk::SimpleFuzzyConnectednessScalarImageFilter<ImageType, BinaryImageType> FuzzySegmentationFilterType;
	FuzzySegmentationFilterType::Pointer fuzzysegmenter_1 = FuzzySegmentationFilterType::New();
	fuzzysegmenter_1->SetInput(image);
	fuzzysegmenter_1->SetObjectSeed(index_1);
	fuzzysegmenter_1->SetMean(meanEstimation_1);
	fuzzysegmenter_1->SetVariance(varianceEstimation_1);
	fuzzysegmenter_1->SetThreshold(0.00001);
	fuzzysegmenter_1->Update();

	typedef unsigned char                            OutputPixelType;
	typedef itk::Image< OutputPixelType, 3 > OutputImageType;

	/*
	double meanTolerance = 0.4;
	double stdTolerance = 0.002;
	typedef  itk::VoronoiSegmentationImageFilter<ImageType, OutputImageType, BinaryImageType> VoronoiSegmentationFilterType;
	VoronoiSegmentationFilterType::Pointer voronoisegmenter = VoronoiSegmentationFilterType::New();
	voronoisegmenter->SetInput(image);
	voronoisegmenter->TakeAPrior(fuzzysegmenter_1->GetOutput());
	voronoisegmenter->SetMeanPercentError(meanTolerance);
	voronoisegmenter->SetSTDPercentError(stdTolerance);
	voronoisegmenter->SetMinRegion(5);
	voronoisegmenter->Update();
	*/
	


	
	typedef itk::RescaleIntensityImageFilter< OutputImageType, OutputImageType > ScalerFilterType;
	ScalerFilterType::Pointer scaler = ScalerFilterType::New();
	scaler->SetOutputMinimum(0);
	scaler->SetOutputMaximum(255);
	scaler->SetInput(fuzzysegmenter_1->GetOutput());
	scaler->Update();



	// 高光谱类型TIF输出
	/*
	typedef  itk::ImageFileWriter< OutputImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName("test.tif");
	writer->SetInput(scaler->GetOutput());
	writer->Update();
	*/

	typedef itk::Image< unsigned char, 3 > ImageType1;



	//Mat格式输出待完善

	//double *outA = new double[M*N*P];
	//for (int k = 0; k < P; k++) {
	//	for (int i = 0; i < M; i++) {
	//		for (int j = 0; j < N; j++) {
	//			ImageType1::IndexType pixelIndex;
	//			pixelIndex[0] = i; // x position 
	//			pixelIndex[1] = j; // y position 
	//			pixelIndex[2] = k; // z position
	//			ImageType1::PixelType pixelValue = scaler->GetOutput()->GetPixel(pixelIndex);

	//			//outA[M*j + i] = A[i][j];
	//		}
	//	}
	//}
	//
	//	
	//		
	//pmatFile = matOpen("A.mat", "w");
	//pMxArray = mxCreateDoubleMatrix(M, N, mxREAL);
	//mxSetData(pMxArray, outA);
	//matPutVariable(pmatFile, "A", pMxArray);
	//matClose(pmatFile);


	
	typedef itk::Image< unsigned char, 2 >     Image2DType;
	typedef itk::ImageSeriesWriter< ImageType1, Image2DType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(scaler->GetOutput());

	typedef itk::NumericSeriesFileNames    NameGeneratorType;
	NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
	
	std::string format = "E:\\test\\Project5\\Result2\\";
	format += "%d.";
	format += "tif";   // filename extension
	nameGenerator->SetSeriesFormat(format.c_str());
	
	ImageType1::ConstPointer inputImage = scaler->GetOutput();
	ImageType1::RegionType   region1 = inputImage->GetLargestPossibleRegion();
	ImageType1::IndexType    start1 = region1.GetIndex();
	
	const unsigned int firstSlice = start[2];
	const unsigned int lastSlice = start[2] + size[2] - 1;

	nameGenerator->SetStartIndex(firstSlice);
	nameGenerator->SetEndIndex(lastSlice);
	nameGenerator->SetIncrementIndex(1);
	
	writer->SetFileNames(nameGenerator->GetFileNames());
	
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & excp)
	{
		std::cerr << "Exception thrown while reading the image" << std::endl;
		std::cerr << excp << std::endl;
	}
	


	return 0;
}