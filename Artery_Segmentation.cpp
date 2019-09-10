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
	// 因为输出成TIFF图，所以需要注册TIFF图像组件
	itk::TIFFImageIOFactory::RegisterOneFactory();

	// 定义double类型3维图像数据
	typedef itk::Image<double, 3> ImageType;

	// 定义mat的file
	MATFile *pmatFile = NULL;
	// 定义mat的矩阵
	mxArray *pMxArray = NULL;
	double *S1;
	// 定义长宽高
	_int64 M, N, P;
	//打开mat文件
	pmatFile = matOpen("WWZResult3.mat", "r");//Hessian matrix from matlab
	// 判断文件是否正常打开
	if (pmatFile == NULL) {
		cout << "MatOpen error!!!" << endl;
	}
	// 把mat文件写入数组
	pMxArray = matGetVariable(pmatFile, "t");

	// 转成double类型后续操作
	S1 = (double*)mxGetData(pMxArray);
	// 得到数组的长宽高
	const mwSize *distrSize = mxGetDimensions(pMxArray);
	M = distrSize[0];
	N = distrSize[1];
	P = distrSize[2];

	
	// 重建3D图像
	// 建立itk的图像类型
	ImageType::Pointer image = ImageType::New();
	//定义图像的开始坐标
	ImageType::IndexType start;
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	//定义每个轴的大小
	ImageType::SizeType size;
	size[0] = M;
	size[1] = N;
	size[2] = P;
	// 定义图像的区域
	ImageType::RegionType region;

	// 设置每个轴的大小
	region.SetSize(size);
	//设置起始坐标
	region.SetIndex(start);
	// 设置图像区域
	image->SetRegions(region);
	image->Allocate();
	//把数组的每一个点的数据写入到图像
	for (int k = 0;k < size[2];k++)
	{
		for (int j = 0;j < size[1];j++)
		{
			for (int i = 0;i < size[0];i++)
			{
				//定义图像的点类型
				ImageType::IndexType point_temp;
				// 赋值
				point_temp[0] = i;//x
				point_temp[1] = j;//y
				point_temp[2] = k;//z
				image->SetPixel(point_temp, S1[k*M*N + j*N + i]);
			}
		}
	}
	// 激活上述操作
	image->Update();



	/*FuzzyConnected,rough Segmentation*/
	
	/*typedef unsigned char LabelPixelType;
	typedef itk::Image<LabelPixelType, 3>  LabelImageType;*/

	// 定义uint8 的二值化图像
	typedef itk::Image<unsigned char, 3>  BinaryImageType;
	// 定义itk的置信滤波器
	typedef itk::ConfidenceConnectedImageFilter<ImageType, BinaryImageType> ConfidenceConnectedFilterType;
	ConfidenceConnectedFilterType::Pointer confidenceConnectedFilter_1 = ConfidenceConnectedFilterType::New();
	
	// 定义图像中的一个点
	ImageType::IndexType index_1;
	index_1[0] = 275;
	index_1[1] = 250;
	index_1[2] = 159;

	// 置信连接
	const double varianceMultiplier_1 = 1;
	confidenceConnectedFilter_1->SetInput(image);
	confidenceConnectedFilter_1->SetMultiplier(varianceMultiplier_1);
	confidenceConnectedFilter_1->SetNumberOfIterations(2);
	confidenceConnectedFilter_1->AddSeed(index_1);
	confidenceConnectedFilter_1->SetInitialNeighborhoodRadius(20);
	confidenceConnectedFilter_1->Update();
	confidenceConnectedFilter_1->SetReplaceValue(1);

	// 得到置信滤波的均值和方差
	double meanEstimation_1 = confidenceConnectedFilter_1->GetMean();
	double varianceEstimation_1 = confidenceConnectedFilter_1->GetVariance();
	std::cout << "Mean estimation_1=" << meanEstimation_1 << std::endl;
	std::cout << "Variance eatimation_1=" << varianceEstimation_1 << std::endl;

	//meanEstimation_1 = 0.5;
	//varianceEstimation_1 = 0.002;

	//定义itk的模糊连接
	typedef itk::SimpleFuzzyConnectednessScalarImageFilter<ImageType, BinaryImageType> FuzzySegmentationFilterType;
	FuzzySegmentationFilterType::Pointer fuzzysegmenter_1 = FuzzySegmentationFilterType::New();
	// 设置输入图像
	fuzzysegmenter_1->SetInput(image);
	// 设置种子点
	fuzzysegmenter_1->SetObjectSeed(index_1);
	// 设置均值，方差， 阈值
	fuzzysegmenter_1->SetMean(meanEstimation_1);
	fuzzysegmenter_1->SetVariance(varianceEstimation_1);
	fuzzysegmenter_1->SetThreshold(0.00001);
	// 激活模糊连接器的运行
	fuzzysegmenter_1->Update();

	// 定义输出图像的数据类型维uint8
	typedef unsigned char                            OutputPixelType;
	typedef itk::Image< OutputPixelType, 3 > OutputImageType;

	// 维诺图连接，可以加上，但需要修改
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
	


	// 定义itk的图像double到0-255的归一化拉伸变化
	typedef itk::RescaleIntensityImageFilter< OutputImageType, OutputImageType > ScalerFilterType;
	ScalerFilterType::Pointer scaler = ScalerFilterType::New();
	scaler->SetOutputMinimum(0);
	scaler->SetOutputMaximum(255);
	// 输入为模糊连接的结果图
	scaler->SetInput(fuzzysegmenter_1->GetOutput());
	scaler->Update();



	// 输出3维的tif图，一共一张
	/*
	typedef  itk::ImageFileWriter< OutputImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName("test.tif");
	writer->SetInput(scaler->GetOutput());
	writer->Update();
	*/

	// 定义itk 3维的uint8图像类型
	typedef itk::Image< unsigned char, 3 > ImageType1;



	//直接输出成mat类型，需小改

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


	// 定义itk 2维uint8图像类型
	typedef itk::Image< unsigned char, 2 >     Image2DType;
	// 定义一个itk图像系列数据处理器，作用是3维uint8图像，输出成n个多张的2维uint8
	typedef itk::ImageSeriesWriter< ImageType1, Image2DType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	// 设置这个处理器的输入维归一化后的3维图像数据
	writer->SetInput(scaler->GetOutput());

	// 定义一个输出2维图像名的系列字符串处理器
	typedef itk::NumericSeriesFileNames    NameGeneratorType;
	NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
	
	// 图像文件名
	std::string format = "E:\\test\\Project5\\Result2\\";
	format += "%d.";
	format += "tif";   // filename extension
	// 设置字符串输出处理器的输入字符串为format
	nameGenerator->SetSeriesFormat(format.c_str());
	
	// 赋值图像系列处理器的输入图像， 每张2维子图的区域， 以及每张2维子图开始的索引
	ImageType1::ConstPointer inputImage = scaler->GetOutput();
	ImageType1::RegionType   region1 = inputImage->GetLargestPossibleRegion();
	ImageType1::IndexType    start1 = region1.GetIndex();
	
	const unsigned int firstSlice = start[2];
	const unsigned int lastSlice = start[2] + size[2] - 1;

	// 定义文件名字符串处理器的开始和结束的数字，以及每次累加1
	nameGenerator->SetStartIndex(firstSlice);
	nameGenerator->SetEndIndex(lastSlice);
	nameGenerator->SetIncrementIndex(1);
	
	// 设置输出图像名
	writer->SetFileNames(nameGenerator->GetFileNames());
	
	try
	{
		// 写出每张2维子图
		writer->Update();
	}
	catch (itk::ExceptionObject & excp)
	{
		std::cerr << "Exception thrown while reading the image" << std::endl;
		std::cerr << excp << std::endl;
	}
	


	return 0;
}