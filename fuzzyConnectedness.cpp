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
#include "itkTIFFImageIOFactory.h"
#include "itkPermuteAxesImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkCurvatureFlowImageFilter.h"
#include "vtkAutoInit.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"

using namespace std;

int main(int argc,char * argv[])
{
	VTK_MODULE_INIT(vtkRenderingVolumeOpenGL2);
	VTK_MODULE_INIT(vtkRenderingOpenGL2);
	VTK_MODULE_INIT(vtkInteractionStyle);
	itk::TIFFImageIOFactory::RegisterOneFactory();
	typedef itk::Image<double, 3> ImageType;
	MATFile *pmatFile = NULL;
	mxArray *pMxArray = NULL;
	//open .mat file 
	double *S1;
	_int64 M, N, P;
	pmatFile = matOpen("E:\\WWZResult.mat", "r");//Hessian matrix from matlab
	pMxArray = matGetVariable(pmatFile, "WWZResult");
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
	typedef itk::Image<unsigned char, 3>  BinaryImageType;
	typedef unsigned char LabelPixelType;
	typedef itk::Image<LabelPixelType, 3>  LabelImageType;
 	typedef itk::ConfidenceConnectedImageFilter<ImageType, BinaryImageType> ConfidenceConnectedFilterType;
	ConfidenceConnectedFilterType::Pointer confidenceConnectedFilter_1 = ConfidenceConnectedFilterType::New();

	typedef itk::SimpleFuzzyConnectednessScalarImageFilter<ImageType, LabelImageType> FuzzySegmentationFilterType;
	FuzzySegmentationFilterType::Pointer fuzzysegmenter_1 = FuzzySegmentationFilterType::New();
	
	ImageType::IndexType index_1;//seed point in the pulmonary artery trunk
	index_1[0] = 213;
	index_1[1] = 275;
	index_1[2] = 103;

	const double varianceMultiplier_1 = 3;
	confidenceConnectedFilter_1->SetInput(image);
	confidenceConnectedFilter_1->SetMultiplier(varianceMultiplier_1);
	confidenceConnectedFilter_1->SetNumberOfIterations(2);
	confidenceConnectedFilter_1->AddSeed(index_1);
	confidenceConnectedFilter_1->SetInitialNeighborhoodRadius(1);

	confidenceConnectedFilter_1->Update();
	confidenceConnectedFilter_1->SetReplaceValue(1);
	const double meanEstimation_1 = confidenceConnectedFilter_1->GetMean();
	const double varianceEstimation_1 = confidenceConnectedFilter_1->GetVariance();
	

	fuzzysegmenter_1->SetInput(image);
	std::cout << "Mean estimation_1=" << meanEstimation_1 << std::endl;
	std::cout << "Variance eatimation_1=" << varianceEstimation_1 << std::endl;

	fuzzysegmenter_1->SetObjectSeed(index_1);
	fuzzysegmenter_1->SetMean(meanEstimation_1);
	fuzzysegmenter_1->SetVariance(varianceEstimation_1);
	fuzzysegmenter_1->SetThreshold(0.08);

	fuzzysegmenter_1->Update();
	
	/*Visualize*/
	typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;
	ConnectorType::Pointer connector = ConnectorType::New();
	connector->SetInput(fuzzysegmenter_1->GetOutput());
	connector->Update();//vtkImageData

	vtkImageCast *readerImageCast = vtkImageCast::New();
	readerImageCast->SetInputData(connector->GetOutput());

	vtkMarchingCubes *skinExtractor = vtkMarchingCubes::New();
	skinExtractor->SetInputConnection(readerImageCast->GetOutputPort());
	skinExtractor->SetValue(0, 0.2);

	vtkDecimatePro *deci = vtkDecimatePro::New();
	deci->SetTargetReduction(0.3);
	deci->SetInputConnection(skinExtractor->GetOutputPort());

	vtkSmoothPolyDataFilter *smooth = vtkSmoothPolyDataFilter::New();  //使图像更加光滑
	smooth->SetInputConnection(deci->GetOutputPort());
	smooth->SetNumberOfIterations(500);

	vtkPolyDataNormals *skinNormals = vtkPolyDataNormals::New();//绘制法线
	skinNormals->SetInputConnection(smooth->GetOutputPort());

	vtkPolyDataMapper *skinMapper = vtkPolyDataMapper::New();
	skinMapper->SetInputConnection(skinNormals->GetOutputPort());
	skinMapper->ScalarVisibilityOff();

	vtkActor *skin = vtkActor::New();
	skin->SetMapper(skinMapper); //获得皮肤几何数据的属性 
	skin->GetProperty()->SetDiffuseColor(.56, .38, .29); //设置皮肤颜色的属性 
	skin->GetProperty()->SetSpecular(.4); //设置反射率 
	skin->GetProperty()->SetSpecularPower(20); //设置反射光强度
	skin->GetProperty()->SetOpacity(0.5);//设置透明度

	vtkCamera *aCamera = vtkCamera::New(); //定义摄像机
	aCamera->SetViewUp(0, 0, -1); //取得摄像机方向
	aCamera->SetPosition(0, 1, 0); //光源位置
	aCamera->SetFocalPoint(0, 0, 0); //取焦点坐标
	aCamera->ComputeViewPlaneNormal();

	vtkRenderer *aRenderer = vtkRenderer::New();
	vtkRenderWindow *renWin = vtkRenderWindow::New();
	renWin->AddRenderer(aRenderer);
	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);

	aRenderer->AddActor(skin); //渲染皮肤 此步奏很关键必须有
	aRenderer->SetActiveCamera(aCamera);
	aRenderer->ResetCamera();
	aCamera->Dolly(1.5); //大于1向摄像机焦点移动小于1则向远离焦点的方向移动
	aRenderer->SetBackground(1, 1, 1);
	renWin->SetSize(640, 700);
	aRenderer->ResetCameraClippingRange(); //裁剪

	iren->Initialize();
	iren->Start();

	getchar();
	return EXIT_SUCCESS;
}
