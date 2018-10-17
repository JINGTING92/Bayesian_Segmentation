// ==================================================================== //
//																		//
// Filename: ImageAnalyser.h					        			    //
//																		//
// Author:	Fraunhofer Singapore 			                //
//																		//
// ===================================================================	//
//																		//
// Short description:													//
// Provides functionality for searching points of interest in an image	//
// -------------------------------------------------------------------	//
//																		//
// Creation Date	: 15.06.2009	Marius Erdt							//
// Modified     	: 10.08.2009	Marius Erdt							//
//												
// Modified             : 2017-2018     JINGTING MA
// ===================================================================	//

#pragma once

#include <my_awesomeproject_Segmentation_Export.h>

// mitk includes
#include <mitkSurface.h>
#include <mitkDataNode.h>
#include <mitkImage.h>
#include <mitkTimeGeometry.h>
#include "mitkImageAccessByItk.h"
#include "mitkITKImageImport.h"
#include "mitkDataStorage.h"
#include "mitkRenderingManager.h"
#include <mitkSurfaceToImageFilter.h>
#include <mitkImageTimeSelector.h>
#include <mitkStatusBar.h>
#include <mitkProgressBar.h>
#include <mitkImageCast.h>



// vtk includes
#include <vtkIdList.h>
#include <vtkPolyData.h>
#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataWriter.h>
#include <vtkDoubleArray.h> 
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkPoints.h>
#include <vtkPolyDataNormals.h>
#include <vtkIterativeClosestPointTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkLandmarkTransform.h>
#include <vtkIdList.h>
#include <vtkPolyLine.h>
#include <vtkCell.h>


// itk includes
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include "itkCastImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkRGBPixel.h"
#include "itkScalarToRGBColormapImageFilter.h"
#include "itkCustomColormapFunction.h"
#include "itkPasteImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkLabelShapeKeepNObjectsImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkBinaryFillholeImageFilter.h"
#include "itkBinaryCrossStructuringElement.h"
#include "itkBinaryThresholdImageFilter.h"
#include <itkMedianImageFilter.h>
#include "itkResampleImageFilter.h"
#include <itkNearestNeighborInterpolateImageFunction.h>
#include "itkVolumeGreyLevelStatisticsFilter.h"
#include <itkRegionOfInterestImageFilter.h>
#include <itkCropImageFilter.h>
#include "itkShiftScaleImageFilter.h"
#include "itkSobelEdgeDetectionImageFilter.h"
#include <itkTimeProbe.h>
#include "itkGradientMagnitudeImageFilter.h"


using namespace std; 


#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif


typedef short InternalPixelType;
typedef itk::Image<InternalPixelType, 3> AtlasImageType;
typedef itk::Image<unsigned char, 3> CharImageType; 
typedef itk::StatisticsImageFilter< AtlasImageType > StatisticsImageFilterType; 
typedef itk::ImageFileWriter< AtlasImageType > AtlasImageWriterType; 
typedef itk::ImageDuplicator< AtlasImageType > ImageDuplicatorType; 
typedef itk::VolumeGreyLevelStatisticsFilter< AtlasImageType, AtlasImageType> statisticFilterType;

class AWESOME_EXPORTS ImageAnalyser{

public:

	ImageAnalyser(mitk::DataNode::Pointer imageNode, AtlasImageType::Pointer probMapImage);


	~ImageAnalyser();


	enum OrganType { ABDOMINAL, LUNG, LIVER, KIDNEYS, KIDNEY_LEFT, KIDNEY_RIGHT, HEAD, BONE, SPLEEN, KNEE, BLADDER, STOMACH };


	/**
	 *  @brief    landmark deformation considering of the probability inside and outside
	 *  @param    pass the current shape, probability map from NN, output the deformed shape 
	 *            the probability of ROI and NOI needs pre-training 
	 *            prob_i = p_i + lambda * m_i 
	 *  @comment  Prior knowledge from preTraining: InNOI = 1.30 +/- 0.08 for ground truth + NN-map 
	 *                                             OutROI = 0.02 +/- 0.04 for ground truth + NN-map 
	 */
	void landmarkDeformation(vtkPolyData* _inputPoly, vtkPolyData* outputPoly, double _lambdaNN, int _radius, int & _onBoundary);


	void trainProbsROIandNOI(vtkPolyData* _inputPoly, double _lambdaNN, int _radius); 


	/**
	 *  @brief   compute InROI and InNOI interior organ as threshold initialization (prior knowledge)
	 *  @param   input  : the current shape, probability map from NN
	 *           output : InROI (mean) + InNOI (mean)
	 */
	void computePriorProbs(vtkPolyData* _inputPoly, double _lambdaNN, int _radius);


	/**
	 *  @brief   compute GMM for image intensity of ROI (interior organ) and NOI (outerior organ)
	 *  @param   input  : the current shape, probability map from NN, bool useSurface = true (compute with the input poly)
	 *           output : ROI_MEAN + ROI_STD, NOI_MEAN + NOI_STD 
	 */
	void trainPriorGMM(vtkPolyData* _inputPoly); 



	// getter for the mitk image holding the original image
	mitk::Image::Pointer GetMitkImage()
	{
		return m_mitkImage; 
	}


	// getter for the mitk image holding the original image
	mitk::DataNode::Pointer GetImageNode()
	{
		m_imageNode->SetData(m_mitkImage); 
		m_imageNode->Update(); 
		return m_imageNode; 
	}


	AtlasImageType::Pointer GetNOIImage()
	{
		return m_NOIImage; 
	}


protected:

	bool checkBoundary(itk::Point< double, 3 > _point, double _lambdaNN, int _radius); 

	void calculateProbInside(AtlasImageType::IndexType _baseIndex, double _lambdaNN, int _radius, double& ProbROI); 

	void calculateProbOutside(AtlasImageType::IndexType _baseIndex, double _lambdaNN, int _radius, double & ProbNOI);

	void calculateProbOfNeighborhood(AtlasImageType::IndexType _baseIndex, double _lambda, int _radius, double & ProbROI, double & ProbNOI); 

	
	// >>>>>>>>>>>>>>>>>>>>>>


	/**
	 *  @brief    nonlinear mapping intensity in probability map 
	 *  @param    if x < lowThres,   p = 0 
	 *            if x >= highThres, p = 1
	 *            if lowThres <= x < highThres,  p = (x - lowThres) / (highThres - lowThres) * (1.0 - 0.0) 
	 */
	void RescaleIntensityProbMap();

	mitk::Image::Pointer	m_mitkImage;	// original image in mitk format

	mitk::DataNode::Pointer m_imageNode;	// node containing the original image

	AtlasImageType::Pointer m_itkImage;	    // itk version of the original image casted to the internal pixel type of this class
	AtlasImageType::Pointer m_edgeImage;	// Used in searchInterestPoints : edge image of m_itkImage in itk format

	AtlasImageType::Pointer m_probMapImage;  // probability map image used in searchInterestPoints 
	AtlasImageType::Pointer m_maskImage;     // extracted image from cropped mask 

	AtlasImageType::Pointer m_ROIImage; 
	AtlasImageType::Pointer m_NOIImage; 

	std::string m_nameSegImg;    // node name for the image to segment 
	int m_closeRadius; 


	static const int s = 3; // for edge detection mask

	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	double PriorProb_ROI; 
	double PriorProb_ROI_STD;
	double PriorProb_NOI; 
	double PriorProb_NOI_STD; 


	double ROI_MEAN; 
	double ROI_STD; 
	double NOI_MEAN; 
	double NOI_STD; 


	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


	double getGaussianProbOfROI(InternalPixelType _pIntens);


	double getGaussianProbOfNOI(InternalPixelType _pIntens);


	double getPriorProbBasedOnNNMap(AtlasImageType::IndexType _inputIndex); 

};
