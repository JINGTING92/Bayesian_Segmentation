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
	 *  @brief                  search boundary landmark within candidate boundary or normal profiles
	 *                          reweights the probability of each voxel iteratively with Maximum A Posterior 
	 *                          p( label | I ) -> p( I | label ) * p( label ) where the prior p( label ) is fixed
	 *  @param                  do not use any prior parameters
	 */
	void SearchPointsOfInterest(vtkPolyData* modelPoly, vtkPolyData* structurePoly); 

	
   /**
	* @brief					calls different methods of boundary detection depending on the organID and the number of layers.
	* @param modelPoly			the reference polydata of the model
	*        structurePoly		the polydata to store the boundary points          (point ids and number of points must be the same as in the reference polydata)
	*        alphaWeights		external point weights for each point in the model (point ids must be the same as in the reference polydata)
	*        internalWeights	internal point weights for each point in the model (point ids must be the same as in the reference polydata)
	*        organID			organ identifier
	*/
	void SearchPointsOfInterestKidney(vtkPolyData* modelPoly, vtkPolyData* structurePoly, vtkFloatArray* alphaWeights, vtkFloatArray* internalWeights);


	/**
	* @brief				Sets model constraints automatically by cutting out the image under the current position of the model
	*						and calculating the 25% and 75% quantiles.
	*                       These quantiles are set as lower and upper HU values in the models.
	* @param     surface		   the model (layer = 1 by default)
	*            organID		   the id of the organ
	*/
	void AutoSetKidneyModelConstraints(mitk::Surface::Pointer surface);


/////////////////////////////////////////////////////////////////////////////////////////////////////////


	/**
	* @brief					calls different methods of boundary detection depending on the organID and the number of layers.
	* @param modelPoly			the reference polydata of the model
	*        structurePoly		the polydata to store the boundary points          (point ids and number of points must be the same as in the reference polydata)
	*        alphaWeights		external point weights for each point in the model (point ids must be the same as in the reference polydata)
	*        internalWeights	internal point weights for each point in the model (point ids must be the same as in the reference polydata)
	*        organID			organ identifier
	*/
	void SearchPointsOfInterestPancreas(vtkPolyData* modelPoly, vtkPolyData* structurePoly, vtkFloatArray* alphaWeights, vtkFloatArray* internalWeights);


	/**
	* @brief				Sets model constraints automatically by cutting out the image under the current position of the model
	*						and calculating the 25% and 75% quantiles.
	*                       These quantiles are set as lower and upper HU values in the models.
	* @param     surface		   the model (layer = 1 by default)
	*            organID		   the id of the organ
	*/
	void AutoSetPancreasModelConstraints(mitk::Surface::Pointer surface);


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



	/**
  	 * @brief    Enhance intensity within region of prob_map
	 *           new_intensity = intensity + 300
	 * @param    Previous kidney intensity range : -30 - 170
	 */
	void EnhanceIntensityWithProbMap();



	/**
	 * @brief	Computes an edge image based on the m_itkImage. 
	 * @param   A sobel filter mask is used.
	 */
	void ComputeEdgeImage(AtlasImageType::Pointer _img);



	/**
	 * @brief	Initializes the internal filter mask with sobel values.
	 */
	void InitEdgeMatrix();


	/**
	 * @brief			    Gets the image intensity at a given point in physical space
	 * @param  point		the point coordinates in itk physical space (in mm)
	 *         intensity	the intensity at the point position
	 *         image		the itk input image
	 */
	void getPixelIntensity(double point[3], double& intensity, AtlasImageType::Pointer image);


	/**
	 * @brief	  Iterates over a neighborhood of a given voxel and checks whether at least 
	 *                     1/4 of all voxels are between HUmin and HUmax.
	 *				       FIX: Percentage should be a parameter of this function
	 * @param     point	   the point coordinates in itk physical space (in mm)
	 * @return		       true if the condition is met, false otherwise
	 */
	bool pointIsInside(double point[3], float& HUmin, float& HUmax);


	bool pointIsOnBoundary(double point[3]); 



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

	// edge image of m_itkImage in itk format
	float mX[s][s][s];
	float mY[s][s][s];
	float mZ[s][s][s];


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
