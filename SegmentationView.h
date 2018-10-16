
#ifndef SegmentationView_h
#define SegmentationView_h

// class / plugin_related 
#include "ui_SegmentationViewControls.h"

#include "AffineSurfaceInteractor.h"
#include "SurfaceDeformInteractor.h"
#include "ImageAnalyser.h"
#include "ModelDeformator.h"
#include "ShapeModel.h"
#include "ShapeModelProjection.h"

#include "RobustKPCA.h"
#include "StatisticalUtilities.h"


// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>

// Qmitk
#include "QmitkStdMultiWidget.h"
#include <QmitkAbstractView.h>


// qt includes
#include <qdir.h>
#include <QMessageBox>
#include <qaction.h>
#include <QMessageBox.h>
#include <qapplication.h>
#include <qcombobox.h>
#include <qcheckbox.h>
#include <qslider.h>
#include <qspinbox.h>
#include <qradiobutton.h>
#include <qpushbutton.h>
#include <QFileDialog>
#include <qlistview.h>
#include <qtreeview.h>

// vtk include
#include <vtkPolyData.h>
#include <vtkMatrix4x4.h>
#include <vtkSmartPointer.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkLinearTransform.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h> 
#include <vtkCurvatures.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataWriter.h>
#include <vtkCleanPolyData.h>
#include <vtktransform.h>
#include <vtkPolyDataReader.h>
#include <vtkReflectionFilter.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPolyDataConnectivityFilter.h> 
#include <vtkProbeFilter.h> 
#include <vtkDataSet.h> 
#include <vtkDataSetReader.h> 
#include <vtkCenterOfMass.h>
#include <vtkPointLocator.h> 
#include <vtkIterativeClosestPointTransform.h>
#include <vtkMatrix4x4.h> 
#include <vtkSmoothPolyDataFilter.h> 



  // new added for remeshing 
//#include "Common\vtkSurface.h"
#include "DiscreteRemeshing/vtkIsotropicDiscreteRemeshing.h"


// itk include 
#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageRegistrationMethod.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSpatialObjectToImageFilter.h"
#include "itkTranslationTransform.h"
#include "itkThresholdImageFilter.h"
#include "itkPasteImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageSeriesWriter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkAffineTransform.h"
#include "itkMutualInformationImageToImageMetric.h"
#include "itkImageRegistrationMethodv4.h"
#include "itkCenteredTransformInitializer.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkVersorRigid3DTransform.h"
#include "itkVersorRigid3DTransformOptimizer.h"
#include "itkFlatStructuringElement.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include <itkBinaryThresholdImageFilter.h>
#include "itkTextOutput.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkLabelOverlapMeasuresImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelToRGBImageFilter.h"
#include "itkLabelShapeKeepNObjectsImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkLabelMapToBinaryImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkGeometryUtilities.h"
#include "itkHausdorffDistanceImageFilter.h"
#include "itkImageMomentsCalculator.h"
#include "itkChangeInformationImageFilter.h"
#include "itkPasteImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkBinaryContourImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkNormalizeImageFilter.h"
#include "itkCurvatureFlowImageFilter.h"
#include "itkIsolatedConnectedImageFilter.h"
#include "src/itkVolumeGreyLevelStatisticsFilter.h"
#include "itkAddImageFilter.h"


// mitk includes
#include "mitkGlobalInteraction.h"
#include "mitkSurfaceToImageFilter.h"
#include "mitkImageTimeSelector.h"
#include "mitkITKImageImport.h"
#include "mitkLookupTable.h"
#include "mitkLookUpTableProperty.h"
#include "mitkMaterial.h"
#include "mitkDataNodeFactory.h"
#include "mitkImageToItk.h"
#include "mitkImageCast.h"
#include "mitkNodePredicateDimension.h"
#include "mitkNodePredicateDataType.h"
#include "mitkNodePredicateAnd.h"
#include "mitkNodePredicateNot.h"
#include "mitkNodePredicateProperty.h"
#include "mitkStatusBar.h"
#include "mitkProgressBar.h"
#include "mitkProperties.h"
#include "mitkVtkRepresentationProperty.h"
#include <mitkStandardFileLocations.h>
#include "mitkLog.h"
#include "mitkImageToSurfaceFilter.h"
#include "J:/MITK_SOURCE\MITK-2014.10.0/Modules/Segmentation/Algorithms/mitkManualSegmentationToSurfaceFilter.h"
#include "mitkBoundingObject.h"
#include "mitkBoundingObjectCutter.h"
#include "mitkCuboid.h"
#include "mitkLookupTableProperty.h"


// general includes 
#include <random>
#include <stdlib.h> 
#include <time.h>

#define USE_PCA


/*!
\brief functionality that provides methods for model based organ segmentation

\sa QmitkFunctionality
\ingroup Functionalities
*/
class SegmentationView : public QmitkAbstractView
{
	// this is needed for all Qt objects that should have a Qt meta-object
	// (everything that derives from QObject and wants to have signal/slots)
	Q_OBJECT

public:

	typedef itk::Image< short, 3 >	AtlasImageType;	
	typedef itk::Image< short, 2 > Atlas2DImageType; 
	typedef itk::Image< unsigned char, 2 > Image2DType;
	typedef itk::Image< unsigned char, 3 > Image3DType; 
	typedef itk::RGBPixel< unsigned char > RGBPixelType; 
	typedef itk::Image< RGBPixelType, 3 > RGBImageType; 

	typedef vnl_vector<double> VectorOfDouble;
	typedef vnl_vector<int>  VectorOfInteger; 
	typedef vnl_matrix<double>	MatrixOfDouble;
	typedef vnl_matrix<int>  MatrixOfInteger; 

	typedef itk::ImageFileWriter< AtlasImageType > WriterFilterType; 
	typedef itk::ImageDuplicator< AtlasImageType > DuplicatorType;
	typedef itk::ThresholdImageFilter< AtlasImageType > ThresholdFilterType; 
	typedef itk::ExtractImageFilter< AtlasImageType, AtlasImageType > ExtractFilterType;

	static const std::string VIEW_ID;

	SegmentationView();

	virtual ~SegmentationView();

	virtual void CreateQtPartControl(QWidget *parent);


// ***** Deep Learning Utilities ****** ///
protected slots:


	/**
	 * @brief		Select an image node from storage to use as fix image in registration
	 *				Checks whether the node already has an "id" (integer property of the node).
	 *				The id is necessary for remembering the processing state of each node.
	 * @param node_	the new node that has been selected
	 */
	void FixImageSelectedForReg(const mitk::DataNode * node_);


	/**
	 * @brief		Resampling input images to make the same spacing
	 *				Comment code is used for registration
	 * @param       Resampling the pixel with 2.0 by default
	 *              Rescaling the intensity to [0, 255]
	 *
	 */
	void LoadAndResampleImage();


	/**
	 *  @brief   compute GMM of the ROI and NOI 
	 */
	void ComputeGMM(); 


     // ****** Deep Learning Utilities ****** //
     /**
      *  @brief   load 3d images to generate 2d slices of jpg 
      *  @param   select the mode of slices 
      */
     void LoadImagesForSliceGeneration(); 


	 /**
	  * @brief   Normalize the mask images to 0 / 1
	  */
	 void NormalizeMask();


	 /**
	  *  @brief   Load SA Slices and Merge to a 3D-ITK image with dimension 256 x 256 x 256  
	  */
	 void LoadAndMergeSASlices(); 


	 /**
	  *  @brief   Load SC Slices and Merge to a 3D-ITK image with dimension 256 x 256 x 256
	  */
	 void LoadAndMergeSCSlices();


	 /**
	  *  @brief   Load CA Slices and Merge to a 3D-ITK image with dimension 256 x 256 x 256
	  */
	 void LoadAndMergeCASlices();


	 /**
	  *  @brief   merge kidney slices (separate)
	  */
	 void MergeKidneyMap(); 


	 /**
	  *  @brief   Merge the 3D-probability map 
	  */
	 void MergePancreasMap(); 



// ***** Modelling Utilities ****** ///
protected slots:


     // ******** Modelling ************** // 
    
    /**
     * @brief   load a set of polys to train the model 
	 * @param   save polys as m_polysForModelling
     */
    void LoadDatasetsForModelling(); 


	/**
  	 * @brief  select model type
	 * @param  m_modelType = "PCA | KPCA | RPCA | RKPCA"
	 */
	void SelectModelType();


	/**
	 * @brief  train m_shapeModel from m_polysForModelling 
	 *         Set sigma and select model type before modelling 
	 */
	void TrainShapeModel(); 



// ***** Evaluation Utilities ****** ///
protected slots:


     /**
      * @brief       Select the ground truth mask 
      */
     void SelectGroundTruthMask(const mitk::DataNode * node);


	 /** 
	  *  @brief      Select the compared segment 
	  */
	 void SelectSegmentToEvaluate(const mitk::DataNode * node); 


	 /**
	  *  @brief      Start evaluation
	  *              Compare the ground truth mask and segment 
	  *  @param      return: VOE , DSC , Hauffdistance , Jaccard Index
	  */
	 void StartSegmentationEvaluation(); 


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/**
	 * @brief		called when a selection in the drop down image selection menu changes. 
	 *				Checks whether the node already has an "id" (integer property of the node).
	 *				The id is necessary for remembering the processing state of each node.
	 * @param node_	the new node that has been selected
	 */
    void ImageSelected( const mitk::DataNode * node_ );


	void ProbMapSelected(const mitk::DataNode * mapNode_); 


	void ProbShapeSelected(const mitk::DataNode* shapeNode_); 


	void Deform(); 


   /**
	*  @brief	Performs the model adaptation. Checks which organ has to be segmented and
	*			sets parameters accordingly. Calls AutoOrientModel and FineAdaptation for
	*			model adaptation.
	*/
	void Segment();


	/**
	 * @brief  Perform Auto-Orientation without model constraint 
	 */
	void AutoOrientationWithoutModel(); 


	/**
 	 * @brief  Perform Auto-Orientation without model constraint
	 */
	void AutoOrientationWithModel();


/////////////////////////////////////// TOOLBOX ///////////////////////////////////////////////////////////////////////


	/**
	 * @brief   create mesh with defined number of landmarks from binary image
	 * @param   load a set of binary images 
	 *          introduce inline functions: 
	 *                vtkPolyData* createMeshFromBinaryImage( AtlasImageType::Pointer mask )
	 *                vtkPolyData* remeshSurface( vtkPolyData* surfacePoly, ... params )
	 * @param   migrate from mitkACVD
	 */
	void GenerateCorrespondingMeshFromBinaryImage(); 


	/**
	 * @brief   crop kidneys from label data   
	 * @param   create the labels for region of interest 
	 *          save the bounding index to the file name 
	 */
	void FilterRegionOfInterest(); 


	/**
	 * @brief  crop the image with the defined axes - start index
	 * @param  load image from  fixImg_TreeNode
	 * @param     user define x_lower / x_upper ; y_lower / y_upper; z_lower / z_upper
	 *            write output image name with "_cp_"
	 */
	void CropImageWithDefinedAxes(); 


	/**
	 * @brief  reflect mesh for kidney_correspondence establishment 
	 * @param  input mesh_samples 
	 *         output mirrored mesh_samples
	 */
	void PolyDataReflectionX(); 


	/**
	 * @brief  reflect mesh for kidney_correspondence establishment
	 * @param  input mesh_samples
	 *         output mirrored mesh_samples
	 */
	void PolyDataReflectionY();


	/**
	 * @brief  reflect mesh for kidney_correspondence establishment
	 * @param  input mesh_samples
	 *         output mirrored mesh_samples
	 */
	void PolyDataReflectionZ();


	/**
	 * @brief  separate left & right kidneys of mask image 
	 * @param  output: maskName_left.nrrd with the same size, origin, spacing Using connectedComponentFilter 
	 */
	void SeparateMasks(); 


	/** 
	 * @brief    check the GT boundary size and index 
	 */
	void CheckMaskBoundary(); 


	/**
	 * @brief  compute the ROI boundary from mask and then cropping 
	 */
	void ComputeAndCropROI(); 



protected:


	virtual void SetFocus();


	/**
	 *  @brief   load probability map vtk for shape model localization - ICP transform 
	 */
	void loadProbShapeForLocalization(); 


	/**
	 *  @brief  load mean shape from local directory
	 */
	void loadDeformShapeFromLocal();


	/**
	 *  @brief  load ground truth shape for pre-training GMM and PriorProbs in initial Deform()
	 */
	void loadGroundTruthShapeForPreTraining(); 


	/**
	 *  @brief  align the m_deformShape with m_probShape with ICP transform 
	 */
	void alignDeformShapeWithProbShape(vtkPolyData* _sourcePoly, vtkPolyData* _targetPoly);


	/**
	 *  @brief   transform m_deformedShape to m_probShape 
	 */
	void transformDeformShapeToProbMap(vtkPolyData* _inputPoly); 


	void reComputeNormals(vtkPolyData* _poly);


	mitk::Surface::Pointer SegmentationView::GetTransformedGeometrySurface(mitk::Surface::Pointer surface); 
	
	/**
	 *  @brief   create mesh with defined landmark from binary image, then use in correspondence 
	 *           migration from MITK_ACVD 
	 */
	void createMeshFromBinaryImage(Image3DType::Pointer _mask, std::string imageName, int numVertices);


	vtkPolyData* remeshSurface(vtkPolyData* surfacePoly, int numVertices);
	

	void addImageToNode(AtlasImageType::Pointer img, std::string imgName, bool visable); 


	void addCharImageToNode(Image3DType::Pointer img, std::string imgName, bool visable);


	void addPolyToNode(vtkPolyData* poly, std::string meshName, bool visable); 


	void CenterOfMassToOrigin(std::vector< vtkSmartPointer< vtkPolyData >> &m_meshes, bool weighting);


	void visualizeDatasets(std::vector< vtkSmartPointer< vtkPolyData > > &m_datasets, std::string name0, bool visible);


	void writeAtlasImage(AtlasImageType::Pointer img, std::string writeName); 


	void writeChar3DImage(Image3DType::Pointer img, std::string writeName); 


////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 
	Ui::SegmentationViewControls m_Controls;


	void AutoSearchBasedOnIntensModel(int iterationCount, double _lambda, int _radius); 


	void AutoRefineWithModel(int iterationCount, bool _useWeights);


	/**
	 * @brief   reset poly data 
	 * @note    delete and create a new polyData 
	 */
	void ResetPolyData( vtkSmartPointer< vtkPolyData > poly ); 


	/**
	 * @brief   write out poly data 
	 */
	void WriteOutPolyData(vtkSmartPointer< vtkPolyData > poly, std::string fileName); 


////////////////////////////////////////////////////////////////////////////////////////////////////////

	/**
	 * @brief   generate slices with the fimeName, sliceMode, write into the folder 
	 */
	void generateSlicesWithMode(std::string fileName, std::string sliceMode, Image3DType::Pointer img);


	//************ segmentation member ***********//

	
	mitk::DataNode::Pointer m_currentSurfaceNode;	/**> this nodes holds the latest (adapted) version of the model **/

	ImageAnalyser*   m_imageAnalyser;				/**> upon image : boundary search **/
	AffineSurfaceInteractor::Pointer m_AffineSurfaceInteractor;  /**> Interactor for rigid transformation **/

	vtkSmartPointer<vtkPolyData>	m_icpPoly;		  /**> reference shape - an affine transformed version of the initial model **/
	vtkSmartPointer<vtkPolyData>	m_structurePoly;  /**> stores the found boundary candidates  **/
	vtkSmartPointer<vtkPolyData>    m_reconsRightKidney;  /**> stores the second data as the template for right kidney segmentation **/

	mitk::Surface::Pointer m_icpSurface;			/**> mitk surface to store m_icpPoly **/ 
	mitk::Surface::Pointer m_structureSurface;		/**> mitk surface to store m_structurePoly **/

	AtlasImageType::Pointer m_GTMask;               /**> ground truth mask for evaluation **/
	AtlasImageType::Pointer m_SegmentMask;          /**> segment for evaluation **/

	AtlasImageType::Pointer m_segImage;             /**> loaded seg image  **/
	AtlasImageType::Pointer m_probMap;              /**> current probability map for segmentation **/
	vtkSmartPointer< vtkPolyData > m_probShape;     /**> probability shape for initial localization , but not necessary**/
	vtkSmartPointer< vtkPolyData > m_deformedShape; /**> deformed shape**/
	vtkSmartPointer< vtkPolyData > m_preTrainShape; /**> GT shape for pre-training **/

	QDir* m_modelsDir;	                            /**> J:/IGDMedApp_VS12/build/bin/Release/Models **/

	mitk::DataNode::Pointer m_visualizationSurfaceNode;	// node to store the adapted model. This node can then be used by other plugins.


	//************ Added member ***********//	

	std::string m_recentDir;

	mitk::Image* m_fixedImageForReg; 

	std::string m_fixedImgName;  
	std::string m_segImageName; 
	std::string m_probMapName; 

	int m_shapeOrder; 


	RobustKPCA* m_shapeModel; /**> current Shape Model **/

	RobustKPCA* m_PCAModel;   /**> always exists for shape constraint **/

	std::string m_modelType;  /**> selected model type **/

	double m_sigmaKPCA;       /**> sigma for kernel models**/

	std::vector< vtkSmartPointer< vtkPolyData > > m_trainingDatasets;  /**> polys for training models **/

	mitk::LookupTableProperty::Pointer m_lookupTableProp;

};

#endif // _SegmentationView_H_INCLUDED


/** Image Gravity Load Shape **/
//typedef itk::ImageMomentsCalculator< AtlasImageType > ImageCalculatorType;
//ImageCalculatorType::Pointer calculator = ImageCalculatorType::New();
//calculator->SetImage(m_imageAnalyser->GetMaskImage());
//calculator->Compute();
//double gravity[3];
//gravity[0] = calculator->GetCenterOfGravity()[0];
//gravity[1] = calculator->GetCenterOfGravity()[1];
//gravity[2] = calculator->GetCenterOfGravity()[2];
//std::cout << __FUNCTION__ << " gravity from calculator = " << gravity[0] << " , " << gravity[1] << " , " << gravity[2] << std::endl;
//mitk::BaseGeometry::Pointer modelGeometry = m_currentSurfaceNode->GetData()->GetTimeGeometry()->GetGeometryForTimeStep(0);
//mitk::Point3D x;
//x[0] = gravity[0];
//x[1] = gravity[1];
//x[2] = gravity[2];
//GetRenderWindowPart()->SetSelectedPosition(x);
//mitk::Vector3D t;
//t = x - modelGeometry->GetCenter();
//modelGeometry->Translate(t);
//modelGeometry->Modified();

