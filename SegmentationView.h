
#ifndef SegmentationView_h
#define SegmentationView_h

// class / plugin_related 
#include "ui_SegmentationViewControls.h"


#include "ImageAnalyser.h"
#include "ModelDeformator.h"
#include "ShapeModel.h"


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
	

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 
	Ui::SegmentationViewControls m_Controls;


	void AutoSearchBasedOnIntensModel(int iterationCount, double _lambda, int _radius); 


	void AutoRefineWithModel(int iterationCount, bool _useWeights);


	/**
	 * @brief   reset poly data 
	 * @note    delete and create a new polyData 
	 */
	void ResetPolyData( vtkSmartPointer< vtkPolyData > poly ); 


//////////////////////////////////////////////////////////////////////////////////////////

	
	//************ segmentation member ***********//

	
	mitk::DataNode::Pointer m_currentSurfaceNode;	/**> this nodes holds the latest (adapted) version of the model **/

	ImageAnalyser*   m_imageAnalyser;				/**> upon image : boundary search **/
	AffineSurfaceInteractor::Pointer m_AffineSurfaceInteractor;  /**> Interactor for rigid transformation **/

	vtkSmartPointer<vtkPolyData>	m_icpPoly;		  /**> reference shape - an affine transformed version of the initial model **/
	vtkSmartPointer<vtkPolyData>	m_structurePoly;  /**> stores the found boundary candidates  **/
	
	mitk::Surface::Pointer m_icpSurface;			/**> mitk surface to store m_icpPoly **/ 
	mitk::Surface::Pointer m_structureSurface;		/**> mitk surface to store m_structurePoly **/

	AtlasImageType::Pointer m_segImage;             /**> loaded seg image  **/
	AtlasImageType::Pointer m_probMap;              /**> current probability map for segmentation **/
	vtkSmartPointer< vtkPolyData > m_probShape;     /**> probability shape for initial localization , but not necessary**/
	vtkSmartPointer< vtkPolyData > m_deformedShape; /**> deformed shape**/
	
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

