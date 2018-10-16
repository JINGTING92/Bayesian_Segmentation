// Main Entry of Plugin 
#include "SegmentationView.h"

// cpp files 
#include "J:/MITK_MSVC2013_X64/MITK-2014.10.0_VS12/ITK-src/Modules/Numerics/Optimizers/src/itkRegularStepGradientDescentOptimizer.cxx"
#include "J:/MITK_MSVC2013_X64/MITK-2014.10.0_VS12/ITK-src/Modules/Numerics/Optimizers/src/itkVersorRigid3DTransformOptimizer.cxx"
#include "J:/MITK_MSVC2013_X64/MITK-2014.10.0_VS12/ITK-src/Modules/Numerics/Optimizers/src/itkOptimizer.cxx"
#include "J:/MITK_MSVC2013_X64/MITK-2014.10.0_VS12/ITK-src/Modules/Numerics/Optimizers/src/itkRegularStepGradientDescentBaseOptimizer.cxx"
#include "J:/MITK_MSVC2013_X64/MITK-2014.10.0_VS12/ITK-src/Modules/Numerics/Optimizers/src/itkGradientDescentOptimizer.cxx"
#include "J:/MITK_MSVC2013_X64/MITK-2014.10.0_VS12/ITK-src/Modules/Numerics/Optimizers/src/itkSingleValuedNonLinearOptimizer.cxx"
#include "J:/MITK_MSVC2013_X64/MITK-2014.10.0_VS12/ITK-src/Modules/Filtering/LabelMap/src/itkGeometryUtilities.cxx"
#include "J:/MITK_SOURCE/MITK-2014.10.0/Modules/Segmentation/Algorithms/mitkManualSegmentationToSurfaceFilter.cpp"


const std::string SegmentationView::VIEW_ID = "igd.a7.segmentation.segmentationview";


struct ClustersQuadrics
{
	explicit ClustersQuadrics(int size)
		: Elements(new double*[size]),
		Size(size)
	{
		for (int i = 0; i < size; ++i)
		{
			Elements[i] = new double[9];

			for (int j = 0; j < 9; ++j)
				Elements[i][j] = 0.0;
		}
	}

	~ClustersQuadrics()
	{
		for (int i = 0; i < Size; ++i)
			delete[] Elements[i];

		delete Elements;
	}

	double** Elements;
	int Size;

private:

	ClustersQuadrics(const ClustersQuadrics&);

	ClustersQuadrics& operator=(const ClustersQuadrics&);
};


SegmentationView::SegmentationView()
	: QmitkAbstractView(), m_shapeModel(nullptr)
{
	m_imageAnalyser   = 0;
	m_AffineSurfaceInteractor = 0;
	m_currentSurfaceNode = 0;

	m_icpPoly = 0;
	m_structurePoly = 0;
	m_icpSurface = 0;
	m_structureSurface = 0;
	m_reconsRightKidney = 0; 

	m_GTMask = AtlasImageType::New(); 
	m_SegmentMask = AtlasImageType::New(); 
	m_segImage = AtlasImageType::New();
	m_probMap = AtlasImageType::New(); 

	m_fixedImgName = ""; 
	m_segImageName = ""; 
	m_probMapName = ""; 
	m_shapeOrder = 0; 

	m_visualizationSurfaceNode = mitk::DataNode::New();
	m_visualizationSurfaceNode->SetBoolProperty("produced", true ); // the visualization node is intermediate data
	m_visualizationSurfaceNode->SetColor( 1.0, 0.66, 0.49 );
	m_visualizationSurfaceNode->SetOpacity( 0.3 );
	m_visualizationSurfaceNode->SetIntProperty( "layer", 99 );
	m_visualizationSurfaceNode->SetVisibility( false );

	// disable itk warnings
	itk::OutputWindow::GlobalWarningDisplayOff();

	m_recentDir = "J:\\Abdomen Datasets";
}


SegmentationView::~SegmentationView()
{
}


void SegmentationView::SetFocus()
{
  
}


void SegmentationView::CreateQtPartControl( QWidget *parent )
{
	// build up qt view, unless already done
	if ( true )
	{
		m_Controls.setupUi( parent );

		connect( (QObject*)(m_Controls.loadImageForSlicing), SIGNAL(clicked()), (QObject*)this, SLOT(LoadImagesForSliceGeneration())); 
		connect( (QObject*)(m_Controls.normalizeMask), SIGNAL(clicked()), (QObject*)this, SLOT(NormalizeMask()));
		connect( (QObject*)(m_Controls.loadSASlicesAndMerge), SIGNAL(clicked()), (QObject*)this, SLOT(LoadAndMergeSASlices())); 
		connect( (QObject*)(m_Controls.loadSCSlicesAndMerge), SIGNAL(clicked()), (QObject*)this, SLOT(LoadAndMergeSCSlices()));
		connect( (QObject*)(m_Controls.loadCASlicesAndMerge), SIGNAL(clicked()), (QObject*)this, SLOT(LoadAndMergeCASlices()));
		connect( (QObject*)(m_Controls.mergePancreas), SIGNAL(clicked()), (QObject*)this, SLOT(MergePancreasMap()));
		connect( (QObject*)(m_Controls.mergeKidney), SIGNAL(clicked()), (QObject*)this, SLOT(MergeKidneyMap())); 

		connect( (QObject*)(m_Controls.m_loadPolysModeling), SIGNAL(clicked()), (QObject*)this, SLOT(LoadDatasetsForModelling())); 
		connect( (QObject*)(m_Controls.selectModelType), SIGNAL(activated(int)), (QObject*)this, SLOT(SelectModelType()));
		connect( (QObject*)(m_Controls.trainingModel), SIGNAL(clicked()), (QObject*)this, SLOT(TrainShapeModel())); 

		connect( (QObject*)(m_Controls.m_segImageSelect), SIGNAL(OnSelectionChanged(const mitk::DataNode *)), (QObject*) this, SLOT(ImageSelected(const mitk::DataNode *)));
		connect( (QObject*)(m_Controls.m_segProbMapSelect), SIGNAL(OnSelectionChanged(const mitk::DataNode *)), (QObject*)this, SLOT(ProbMapSelected(const mitk::DataNode *))); 
		connect( (QObject*)(m_Controls.m_probShapeSelect), SIGNAL(OnSelectionChanged(const mitk::DataNode *)), (QObject*)this, SLOT(ProbShapeSelected(const mitk::DataNode *)));
		connect( (QObject*)(m_Controls.m_deform), SIGNAL(clicked()), (QObject*)this, SLOT(Deform())); 
		connect( (QObject*)(m_Controls.m_btnDeformModel), SIGNAL(clicked()),(QObject*) this, SLOT(Segment()));
		connect( (QObject*)(m_Controls.autoOrient), SIGNAL(clicked()), (QObject*)this, SLOT(AutoOrientationWithoutModel())); 
		connect( (QObject*)(m_Controls.autoModel), SIGNAL(clicked()), (QObject*)this, SLOT(AutoOrientationWithModel())); 

		connect( (QObject*)(m_Controls.eval_GTmaskSelect), SIGNAL(OnSelectionChanged(const mitk::DataNode *)), (QObject*)this, SLOT(SelectGroundTruthMask(const mitk::DataNode *)));
		connect( (QObject*)(m_Controls.eval_SegmentMaskSelect), SIGNAL(OnSelectionChanged(const mitk::DataNode *)), (QObject*)this, SLOT(SelectSegmentToEvaluate(const mitk::DataNode *)));
		connect( (QObject*)(m_Controls.startEvaluation), SIGNAL(clicked()), (QObject*)this, SLOT(StartSegmentationEvaluation())); 

		connect( (QObject*)(m_Controls.fixImg_TreeNode), SIGNAL(OnSelectionChanged(const mitk::DataNode *)), (QObject*)this, SLOT(FixImageSelectedForReg(const mitk::DataNode *))); 
		connect( (QObject*)(m_Controls.loadImgFolder_Registration), SIGNAL(clicked()), (QObject*)this, SLOT(LoadAndResampleImage())); 
		connect( (QObject*)(m_Controls.computeGMM), SIGNAL(clicked()), (QObject*)this, SLOT(ComputeGMM())); 

		connect( (QObject*)(m_Controls.loadLabelsForKidneys), SIGNAL(clicked()), (QObject*)this, SLOT(FilterRegionOfInterest()));
		connect( (QObject*)(m_Controls.cropByAxes), SIGNAL(clicked()), (QObject*)this, SLOT(CropImageWithDefinedAxes())); 
		connect( (QObject*)(m_Controls.generateCorrespondingMeshFRomBinaryImage), SIGNAL(clicked()), (QObject*)this, SLOT(GenerateCorrespondingMeshFromBinaryImage()));
		connect( (QObject*)(m_Controls.meshReflectionX), SIGNAL(clicked()), (QObject*)this, SLOT(PolyDataReflectionX()));
		connect( (QObject*)(m_Controls.meshReflectionY), SIGNAL(clicked()), (QObject*)this, SLOT(PolyDataReflectionY()));
		connect( (QObject*)(m_Controls.meshReflectionZ), SIGNAL(clicked()), (QObject*)this, SLOT(PolyDataReflectionZ()));
		connect( (QObject*)(m_Controls.separateMask), SIGNAL(clicked()), (QObject*)this, SLOT(SeparateMasks())); 
		connect( (QObject*)(m_Controls.maskBoundCheck), SIGNAL(clicked()), (QObject*)this, SLOT(CheckMaskBoundary())); 
		connect( (QObject*)(m_Controls.computeAndCrop), SIGNAL(clicked()), (QObject*)this, SLOT(ComputeAndCropROI())); 

		mitk::NodePredicateDimension::Pointer dimensionPredicate	= mitk::NodePredicateDimension::New(3);
		mitk::NodePredicateDataType::Pointer imagePredicate			= mitk::NodePredicateDataType::New("Image");
		mitk::NodePredicateProperty::Pointer predicateProp			= mitk::NodePredicateProperty::New("produced");

		m_Controls.m_segImageSelect->SetDataStorage(GetDataStorage());
		m_Controls.m_segImageSelect->SetPredicate(mitk::NodePredicateAnd::New(mitk::NodePredicateAnd::New(dimensionPredicate, imagePredicate), mitk::NodePredicateNot::New(predicateProp)));
		m_Controls.m_segProbMapSelect->SetDataStorage(GetDataStorage()); 
//		m_Controls.m_segProbMapSelect->SetAutoSelectNewItems(true); 
		m_Controls.m_segProbMapSelect->SetPredicate(mitk::NodePredicateAnd::New(mitk::NodePredicateAnd::New(dimensionPredicate, imagePredicate), mitk::NodePredicateNot::New(predicateProp))); 
		m_Controls.m_probShapeSelect->SetDataStorage(GetDataStorage());
		m_Controls.m_probShapeSelect->SetAutoSelectNewItems(true); 
		m_Controls.m_probShapeSelect->SetPredicate(mitk::NodePredicateDataType::New("Surface"));


		m_Controls.eval_GTmaskSelect->SetDataStorage(GetDataStorage());
		m_Controls.eval_GTmaskSelect->SetPredicate(mitk::NodePredicateAnd::New(mitk::NodePredicateAnd::New(dimensionPredicate, imagePredicate), mitk::NodePredicateNot::New(predicateProp)));
		m_Controls.eval_SegmentMaskSelect->SetDataStorage(GetDataStorage());
		m_Controls.eval_SegmentMaskSelect->SetPredicate(mitk::NodePredicateAnd::New(mitk::NodePredicateAnd::New(dimensionPredicate, imagePredicate), mitk::NodePredicateNot::New(predicateProp)));

		m_Controls.fixImg_TreeNode->SetDataStorage(GetDataStorage());
		m_Controls.fixImg_TreeNode->SetAutoSelectNewItems(true); 
		m_Controls.fixImg_TreeNode->SetPredicate(mitk::NodePredicateAnd::New(mitk::NodePredicateAnd::New(dimensionPredicate, imagePredicate), mitk::NodePredicateNot::New(predicateProp)));

		m_Controls.selectModelType->setAcceptDrops(true);
		m_Controls.selectModelType->setMaxVisibleItems(7);
		m_Controls.selectModelType->insertItem(0, "PCA");
		m_Controls.selectModelType->insertItem(1, "RPCA");
		m_Controls.selectModelType->insertItem(2, "KPCA");
		m_Controls.selectModelType->insertItem(3, "NIPS09"); 
		m_Controls.selectModelType->insertItem(4, "MICCAI17"); 
		m_Controls.selectModelType->insertItem(5, "RKPCA");
		m_Controls.selectModelType->setCurrentIndex(5);

		m_Controls.m_modelNameSelect->setEnabled(true); 
		m_Controls.m_modelNameSelect->setAcceptDrops(true);
		m_Controls.m_modelNameSelect->setMaxVisibleItems(7); 
		m_Controls.m_modelNameSelect->insertItem(0, "Kidneys"); 
		m_Controls.m_modelNameSelect->insertItem(1, "Pancreas"); 
		m_Controls.m_modelNameSelect->setCurrentIndex(1);

		m_Controls.selectSliceMode->setAcceptDrops(true); 
		m_Controls.selectSliceMode->setMaxVisibleItems(4); 
		m_Controls.selectSliceMode->insertItem(0, "SC");
		m_Controls.selectSliceMode->insertItem(1, "SA"); 
		m_Controls.selectSliceMode->insertItem(2, "CA"); 
		m_Controls.selectSliceMode->itemText(1); 

		m_Controls.closeRadius->setEnabled(true); 
		m_Controls.closeRadius->setMinimum(0); 
		m_Controls.closeRadius->setMaximum(512);
		m_Controls.closeRadius->setValue(2); 

		m_Controls.setSigma->setEnabled(true); 
		m_Controls.setSigma->setMinimum(0.001); 
		m_Controls.setSigma->setMaximum(10000); 
		m_Controls.setSigma->setValue(100); 

		m_Controls.HU_leftKidney->setValue(193); 
		m_Controls.HU_rightKidney->setValue(193); 
		m_Controls.HU_leftKidney->setMaximum(1000); 
		m_Controls.HU_rightKidney->setMaximum(1000); 
		m_Controls.HU_leftKidney->setMinimum(0);
		m_Controls.HU_rightKidney->setMinimum(0); 

		m_Controls.cropX_lower->setEnabled(true); 
		m_Controls.cropX_upper->setEnabled(true); 
		m_Controls.cropY_lower->setEnabled(true); 
		m_Controls.cropY_upper->setEnabled(true); 
		m_Controls.cropZ_lower->setEnabled(true); 
		m_Controls.cropZ_upper->setEnabled(true); 

		m_Controls.cropX_lower->setMinimum(0); 
		m_Controls.cropX_upper->setMaximum(1000); 
		m_Controls.cropY_lower->setMinimum(0);
		m_Controls.cropY_upper->setMaximum(1000); 
		m_Controls.cropZ_lower->setMinimum(0);
		m_Controls.cropZ_upper->setMaximum(1000); 

		m_Controls.resample_x->setMinimum(0); 
		m_Controls.resample_x->setMaximum(512); 
		m_Controls.resample_y->setMinimum(0); 
		m_Controls.resample_y->setMaximum(512); 
		m_Controls.resample_z->setMinimum(0); 
		m_Controls.resample_z->setMaximum(512); 

		m_Controls.check_mask->setEnabled(true); 
		m_Controls.check_mask->setChecked(false); 

		m_Controls.numVertices->setEnabled(true); 
		m_Controls.numVertices->setValue(2500); 
		m_Controls.numVertices->setMinimum(500); 
		m_Controls.numVertices->setMaximum(90000); 

		m_Controls.m_segImageSelect->setEnabled(true);
		m_Controls.m_segImageSelect->SetAutoSelectNewItems(true); 
		m_Controls.m_segProbMapSelect->setEnabled(true);
		m_Controls.m_segProbMapSelect->SetAutoSelectNewItems(true); 
		m_Controls.m_probShapeSelect->setEnabled(true);

		m_Controls.fixImg_TreeNode->setEnabled(true);
		m_Controls.eval_SegmentMaskSelect->setEnabled(true);
		m_Controls.eval_GTmaskSelect->setEnabled(true);

		m_Controls.m_cropRadius->setEnabled(true); 
		m_Controls.m_cropRadius->setValue(1); 
		m_Controls.m_cropRadius->setMinimum(0); 
		m_Controls.m_cropRadius->setMaximum(50); 

		m_Controls.m_lambdaNN->setEnabled(true); 
		m_Controls.m_lambdaNN->setValue(1.0); 
		m_Controls.m_lambdaNN->setMinimum(0); 
		m_Controls.m_lambdaNN->setMaximum(10.0); 

		m_Controls.m_iterNum->setEnabled(true);
		m_Controls.m_iterNum->setValue(5.00);
		m_Controls.m_iterNum->setMinimum(0.0);
		m_Controls.m_iterNum->setMaximum(1000.0);

		m_Controls.m_checkAutoOrient->setEnabled(true); 
		m_Controls.m_checkAutoOrient->setChecked(true); 

		m_Controls.m_scalingModel->setEnabled(true); 
		m_Controls.m_scalingModel->setChecked(false); 
	}


	// load statemachine for user interaction
	QDir exePath = QDir(qApp->applicationDirPath() + "/ModelBasedOrganSegmentationStateMachines.xml");
	if( ! mitk::GlobalInteraction::GetInstance()->GetStateMachineFactory()->LoadBehavior( exePath.absolutePath().toStdString() ) ) 
	{
		mitk::GlobalInteraction::GetInstance()->GetStateMachineFactory()->LoadBehavior("J:/MITK_SOURCE/MITK-2014.10.0/Core/Code/Resources/ModelBasedOrganSegmentationStateMachines.xml"); 
	}
}


///////////////////////// --- Modelling Procedure : m_shapeModel ---- //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void SegmentationView::LoadDatasetsForModelling()
{
	if (this->m_trainingDatasets.size() > 0)
	{
		while (!m_trainingDatasets.empty())
		{
			m_trainingDatasets.pop_back();
		}
	}

	QFileDialog fd;
	fd.setDirectory(m_recentDir.c_str());
	fd.setFileMode(QFileDialog::Directory);
	QString selected;

	if (fd.exec() == QDialog::Accepted)
	{
		m_recentDir = fd.directory().absolutePath().toAscii();
		QStringList myfiles = fd.selectedFiles();
		if (!myfiles.isEmpty())
			selected = myfiles[0];
	}

	QDir* meshDir = new QDir(selected);
	QStringList::Iterator it;
	QStringList files = meshDir->entryList();
	it = files.begin();

	while (it != files.end())
	{
		if (QFileInfo(*meshDir, *it).isFile() && ((*it).endsWith)(".vtk"))
		{
			vtkSmartPointer< vtkPolyDataReader > pReader = vtkSmartPointer<vtkPolyDataReader>::New();

			pReader->SetFileName(QFileInfo(*meshDir, *it).absoluteFilePath().toAscii());
			pReader->Update();

			vtkSmartPointer< vtkPolyData > loadedPoly = pReader->GetOutput();
			this->m_trainingDatasets.push_back(loadedPoly);
		}
		++it;
	}

	delete meshDir;

	m_PCAModel = new RobustKPCA(m_Controls.setSigma->value(), m_Controls.m_scalingModel->isChecked());
	m_PCAModel->ReadDataMatrix(this->m_trainingDatasets);
	m_PCAModel->performPCA();

	std::cout <<  __FUNCTION__ << " : " << this->m_trainingDatasets.size() << " meshes are loaded for model training." << std::endl;

}


void SegmentationView::SelectModelType()
{
	this->m_modelType = m_Controls.selectModelType->currentText().toStdString(); 
	std::cout << __FUNCTION__ << " : current selected model type : " << this->m_modelType << std::endl;
}


void SegmentationView::TrainShapeModel()
{
	if (!this->m_deformedShape) this->m_deformedShape = vtkPolyData::New();

	if (this->m_trainingDatasets.size() == 0)
	{
		MITK_WARN << " Please load datasets ahead! " << std::endl;
		return;
	}

	std::vector< vtkSmartPointer< vtkPolyData > > reconstructedDatasets;
	std::vector< vtkSmartPointer< vtkPolyData > > pcaDatasets;
	for (int i = 0; i < this->m_trainingDatasets.size(); i++)
	{
		vtkSmartPointer< vtkPolyData > poly = vtkSmartPointer< vtkPolyData >::New();
		poly->DeepCopy(this->m_trainingDatasets.at(i));
		reconstructedDatasets.push_back(poly);
		pcaDatasets.push_back(poly);
	}

	/********************************************************************************************************/
	itk::TimeProbe itkClock;
	itkClock.Start();

	if (this->m_shapeModel)
	{
		std::cout << __FUNCTION__ << " Deleting the current model and create a new model ... " << std::endl;
		delete this->m_shapeModel;
	}
	m_shapeModel = new RobustKPCA(m_Controls.setSigma->value(), this->m_Controls.m_scalingModel->isChecked());
	m_shapeModel->ReadDataMatrix(this->m_trainingDatasets);


	if (this->m_modelType == "PCA")
	{
		m_shapeModel->performPCA();
		reconstructedDatasets = m_shapeModel->getConstructed_PCA();
	}
	else if (this->m_modelType == "RPCA")
	{
		m_shapeModel->RobustPCA();
		reconstructedDatasets = m_shapeModel->getConstructed_RobustPCA();
	}
	else if (this->m_modelType == "KPCA")
	{
		m_shapeModel->performKPCA();
		reconstructedDatasets = m_shapeModel->getConstructed_KPCA();
	}
	else if (this->m_modelType == "RKLRR")
	{
		m_shapeModel->RobustKernelLRR(); 
		reconstructedDatasets = m_shapeModel->getConstructed_RobustKernelLRR(); 
	}
	else if (this->m_modelType == "NIPS09")
	{
		m_shapeModel->performNIPS09(); 
		reconstructedDatasets = m_shapeModel->getConstructed_NIPS09(); 
	}
	else
	{
		m_shapeModel->RobustKernelPCA();
		reconstructedDatasets = m_shapeModel->getConstructed_RobustKernelPCA();
	}

	this->m_deformedShape->DeepCopy(this->m_shapeModel->GetMeanShape());

	if (!this->m_deformedShape->GetPointData()->GetArray("weights"))
	{
		// set new weights 
		vtkSmartPointer< vtkDataArray > weights = vtkSmartPointer< vtkDoubleArray >::New();
		weights->SetName("weights");
		weights->SetNumberOfTuples(this->m_deformedShape->GetNumberOfPoints());
		for (int i = 0; i < weights->GetNumberOfTuples(); i++) weights->SetTuple1(i, 1.0);
		weights->Modified();
		this->m_deformedShape->GetPointData()->SetScalars(weights);
		this->m_deformedShape->Modified();
	}

	// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	itkClock.Stop();
	MITK_INFO << __FUNCTION__ << " : Time for Training ( " << reconstructedDatasets.size() << " datasets ): " << itkClock.GetMean() << "seconds ." << std::endl;
	std::cout << std::endl;

}


//////////////////////// ---- Segmentation Tools ---- ////////////////////////////////////////////////////////////////////////////////////////////////


void SegmentationView::ImageSelected(const mitk::DataNode * node_)
{
	if (node_)
	{
		std::string segImgName; 
		if (node_->GetName(segImgName))
		{
			this->m_segImageName = segImgName;
			std::cout << __FUNCTION__ << " : Current Selected Image to Segment : " << m_segImageName ;
			
			mitk::BaseData* data = node_->GetData();
			if (data)
			{
				mitk::Image::Pointer mitkImage = dynamic_cast<mitk::Image*>(data);
				if (mitkImage)
				{
					mitk::CastToItkImage(mitkImage, this->m_segImage);
				}
			}

			if (m_Controls.m_modelNameSelect->currentText() == "Pancreas")
			{
				// find the shape order , remove the substring "Pancrea" 
				std::string subString = "Pancrea";
				std::string::size_type foundPos = segImgName.find(subString);

				while (foundPos != std::string::npos)
				{
					segImgName.erase(foundPos, subString.length());
					foundPos = segImgName.find(subString, foundPos);
				}

				this->m_shapeOrder = 0;
				std::stringstream geek(segImgName);

				geek >> this->m_shapeOrder;

				if (this->m_shapeOrder > 5) this->m_shapeOrder = this->m_shapeOrder - 1;
				else if (this->m_shapeOrder == 5) this->m_shapeOrder = 0;
			}


			if (m_Controls.m_modelNameSelect->currentText() == "Kidneys")
			{
				this->m_shapeOrder = 0;
				std::stringstream geek(segImgName);

				geek >> this->m_shapeOrder;
			}

			cout << "  with m_shapeOrder = " << m_shapeOrder << " of " << m_Controls.m_modelNameSelect->currentText().toStdString() << endl;

		}
	}
}


void SegmentationView::ProbMapSelected(const mitk::DataNode * probNode)
{
	if (probNode)
	{
		mitk::BaseData* data = probNode->GetData();
		this->m_probMapName = probNode->GetName(); 
		
		if (data)
		{
			mitk::Image* image = dynamic_cast<mitk::Image*>(data);
			if (image)
			{
				std::cout << __FUNCTION__ << " : the loaded probability map is " << this->m_probMapName << std::endl;
				mitk::CastToItkImage(image, this->m_probMap); 
			}
		}
	}
	 
}


void SegmentationView::ProbShapeSelected(const mitk::DataNode * shapeNode)
{
	/*if (!this->m_probShape) this->m_probShape = vtkSmartPointer< vtkPolyData >::New();

	if (shapeNode)
	{
		mitk::BaseData* data = shapeNode->GetData();
		
		if (data)
		{
			mitk::Surface::Pointer surface = dynamic_cast< mitk::Surface* >(data); 
			if (surface)
			{
				std::cout << __FUNCTION__ << " m_probShape is selected " << std::endl;
				this->m_probShape->DeepCopy(surface->GetVtkPolyData());
			}
		}
	}*/
}


void SegmentationView::Deform()  /* pass m_deformedShape to m_surfaceNode */
{
	this->loadProbShapeForLocalization();   // return m_probShape

	if (!this->m_shapeModel)                // return m_deformShape , not in use if model is created 
	{
		this->loadDeformShapeFromLocal();   
	}

	/**> initialize weights <**/
	vtkSmartPointer< vtkDataArray > weights = vtkSmartPointer< vtkDoubleArray >::New();
	weights->SetName("weights");
	weights->SetNumberOfTuples(this->m_deformedShape->GetNumberOfPoints());
	for (int i = 0; i < weights->GetNumberOfTuples(); i++) 
		weights->SetTuple1(i, 0.0); 
	weights->Modified();
	this->m_deformedShape->GetPointData()->SetScalars(weights); 
	this->m_deformedShape->Modified(); 
	
	this->alignDeformShapeWithProbShape(this->m_deformedShape, this->m_probShape);   // align deformation shape (mean by default) with probability map with ICP 
	
	this->transformDeformShapeToProbMap(this->m_deformedShape);   // transform with closest Cell locator 

	if (this->m_shapeModel)
	{
		for (int i = 0; i < 2; i++)
		{
			if (this->m_modelType == "PCA" || this->m_modelType == "RPCA")
			{
				m_deformedShape->DeepCopy(this->m_shapeModel->ProjectShapeLinearModel(m_deformedShape));
			}
			else if (this->m_modelType == "KPCA" || this->m_modelType == "RKPCA")
			{
				m_deformedShape->DeepCopy(this->m_shapeModel->ProjectShapeNonlinearModel(this->m_deformedShape)); 
			}
		}

		if (this->m_modelType == "PCA" || this->m_modelType == "RPCA")
		{
			m_deformedShape->DeepCopy(this->m_shapeModel->ProjectShapeLinearModel(m_deformedShape));
		}
		else if (this->m_modelType == "KPCA" || this->m_modelType == "RKPCA")
		{
			m_deformedShape->DeepCopy(this->m_shapeModel->ProjectShapeNonlinearModel(this->m_deformedShape));
//			m_deformedShape->DeepCopy(this->m_shapeModel->ProjectShapeNonlinearModelWithWeights(m_deformedShape, this->m_shapeOrder, true));
		}
	}


/// >>>>>>>>>>>>>>>>>>>>> create m_imageAnalyser && Localization >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	mitk::DataNode::Pointer nodeImg = m_Controls.m_segImageSelect->GetSelectedNode();
	mitk::BaseData* data = nodeImg->GetData();
	mitk::Image* image = dynamic_cast<mitk::Image*>(data);

	m_imageAnalyser = new ImageAnalyser(nodeImg, this->m_probMap);

//	this->loadGroundTruthShapeForPreTraining();   // later use m_deformedShape || m_probShape 

	m_imageAnalyser->trainPriorGMM(this->m_probShape);  // m_probImage -> m_currentPoly iteratively this->m_preTrainShape

	m_imageAnalyser->computePriorProbs(this->m_probShape, 0.0, 3);   // m_probImage -> m_currentPoly iteratively this->m_preTrainShape

//>>>>>>>>>>>>>>>>>>>>>>>> Manual - deformation iteratively >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	mitk::RenderingManager::GetInstance()->RequestUpdateAll();
	mitk::LevelWindowProperty::Pointer levWinProp = mitk::LevelWindowProperty::New();
	mitk::LevelWindow levelWindow;
	levelWindow.SetLevelWindow(50, 350);
	levWinProp->SetLevelWindow(levelWindow);
	m_imageAnalyser->GetImageNode()->GetPropertyList()->SetProperty("levelwindow", levWinProp);
	
	mitk::Surface::Pointer surface = mitk::Surface::New();
	surface->SetVtkPolyData(this->m_deformedShape);
	surface->CalculateBoundingBox();

	m_currentSurfaceNode = mitk::DataNode::New();
	m_currentSurfaceNode->SetData(surface);
	m_currentSurfaceNode->SetColor(0.0, 1.0, 1.0);
	m_currentSurfaceNode->SetName(m_Controls.m_modelNameSelect->currentText().toStdString().c_str());
	m_currentSurfaceNode->SetProperty("LookupTable", m_lookupTableProp);
	m_currentSurfaceNode->GetPropertyList()->SetProperty("scalar visibility", mitk::BoolProperty::New(true));
	m_currentSurfaceNode->GetPropertyList()->SetProperty("ScalarsRangeMinimum", mitk::FloatProperty::New(0.0));
	m_currentSurfaceNode->GetPropertyList()->SetProperty("ScalarsRangeMaximum", mitk::FloatProperty::New(100.0));
	mitk::Material* mp2 = dynamic_cast<mitk::Material*>(m_currentSurfaceNode->GetProperty("material"));
	if (mp2)
	{
		mp2->SetColorCoefficient(0.001);
		mp2->SetRepresentation(mitk::Material::Wireframe);
	}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	m_AffineSurfaceInteractor = AffineSurfaceInteractor::New("AffineInteractions ctrl-drag", m_currentSurfaceNode);
	mitk::GlobalInteraction::GetInstance()->AddInteractor(m_AffineSurfaceInteractor);

	m_currentSurfaceNode->Modified();
	m_currentSurfaceNode->SetColor(0.0, 1.0, 1.0);

	GetDataStorage()->Add(m_currentSurfaceNode, m_imageAnalyser->GetImageNode());
	mitk::RenderingManager::GetInstance()->RequestUpdateAll();

	std::cout << __FUNCTION__ << " all procedure has finished, start segmentation ! " << std::endl;

	return;

} 


void SegmentationView::AutoOrientationWithoutModel()
{
	std::cout << std::endl;
	std::cout << __FUNCTION__ << " Start ... " << std::endl;
	std::cout << std::endl;

	this->m_currentSurfaceNode->Update();

	dynamic_cast<mitk::Surface*>(this->m_currentSurfaceNode->GetData())->CalculateBoundingBox();
	vtkSmartPointer<vtkPolyData> currentModelPoly = vtkSmartPointer< vtkPolyData>::New();
	currentModelPoly->DeepCopy(dynamic_cast<mitk::Surface*>(this->m_currentSurfaceNode->GetData())->GetVtkPolyData());


	if (m_icpSurface.IsNull())  m_icpSurface = mitk::Surface::New();
	if (!m_icpPoly) m_icpPoly = vtkSmartPointer<vtkPolyData>::New();
	else this->ResetPolyData(m_icpPoly);
	m_icpPoly->DeepCopy(currentModelPoly);
	m_icpSurface->SetVtkPolyData(m_icpPoly);
	m_icpSurface->CalculateBoundingBox();


	if (m_structureSurface.IsNull())  m_structureSurface = mitk::Surface::New();
	if (!m_structurePoly) m_structurePoly = vtkSmartPointer<vtkPolyData>::New();
	else this->ResetPolyData(m_structurePoly);
	m_structurePoly->DeepCopy(currentModelPoly);
	m_structureSurface->SetVtkPolyData(m_structurePoly);
	m_structureSurface->CalculateBoundingBox();

	if (m_Controls.m_checkAutoOrient->isChecked() == false)
	{
		this->AutoSearchBasedOnIntensModel(m_Controls.m_iterNum->value(), m_Controls.m_lambdaNN->value(), m_Controls.m_cropRadius->value());
	}
	else
	{
		// model deformation  
		if (this->m_shapeModel)
		{
			this->AutoOrientationWithModel();
		}

		cout << endl;
		cout << __FUNCTION__ << " -------- Use Probability Map with Radius = 3, iter = 50 ------------ " << endl;
		this->AutoSearchBasedOnIntensModel(50, 1.0, 3);

		cout << endl;
		cout << __FUNCTION__ << " -------- Use Probability Map with Radius = 2, iter = 50 ------------- " << endl;
		this->AutoSearchBasedOnIntensModel(50, 1.0, 2);

		cout << endl;
		cout << __FUNCTION__ << " -------- Use Probability Map with Radius = 1, iter = 100 ------------ " << endl;
		this->AutoSearchBasedOnIntensModel(50, 1.0, 1);

		cout << endl;
		cout << __FUNCTION__ << " -------- Use Probability Map with Radius = 0, iter = 100 ------------- " << endl;
		this->AutoSearchBasedOnIntensModel(100, 1.0, 0);

		cout << endl;
		cout << __FUNCTION__ << " --------------- Not Use Probability Map with Radius = 0, iter = 100 ----------------------- " << endl;
		this->AutoSearchBasedOnIntensModel(100, 0.0, 0);
	}

	std::cout << " [ " << __FUNCTION__ << " ] Done ! " << std::endl;

	return;
}


void SegmentationView::AutoOrientationWithModel()
{
	std::cout << std::endl;
	std::cout << __FUNCTION__ << " Start ... " << std::endl;
	std::cout << std::endl;

	this->m_currentSurfaceNode->Update();

	dynamic_cast<mitk::Surface*>(this->m_currentSurfaceNode->GetData())->CalculateBoundingBox();
	vtkSmartPointer<vtkPolyData> currentModelPoly = vtkSmartPointer< vtkPolyData>::New();
	currentModelPoly->DeepCopy(dynamic_cast<mitk::Surface*>(this->m_currentSurfaceNode->GetData())->GetVtkPolyData());

	if (m_icpSurface.IsNull())  m_icpSurface = mitk::Surface::New();
	if (!m_icpPoly) m_icpPoly = vtkSmartPointer<vtkPolyData>::New();
	else this->ResetPolyData(m_icpPoly);

	m_icpPoly->DeepCopy(currentModelPoly);
	m_icpSurface->SetVtkPolyData(m_icpPoly);
	m_icpSurface->CalculateBoundingBox();

	if (m_structureSurface.IsNull())  m_structureSurface = mitk::Surface::New();
	if (!m_structurePoly) m_structurePoly = vtkSmartPointer<vtkPolyData>::New();
	else this->ResetPolyData(m_structurePoly);

	m_structurePoly->DeepCopy(currentModelPoly);
	m_structureSurface->SetVtkPolyData(m_structurePoly);
	m_structureSurface->CalculateBoundingBox();

	//>>>>>>>>>>>>>>>>>>>>>>>>

	this->AutoRefineWithModel(3, false);

	//<<<<<<<<<<<<<<<<<<<<<<<

	std::cout << __FUNCTION__ << " Done ! " << std::endl;

	return;
}


void SegmentationView::Segment()
{	
	std::cout << std::endl; 
	std::cout << "--------------- Start Segmentation : Segment -------------" << std::endl; 
	std::cout << std::endl; 

//>>>>>>>>>>>>>>>>>>>>>>>>> AddSegmentationToDataStorage();

	if (!m_currentSurfaceNode) return;
	if (!m_currentSurfaceNode->GetData()) return;

	std::string imageName;
	if (m_imageAnalyser->GetImageNode()->GetName(imageName))
	{
		imageName += "_" + this->m_modelType + "_Segment";
		mitk::DataNode::Pointer segNode = GetDataStorage()->GetNamedNode(imageName);
		if (segNode.IsNotNull())
		{
			GetDataStorage()->Remove(segNode);
			segNode = NULL;
		}
	}

	mitk::Image::Pointer mitkImage = m_imageAnalyser->GetMitkImage();
	mitk::Surface::Pointer surface;

	surface = dynamic_cast<mitk::Surface*>(m_currentSurfaceNode->GetData());

	mitk::SurfaceToImageFilter::Pointer mitkFilter = mitk::SurfaceToImageFilter::New();
	mitkFilter->SetInput(surface);
	mitkFilter->SetImage(mitkImage);
	mitkFilter->SetMakeOutputBinary(true);	// save as binary or not? 

	try
	{
		mitkFilter->Update();
	}
	catch (...)
	{
		QMessageBox::information(NULL, "SurfaceToImageFilter fail", QString("Could not convert segmentation."), QMessageBox::Ok); return;
	}

	mitk::Image::Pointer outputImage = mitkFilter->GetOutput();

	mitk::ImageTimeSelector::Pointer timeSelector = mitk::ImageTimeSelector::New();
	timeSelector->SetInput(outputImage);
	timeSelector->SetTimeNr(0);
	timeSelector->UpdateLargestPossibleRegion();

	outputImage = timeSelector->GetOutput();

	// add to data storage
	mitk::DataNode::Pointer newNode = mitk::DataNode::New();
	newNode->SetData(outputImage);

	AtlasImageType::Pointer outputItkImg = AtlasImageType::New();
	mitk::CastToItkImage(outputImage, outputItkImg);
	this->addImageToNode(outputItkImg, m_probMapName + m_modelType + "_segment", true);

	mitk::RenderingManager::GetInstance()->RequestUpdateAll();
		
	std::cout << " [ " << __FUNCTION__ << " ] Kidney Segmentation Done ! " << std::endl; 

	return;
}


/////////////////////// ---- Internal toolkits ---- ///////////////////////////////////////////////


void SegmentationView::AutoSearchBasedOnIntensModel(int iterationCount, double _lambda, int _radius)
{
	this->m_currentSurfaceNode->Update(); 

	for (int n = 0; n < iterationCount; n++)
	{
		mitk::Surface::Pointer currentSurface = dynamic_cast<mitk::Surface*>(m_currentSurfaceNode->GetData());

		currentSurface = this->GetTransformedGeometrySurface(currentSurface);  // vtkLinearTransform of Surface

		m_currentSurfaceNode->SetData(currentSurface);

		//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Add Color >>>>>>>>>>>>>>>>>>>>>>>>>

		m_currentSurfaceNode->SetProperty("LookupTable", m_lookupTableProp);
		m_currentSurfaceNode->GetPropertyList()->SetProperty("scalar visibility", mitk::BoolProperty::New(true));
		m_currentSurfaceNode->GetPropertyList()->SetProperty("ScalarsRangeMinimum", mitk::FloatProperty::New(0.0));
		m_currentSurfaceNode->GetPropertyList()->SetProperty("ScalarsRangeMaximum", mitk::FloatProperty::New(100.0));
		mitk::Material* mp2 = dynamic_cast<mitk::Material*>(m_currentSurfaceNode->GetProperty("material"));
		if (mp2)
		{
			mp2->SetColorCoefficient(0.001);
			mp2->SetRepresentation(mitk::Material::Wireframe);
		}

		m_currentSurfaceNode->Update();
	
		//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

		vtkSmartPointer<vtkPolyData> currentPoly = currentSurface->GetVtkPolyData();

		this->ResetPolyData(m_structurePoly);
		m_structurePoly->DeepCopy(currentPoly); // structure poly needs constraints, since it is used as a container for the SSM
		m_structurePoly->Modified();

		this->reComputeNormals(currentPoly);

		int onBoundary = 0; 
		m_imageAnalyser->landmarkDeformation(currentPoly, m_structurePoly, _lambda, _radius, onBoundary);
		
		if (onBoundary >= currentPoly->GetNumberOfPoints() - 100)
		{
			cout << __FUNCTION__ << " enough ! " << endl; 
			break; 
		}

		if (this->m_PCAModel && n % 7 == 0 && m_Controls.m_checkAutoOrient->isChecked())
		{
			m_structurePoly->DeepCopy(this->m_PCAModel->ProjectShapeLinearModel(m_structurePoly));
		}

		this->ResetPolyData(currentPoly);
		currentPoly->DeepCopy(m_structurePoly); 
		currentPoly->Modified();

		mitk::RenderingManager::GetInstance()->RequestUpdateAll();
		mitk::ProgressBar::GetInstance()->Progress();

		currentSurface->CalculateBoundingBox(); // poly data has changed -> update surface
		currentSurface->Update();

		this->ResetPolyData(m_icpPoly);
		m_icpPoly->DeepCopy(currentPoly);
		m_icpPoly->Modified();
		m_icpSurface->SetVtkPolyData(m_icpPoly);
		m_icpSurface->CalculateBoundingBox();

		mitk::RenderingManager::GetInstance()->RequestUpdateAll();
		qApp->processEvents();

	}

}


void SegmentationView::AutoRefineWithModel(int iterationCount, bool useWeights)
{
	m_currentSurfaceNode->Modified();

	mitk::Surface::Pointer currentSurface = dynamic_cast<mitk::Surface*>(m_currentSurfaceNode->GetData());

	currentSurface = this->GetTransformedGeometrySurface(currentSurface);

	m_currentSurfaceNode->SetData(currentSurface);

	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	m_currentSurfaceNode->SetProperty("LookupTable", m_lookupTableProp);
	m_currentSurfaceNode->GetPropertyList()->SetProperty("scalar visibility", mitk::BoolProperty::New(true));
	m_currentSurfaceNode->GetPropertyList()->SetProperty("ScalarsRangeMinimum", mitk::FloatProperty::New(0.0));
	m_currentSurfaceNode->GetPropertyList()->SetProperty("ScalarsRangeMaximum", mitk::FloatProperty::New(100.0));
	mitk::Material* mp2 = dynamic_cast<mitk::Material*>(m_currentSurfaceNode->GetProperty("material"));
	if (mp2)
	{
		mp2->SetColorCoefficient(0.001);
		mp2->SetRepresentation(mitk::Material::Wireframe);
	}

	m_currentSurfaceNode->Update();

	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	vtkSmartPointer<vtkPolyData> currentPoly = currentSurface->GetVtkPolyData();

	this->ResetPolyData(m_structurePoly);
	m_structurePoly->DeepCopy(currentPoly); // structure poly needs constraints, since it is used as a container for the SSM
	m_structurePoly->Modified();

	this->reComputeNormals(currentPoly); 
	this->ResetPolyData(currentPoly);

	
	for (int n = 0; n < iterationCount; n++)
	{			
		if (this->m_modelType == "PCA" || this->m_modelType == "RPCA")
		{
			currentPoly->DeepCopy(this->m_shapeModel->ProjectShapeLinearModel(m_structurePoly));
		}
		else if (this->m_modelType == "KPCA" || this->m_modelType == "RKPCA")
		{
			currentPoly->DeepCopy(this->m_shapeModel->ProjectShapeNonlinearModel(m_structurePoly));
		}

		currentPoly->Modified();

		this->alignDeformShapeWithProbShape(currentPoly, m_probShape);  

		this->reComputeNormals(currentPoly); 

		mitk::RenderingManager::GetInstance()->RequestUpdateAll();
		mitk::ProgressBar::GetInstance()->Progress();
	}

	currentSurface->CalculateBoundingBox(); // poly data has changed -> update surface
	currentSurface->Update();

	this->ResetPolyData(m_icpPoly);
	m_icpPoly->DeepCopy(currentPoly);
	m_icpPoly->Modified();
	m_icpSurface->SetVtkPolyData(m_icpPoly);
	m_icpSurface->CalculateBoundingBox();

	mitk::RenderingManager::GetInstance()->RequestUpdateAll();
	qApp->processEvents();

}


//////////////////////// ---- Evaluation Tools ---- ////////////////////////////////////////////////////////////////////////////////////////////////


void SegmentationView::SelectGroundTruthMask(const mitk::DataNode* _GTNode)
{
	if (_GTNode){
		std::string name;
		if (_GTNode->GetName(name)){
			mitk::BaseData* data = _GTNode->GetData();
			if (data){
				mitk::Image* image = dynamic_cast<mitk::Image*>(data);
				if (image){
					mitk::CastToItkImage(image, m_GTMask); 
				}
			}
		}
	}
}


void SegmentationView::SelectSegmentToEvaluate(const mitk::DataNode* _segmentNode)
{
	if (_segmentNode){
		std::string name;
		if (_segmentNode->GetName(name)){
			mitk::BaseData* data = _segmentNode->GetData();
			if (data){
				mitk::Image* image = dynamic_cast<mitk::Image*>(data);
				if (image){
					mitk::CastToItkImage(image, m_SegmentMask);
				}
			}
		}
	}
}


void SegmentationView::StartSegmentationEvaluation()
{

	typedef itk::LabelOverlapMeasuresImageFilter< AtlasImageType > MeasureType; 
	MeasureType::Pointer measure = MeasureType::New(); 


	measure->SetSourceImage(this->m_SegmentMask); 
	measure->SetTargetImage(this->m_GTMask); 
	measure->Update(); 

	std::cout << "  GetJaccardCoefficient : " << measure->GetJaccardCoefficient() << std::endl;
	std::cout << "  GetDiceCoefficient : " << measure->GetDiceCoefficient() << std::endl;


	std::cout << "                                          "
		<< "************ All Labels *************" << std::endl;
	std::cout << std::setw(10) << "   "
		<< std::setw(17) << "Total"
		<< std::setw(17) << "Union (jaccard)"
		<< std::setw(17) << "Mean (dice)"
		<< std::setw(17) << "Volume Similarity."
		<< std::setw(17) << "False negative"
		<< std::setw(17) << "False positive" << std::endl;
	std::cout << std::setw(10) << "   ";
	std::cout << std::setw(17) << measure->GetTotalOverlap();
	std::cout << std::setw(17) << measure->GetUnionOverlap();
	std::cout << std::setw(17) << measure->GetMeanOverlap();
	std::cout << std::setw(17) << measure->GetVolumeSimilarity();
	std::cout << std::setw(17) << measure->GetFalseNegativeError();
	std::cout << std::setw(17) << measure->GetFalsePositiveError();
	std::cout << std::endl;

	/*std::cout << "                                       "
		<< "************ Individual Labels *************" << std::endl;
	std::cout << std::setw(10) << "Label"
		<< std::setw(17) << "Target"
		<< std::setw(17) << "Union (jaccard)"
		<< std::setw(17) << "Mean (dice)"
		<< std::setw(17) << "Volume similarity."
		<< std::setw(17) << "False negative"
		<< std::setw(17) << "False positive" << std::endl;*/
	//MeasureType::MapType labelMap = measure->GetLabelSetMeasures();
	//MeasureType::MapType::const_iterator it;
	//for (it = labelMap.begin(); it != labelMap.end(); ++it)
	//{
	//	if ((*it).first == 0)
	//	{
	//		continue;
	//	}
	//	int label = (*it).first;
	//	std::cout << std::setw(10) << label;
	//	std::cout << std::setw(17) << measure->GetTargetOverlap(label);
	//	std::cout << std::setw(17) << measure->GetUnionOverlap(label);
	//	std::cout << std::setw(17) << measure->GetMeanOverlap(label);
	//	std::cout << std::setw(17) << measure->GetVolumeSimilarity(label);
	//	std::cout << std::setw(17) << measure->GetFalseNegativeError(label);
	//	std::cout << std::setw(17) << measure->GetFalsePositiveError(label);
	//	std::cout << std::endl;
	//}

	/*typedef itk::HausdorffDistanceImageFilter< AtlasImageType, AtlasImageType > HausdorffType; 
	HausdorffType::Pointer hausdorff = HausdorffType::New(); 
	hausdorff->SetInput1(this->m_GTMask); 
	hausdorff->SetInput2(this->m_SegmentMask); 
	hausdorff->Update(); 

	std::cout << "  Hausdorff Distance : " << std::setw(17) << hausdorff->GetHausdorffDistance() << std::endl;*/
	//std::cout << "  Average Hausdorff Distance : " << std::setw(17) << hausdorff->GetAverageHausdorffDistance() << std::endl;

	std::cout << std::endl; 
	std::cout << " ---------------------------------------------------------------------------------------- " << std::endl; 
	std::cout << std::endl; 
}


///////////////////////// ---- Deep Learning Training Utilities ---- //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void SegmentationView::LoadImagesForSliceGeneration()
{
	// load image folder 
	QFileDialog fd;
	fd.setDirectory(m_recentDir.c_str());
	fd.setFileMode(QFileDialog::Directory);
	QString selected;

	if (fd.exec() == QDialog::Accepted)
	{
		m_recentDir = fd.directory().absolutePath().toAscii();

		QStringList myfiles = fd.selectedFiles();

		if (!myfiles.isEmpty())
			selected = myfiles[0];
	}

	QDir* imageDir = new QDir(selected);
	QStringList::Iterator it;
	QStringList files = imageDir->entryList();

	it = files.begin();
	int Cnt = 0;

	while (it != files.end())
	{
		if (QFileInfo(*imageDir, *it).isFile() && ((*it).endsWith)(".pic") || (*it).endsWith(".nrrd") || (*it).endsWith(".mhd") || (*it).endsWith(".nii.gz"))
		{
			Cnt++;

			mitk::DataNodeFactory::Pointer nodeReader = mitk::DataNodeFactory::New();

			QFileInfo& fileInfoIMG_ = QFileInfo(*imageDir, *it);
			std::string imgFileName = fileInfoIMG_.absoluteFilePath().toStdString();
			std::string imgAbsName = QFileInfo(*imageDir, *it).baseName().toStdString();

			Image3DType::Pointer itkImage = Image3DType::New();

			try
			{
				nodeReader->SetFileName(imgFileName);
				nodeReader->Update();

				mitk::DataNode::Pointer imageNode = nodeReader->GetOutput();
				mitk::Image::Pointer mitkImage = dynamic_cast<mitk::Image*>(imageNode->GetData());
				mitk::CastToItkImage(mitkImage, itkImage);
			}
			catch (itk::ExceptionObject & err)
			{
				std::cerr << " Exception Object Caught in reading nrrd image !" << err << std::endl;
				return;
			}

			if (m_Controls.check_mask->isChecked())
			{
				typedef itk::RescaleIntensityImageFilter< Image3DType, Image3DType > RescalerType;
				RescalerType::Pointer rescaler = RescalerType::New();
				rescaler->SetInput(itkImage);
				rescaler->SetOutputMaximum(255);
				rescaler->SetOutputMinimum(0);
				rescaler->Update();

				itkImage = rescaler->GetOutput();
			}

			this->generateSlicesWithMode(imgAbsName, m_Controls.selectSliceMode->currentText().toStdString(), itkImage);

		}
		++it;
	}

	delete imageDir;

	std::cout << " ----------- Finish Generating Image Slices --------------------- " << std::endl;
}


void SegmentationView::LoadAndMergeSASlices()
{
	std::cout << " [ " << __FUNCTION__ << " ] Merge SA Slices .... " << std::endl;

	QFileDialog fd;
	fd.setDirectory(m_recentDir.c_str());
	fd.setFileMode(QFileDialog::DirectoryOnly);
	fd.setOption(QFileDialog::DontUseNativeDialog, true);

	QListView *list = fd.findChild< QListView* >("listView");
	if (list)  list->setSelectionMode(QAbstractItemView::MultiSelection);
	QTreeView *tree = fd.findChild< QTreeView* >();
	if (tree)  tree->setSelectionMode(QAbstractItemView::MultiSelection);

	if (fd.exec() == QDialog::Accepted)
	{
		m_recentDir = fd.directory().absolutePath().toAscii();
		QStringList myfiles = fd.selectedFiles();

		QString selected;

		if (!myfiles.isEmpty() && myfiles.size() > 1)
		{
			for (unsigned int dirCount = 1; dirCount <= myfiles.size() - 1; dirCount++)
			{
				selected = myfiles[dirCount];  // the first folder is empty 

				AtlasImageType::Pointer SAItkImage = AtlasImageType::New();
				AtlasImageType::RegionType region;
				AtlasImageType::SizeType size = { 256, 256, 256 };
				AtlasImageType::IndexType index = { 0, 0, 0 };
				region.SetIndex(index);
				region.SetSize(size);
				SAItkImage->SetRegions(region);
				SAItkImage->Allocate();
				SAItkImage->FillBuffer(itk::NumericTraits< AtlasImageType::PixelType >::Zero);

				QDir* imageDir = new QDir(selected);
				QStringList::Iterator it;
				QStringList files = imageDir->entryList();

				it = files.begin();
				int sliceNum = 0;

				std::stringstream ss; ss << dirCount;
				std::string absName = "SA_" + ss.str();

				while (it != files.end())
				{
					if (QFileInfo(*imageDir, *it).isFile() && ((*it).endsWith)(".jpg"))
					{
						mitk::DataNodeFactory::Pointer nodeReader = mitk::DataNodeFactory::New();

						QFileInfo& fileInfoIMG_ = QFileInfo(*imageDir, *it);
						std::string imgFileName = fileInfoIMG_.absoluteFilePath().toStdString();

						Image2DType::Pointer jpgImg = Image2DType::New();

						try
						{
							nodeReader->SetFileName(imgFileName);
							nodeReader->Update();

							mitk::DataNode::Pointer imageNode = nodeReader->GetOutput();
							mitk::Image::Pointer mitkImage = dynamic_cast<mitk::Image*>(imageNode->GetData());
							mitk::CastToItkImage(mitkImage, jpgImg);

							// copy the intensity to itk image 
							int xDim = jpgImg->GetLargestPossibleRegion().GetSize()[0];
							int zDim = jpgImg->GetLargestPossibleRegion().GetSize()[1];

							AtlasImageType::IndexType currIndex;
							Image2DType::IndexType jpgIndex;
							for (int z = 0; z < zDim; z++)
							{
								for (int x = 0; x < xDim; x++)
								{
									currIndex[0] = x;
									currIndex[1] = sliceNum;
									currIndex[2] = z;    // z flipping 

									jpgIndex[0] = x;
									jpgIndex[1] = z;
									Image2DType::PixelType currPixel = jpgImg->GetPixel(jpgIndex);

									SAItkImage->SetPixel(currIndex, AtlasImageType::PixelType(currPixel));
								}
							}
						}
						catch (itk::ExceptionObject & err)
						{
							std::cerr << " Exception Object Caught in reading nrrd image !" << err << std::endl;
							return;
						}

						sliceNum++;

					}
					++it;
				} // end reading files 

				this->writeAtlasImage(SAItkImage, absName + ".nrrd");

			} // iterate all folders 

		} // if - folder is available

	} // folder is accepted 

	std::cout << " ----------- Writing SA.nrrd End--------------------- " << std::endl;
}


void SegmentationView::LoadAndMergeSCSlices()
{
	std::cout << " [ " << __FUNCTION__ << " ] Merge SC Slices .... " << std::endl;

	QFileDialog fd;
	fd.setDirectory(m_recentDir.c_str());
	fd.setFileMode(QFileDialog::DirectoryOnly);
	fd.setOption(QFileDialog::DontUseNativeDialog, true);

	QListView *list = fd.findChild< QListView* >("listView");
	if (list)  list->setSelectionMode(QAbstractItemView::MultiSelection);
	QTreeView *tree = fd.findChild< QTreeView* >();
	if (tree)  tree->setSelectionMode(QAbstractItemView::MultiSelection);

	if (fd.exec() == QDialog::Accepted)
	{
		m_recentDir = fd.directory().absolutePath().toAscii();
		QStringList myfiles = fd.selectedFiles();

		QString selected;

		if (!myfiles.isEmpty() && myfiles.size() > 1)
		{
			for (unsigned int dirCount = 1; dirCount <= myfiles.size() - 1; dirCount++)
			{
				selected = myfiles[dirCount];  // the first folder is empty 

				AtlasImageType::Pointer SCItkImage = AtlasImageType::New();
				AtlasImageType::RegionType region;
				AtlasImageType::SizeType size = { 256, 256, 256 };
				AtlasImageType::IndexType index = { 0, 0, 0 };
				region.SetIndex(index);
				region.SetSize(size);
				SCItkImage->SetRegions(region);
				SCItkImage->Allocate();
				SCItkImage->FillBuffer(itk::NumericTraits< AtlasImageType::PixelType >::Zero);

				QDir* imageDir = new QDir(selected);
				QStringList::Iterator it;
				QStringList files = imageDir->entryList();

				it = files.begin();
				int sliceNum = 0;

				std::stringstream ss; ss << dirCount;
				std::string absName = "SC_" + ss.str();

				while (it != files.end())
				{
					if (QFileInfo(*imageDir, *it).isFile() && ((*it).endsWith)(".jpg"))
					{
						mitk::DataNodeFactory::Pointer nodeReader = mitk::DataNodeFactory::New();

						QFileInfo& fileInfoIMG_ = QFileInfo(*imageDir, *it);
						std::string imgFileName = fileInfoIMG_.absoluteFilePath().toStdString();

						Image2DType::Pointer jpgImg = Image2DType::New();

						try
						{
							nodeReader->SetFileName(imgFileName);
							nodeReader->Update();

							mitk::DataNode::Pointer imageNode = nodeReader->GetOutput();
							mitk::Image::Pointer mitkImage = dynamic_cast<mitk::Image*>(imageNode->GetData());
							mitk::CastToItkImage(mitkImage, jpgImg);

							// copy the intensity to itk image 
							int xDim = jpgImg->GetLargestPossibleRegion().GetSize()[0];
							int yDim = jpgImg->GetLargestPossibleRegion().GetSize()[1];

							AtlasImageType::IndexType currIndex;
							Image2DType::IndexType jpgIndex;
							for (int y = 0; y < yDim; y++)
							{
								for (int x = 0; x < xDim; x++)
								{
									currIndex[0] = x;
									currIndex[1] = y;
									currIndex[2] = sliceNum;

									jpgIndex[0] = x;
									jpgIndex[1] = y;
									Image2DType::PixelType currPixel = jpgImg->GetPixel(jpgIndex);

									SCItkImage->SetPixel(currIndex, AtlasImageType::PixelType(currPixel));
								}
							}
						}
						catch (itk::ExceptionObject & err)
						{
							std::cerr << " Exception Object Caught in reading nrrd image !" << err << std::endl;
							return;
						}

						sliceNum++;
					}
					++it;
				}

				//	this->addImageToNode(SCItkImage, "SC" + AbsPath, true);

				this->writeAtlasImage(SCItkImage, absName + ".nrrd");

			} // iterate all folders 

		} // if - folder is available

	} // folder is accepted 

	std::cout << " ----------- Writing SC.nrrd End--------------------- " << std::endl;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//AtlasImageType::Pointer SCItkImage = AtlasImageType::New();
	//AtlasImageType::RegionType region;
	//AtlasImageType::SizeType size = { 256, 256, 256 };
	//AtlasImageType::IndexType index = { 0, 0, 0 };
	//region.SetIndex(index);
	//region.SetSize(size);
	//SCItkImage->SetRegions(region);
	//SCItkImage->Allocate();
	//SCItkImage->FillBuffer(itk::NumericTraits< AtlasImageType::PixelType >::Zero);
	//QFileDialog fd;
	//fd.setDirectory(m_recentDir.c_str());
	//fd.setFileMode(QFileDialog::Directory);
	//QString selected;
	//if (fd.exec() == QDialog::Accepted)
	//{
	//	m_recentDir = fd.directory().absolutePath().toAscii();
	//	QStringList myfiles = fd.selectedFiles();
	//	if (!myfiles.isEmpty())
	//		selected = myfiles[0];
	//}
	//QDir* imageDir = new QDir(selected);
	//QStringList::Iterator it;
	//QStringList files = imageDir->entryList();
	//it = files.begin();
	//int sliceNum = 0;
	//while (it != files.end())
	//{
	//	if (QFileInfo(*imageDir, *it).isFile() && ((*it).endsWith)(".jpg"))
	//	{
	//		mitk::DataNodeFactory::Pointer nodeReader = mitk::DataNodeFactory::New();
	//		QFileInfo& fileInfoIMG_ = QFileInfo(*imageDir, *it);
	//		std::string imgFileName = fileInfoIMG_.absoluteFilePath().toStdString();
	//		Image2DType::Pointer jpgImg = Image2DType::New();
	//		try
	//		{
	//			nodeReader->SetFileName(imgFileName);
	//			nodeReader->Update();
	//			mitk::DataNode::Pointer imageNode = nodeReader->GetOutput();
	//			mitk::Image::Pointer mitkImage = dynamic_cast<mitk::Image*>(imageNode->GetData());
	//			mitk::CastToItkImage(mitkImage, jpgImg);
	//			// copy the intensity to itk image 
	//			int xDim = jpgImg->GetLargestPossibleRegion().GetSize()[0];
	//			int yDim = jpgImg->GetLargestPossibleRegion().GetSize()[1];
	//			AtlasImageType::IndexType currIndex;
	//			Image2DType::IndexType jpgIndex;
	//			for (int y = 0; y < yDim; y++)
	//			{
	//				for (int x = 0; x < xDim; x++)
	//				{
	//					currIndex[0] = x;
	//					currIndex[1] = y;
	//					currIndex[2] = sliceNum;
	//					jpgIndex[0] = x;
	//					jpgIndex[1] = y;
	//					Image2DType::PixelType currPixel = jpgImg->GetPixel(jpgIndex);
	//					SCItkImage->SetPixel(currIndex, AtlasImageType::PixelType(currPixel));
	//				}
	//			}
	//		}
	//		catch (itk::ExceptionObject & err)
	//		{
	//			std::cerr << " Exception Object Caught in reading nrrd image !" << err << std::endl;
	//			return;
	//		}
	//		sliceNum++;
	//	}
	//	++it;
	//}
	//delete imageDir;
	//this->addImageToNode(SCItkImage, "SC", true);
	//this->writeAtlasImage(SCItkImage, "SC.nrrd"); 

}


void SegmentationView::LoadAndMergeCASlices()
{
	std::cout << " [ " << __FUNCTION__ << " ] Merge CA Slices .... " << std::endl;

	/*Image2DType::IndexType startIndex = { 71, 29 };
	Image2DType::IndexType endIndex = { 162, 202 }; */

	QFileDialog fd;
	fd.setDirectory(m_recentDir.c_str());
	fd.setFileMode(QFileDialog::DirectoryOnly);
	fd.setOption(QFileDialog::DontUseNativeDialog, true);

	QListView *list = fd.findChild< QListView* >("listView");
	if (list)  list->setSelectionMode(QAbstractItemView::MultiSelection);
	QTreeView *tree = fd.findChild< QTreeView* >();
	if (tree)  tree->setSelectionMode(QAbstractItemView::MultiSelection);

	if (fd.exec() == QDialog::Accepted)
	{
		m_recentDir = fd.directory().absolutePath().toAscii();
		QStringList myfiles = fd.selectedFiles();

		QString selected;

		if (!myfiles.isEmpty() && myfiles.size() > 1)
		{
			for (unsigned int dirCount = 1; dirCount <= myfiles.size() - 1; dirCount++)
			{
				selected = myfiles[dirCount];  // the first folder is empty 

				AtlasImageType::Pointer CAItkImage = AtlasImageType::New();
				AtlasImageType::RegionType region;
				AtlasImageType::SizeType size = { 256, 256, 256 };
				AtlasImageType::IndexType index = { 0, 0, 0 };
				region.SetIndex(index);
				region.SetSize(size);
				CAItkImage->SetRegions(region);
				CAItkImage->Allocate();
				CAItkImage->FillBuffer(itk::NumericTraits< AtlasImageType::PixelType >::Zero);

				QDir* imageDir = new QDir(selected);
				QStringList::Iterator it;
				QStringList files = imageDir->entryList();

				it = files.begin();
				int sliceNum = 0;

				std::stringstream ss; ss << dirCount;
				std::string absName = "CA_" + ss.str();

				while (it != files.end())
				{
					if (QFileInfo(*imageDir, *it).isFile() && ((*it).endsWith)(".jpg"))
					{
						mitk::DataNodeFactory::Pointer nodeReader = mitk::DataNodeFactory::New();

						QFileInfo& fileInfoIMG_ = QFileInfo(*imageDir, *it);
						std::string imgFileName = fileInfoIMG_.absoluteFilePath().toStdString();

						Image2DType::Pointer jpgImg = Image2DType::New();

						try
						{
							nodeReader->SetFileName(imgFileName);
							nodeReader->Update();

							mitk::DataNode::Pointer imageNode = nodeReader->GetOutput();
							mitk::Image::Pointer mitkImage = dynamic_cast<mitk::Image*>(imageNode->GetData());
							mitk::CastToItkImage(mitkImage, jpgImg);

							// copy the intensity to itk image 
							int yDim = jpgImg->GetLargestPossibleRegion().GetSize()[0];
							int zDim = jpgImg->GetLargestPossibleRegion().GetSize()[1];

							AtlasImageType::IndexType currIndex;
							Image2DType::IndexType jpgIndex;
							for (int z = 0; z < zDim; z++)
							{
								for (int y = 0; y < yDim; y++)
								{
									currIndex[0] = sliceNum;
									currIndex[1] = y;
									currIndex[2] = z;  // change, from down to up

									jpgIndex[0] = y;
									jpgIndex[1] = z;
									Image2DType::PixelType currPixel = jpgImg->GetPixel(jpgIndex);

									CAItkImage->SetPixel(currIndex, AtlasImageType::PixelType(currPixel));
								}
							}
						}
						catch (itk::ExceptionObject & err)
						{
							std::cerr << " Exception Object Caught in reading nrrd image !" << err << std::endl;
							return;
						}

						sliceNum++;

					}
					++it;
				}
				//	this->addImageToNode(SCItkImage, "SC" + AbsPath, true);

				// this->writeChar3DImage(CAItkImage, absName + ".nrrd"); 
				this->writeAtlasImage(CAItkImage, absName + ".nrrd");

			} // iterate all folders 

		} // if - folder is available

	} // folder is accepted 

	std::cout << " ----------- Writing CA.nrrd End--------------------- " << std::endl;

	//AtlasImageType::Pointer CAItkImage = AtlasImageType::New();
	//AtlasImageType::RegionType region;
	//AtlasImageType::SizeType size = { 256, 256, 256 };
	//AtlasImageType::IndexType index = { 0, 0, 0 };
	//region.SetIndex(index);
	//region.SetSize(size);
	//CAItkImage->SetRegions(region);
	//CAItkImage->Allocate();
	//CAItkImage->FillBuffer(itk::NumericTraits< AtlasImageType::PixelType >::Zero);
	//QFileDialog fd;
	//fd.setDirectory(m_recentDir.c_str());
	//fd.setFileMode(QFileDialog::Directory);
	//QString selected;
	//if (fd.exec() == QDialog::Accepted)
	//{
	//	m_recentDir = fd.directory().absolutePath().toAscii();
	//	QStringList myfiles = fd.selectedFiles();
	//	if (!myfiles.isEmpty())
	//		selected = myfiles[0];
	//}
	//QDir* imageDir = new QDir(selected);
	//QStringList::Iterator it;
	//QStringList files = imageDir->entryList();
	//it = files.begin();
	//int sliceNum = 0;
	//while (it != files.end())
	//{
	//	if (QFileInfo(*imageDir, *it).isFile() && ((*it).endsWith)(".jpg"))
	//	{
	//		mitk::DataNodeFactory::Pointer nodeReader = mitk::DataNodeFactory::New();
	//		QFileInfo& fileInfoIMG_ = QFileInfo(*imageDir, *it);
	//		std::string imgFileName = fileInfoIMG_.absoluteFilePath().toStdString();
	//		Image2DType::Pointer jpgImg = Image2DType::New();
	//		try
	//		{
	//			nodeReader->SetFileName(imgFileName);
	//			nodeReader->Update();
	//			mitk::DataNode::Pointer imageNode = nodeReader->GetOutput();
	//			mitk::Image::Pointer mitkImage = dynamic_cast<mitk::Image*>(imageNode->GetData());
	//			mitk::CastToItkImage(mitkImage, jpgImg);
	//			// copy the intensity to itk image 
	//			int yDim = jpgImg->GetLargestPossibleRegion().GetSize()[0];
	//			int zDim = jpgImg->GetLargestPossibleRegion().GetSize()[1];
	//			AtlasImageType::IndexType currIndex;
	//			Image2DType::IndexType jpgIndex;
	//			for (int z = 0; z < zDim; z++)
	//			{
	//				for (int y = 0; y < yDim; y++)
	//				{
	//					currIndex[0] = sliceNum;
	//					currIndex[1] = y ;  
	//					currIndex[2] = z ;  // change, from down to up
	//					jpgIndex[0] = y;
	//					jpgIndex[1] = z;
	//					Image2DType::PixelType currPixel = jpgImg->GetPixel(jpgIndex);
	//					CAItkImage->SetPixel(currIndex, AtlasImageType::PixelType(currPixel));
	//				}
	//			}
	//		}
	//		catch (itk::ExceptionObject & err)
	//		{
	//			std::cerr << " Exception Object Caught in reading nrrd image !" << err << std::endl;
	//			return;
	//		}
	//		sliceNum++;
	//	}
	//	++it;
	//}
	//delete imageDir;
	//this->addImageToNode(CAItkImage, "CA", true);
	//this->writeAtlasImage(CAItkImage, "CA.nrrd"); 
	//std::cout << " ----------- Writing CA.nrrd End --------------------- " << std::endl;
}


void SegmentationView::MergePancreasMap()
{
	AtlasImageType::IndexType minBound = { 90, 71, 29 };
	AtlasImageType::IndexType maxBound = { 203, 163, 203 };

	if (!this->m_probMap)
	{
		std::cout << " Please Select One Image To Crop ...  " << std::endl;
		return; 
	}

	cout << " produce image = " << this->m_probMapName << endl; 

	AtlasImageType::Pointer myProbMap = AtlasImageType::New();
	AtlasImageType::SizeType size = { 256, 256, 256 };
	AtlasImageType::IndexType index = { 0, 0, 0 };
	AtlasImageType::RegionType inputRegion(index, size);
	myProbMap->SetRegions(inputRegion);
	myProbMap->SetOrigin(this->m_probMap->GetOrigin());
	myProbMap->SetSpacing(this->m_probMap->GetSpacing());
	myProbMap->Allocate();
	myProbMap->FillBuffer(itk::NumericTraits< AtlasImageType::PixelType >::Zero);

	QFileDialog fd;
	fd.setDirectory(m_recentDir.c_str());
	fd.setFileMode(QFileDialog::Directory);
	QString selected;

	if (fd.exec() == QDialog::Accepted)
	{
		m_recentDir = fd.directory().absolutePath().toAscii();
		QStringList myfiles = fd.selectedFiles();
		if (!myfiles.isEmpty())
			selected = myfiles[0];
	}

	QDir* imageDir = new QDir(selected);
	QStringList::Iterator it;
	QStringList files = imageDir->entryList();
	it = files.begin();

	int Cnt = 0;

	while (it != files.end())
	{
		if (QFileInfo(*imageDir, *it).isFile() && ((*it).endsWith)(".pic") || (*it).endsWith(".nrrd") || (*it).endsWith(".mhd"))
		{
			mitk::DataNodeFactory::Pointer nodeReader = mitk::DataNodeFactory::New();

			QFileInfo& fileInfoIMG_ = QFileInfo(*imageDir, *it);
			std::string imgFileName = fileInfoIMG_.absoluteFilePath().toStdString();
			std::string imgAbsName = QFileInfo(*imageDir, *it).baseName().toStdString();

			AtlasImageType::Pointer itkImage = AtlasImageType::New();
			AtlasImageType::Pointer newImg = AtlasImageType::New();
			newImg->SetRegions(inputRegion);
			newImg->SetSpacing(this->m_probMap->GetSpacing());
			newImg->SetOrigin(this->m_probMap->GetOrigin());
			newImg->Allocate();
			newImg->FillBuffer(itk::NumericTraits<AtlasImageType::PixelType>::Zero);

			Cnt++;

			try{

				nodeReader->SetFileName(imgFileName);
				nodeReader->Update();

				mitk::DataNode::Pointer imageNode = nodeReader->GetOutput();
				mitk::Image::Pointer mitkImage = dynamic_cast<mitk::Image*>(imageNode->GetData());
				mitk::CastToItkImage(mitkImage, itkImage);
			}
			catch (...){}

			itk::ImageRegionIterator< AtlasImageType > it(itkImage, itkImage->GetRequestedRegion());
			itk::ImageRegionIterator< AtlasImageType > prob(myProbMap, myProbMap->GetRequestedRegion());

			for (it.GoToBegin(), prob.GoToBegin(); !it.IsAtEnd() && !prob.IsAtEnd(); ++it, ++prob)
			{
				if (it.GetIndex()[0] >= minBound[0] && it.GetIndex()[0] <= maxBound[0] &&
					it.GetIndex()[1] >= minBound[1] && it.GetIndex()[1] <= maxBound[1] &&
					it.GetIndex()[2] >= minBound[2] && it.GetIndex()[2] <= maxBound[2] && it.Get() > 20)
				{

					prob.Set(prob.Get() + it.Get());

				}
				else
					prob.Set(0);
			}
		}

		++it;
	}

	delete imageDir;

	//	this->addImageToNode(myProbMap, "map", true);
	//	this->writeAtlasImage(myProbMap, this->m_fixedImgName + "_map.nrrd"); 

	//// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Pancreas 


	itk::ImageRegionIterator< AtlasImageType > initProb(myProbMap, myProbMap->GetRequestedRegion());
	for (initProb.GoToBegin(); !initProb.IsAtEnd(); ++initProb)
	{
		if (initProb.GetIndex()[0] >= minBound[0] && initProb.GetIndex()[0] <= maxBound[0] &&
			initProb.GetIndex()[1] >= minBound[1] && initProb.GetIndex()[1] <= maxBound[1] &&
			initProb.GetIndex()[2] >= minBound[2] && initProb.GetIndex()[2] <= maxBound[2] &&
			initProb.Get() >= m_Controls.closeRadius->value())
		{

		}
		else
			initProb.Set(0);
	}

	// get the connected components 
	typedef itk::ConnectedComponentImageFilter< AtlasImageType, AtlasImageType> ConnectedComponentImageFilterType;
	ConnectedComponentImageFilterType::Pointer connectedFilter = ConnectedComponentImageFilterType::New();
	connectedFilter->SetInput(myProbMap);
	connectedFilter->Update();

	typedef itk::RelabelComponentImageFilter< AtlasImageType, AtlasImageType > RelabelFilterType;
	RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
	relabelFilter->SetInput(connectedFilter->GetOutput());
	relabelFilter->Update();

	typedef itk::LabelShapeKeepNObjectsImageFilter< AtlasImageType > LabelShapeKeepNObjectsImageFilterType;
	LabelShapeKeepNObjectsImageFilterType::Pointer labelFilter = LabelShapeKeepNObjectsImageFilterType::New();
	labelFilter->SetInput(relabelFilter->GetOutput());
	labelFilter->SetBackgroundValue(0);
	labelFilter->SetNumberOfObjects(1);   // right 
	labelFilter->SetAttribute(LabelShapeKeepNObjectsImageFilterType::LabelObjectType::PHYSICAL_SIZE);
	labelFilter->Update();

	this->writeAtlasImage(labelFilter->GetOutput(), this->m_probMapName + "_map.nrrd");
	//
	//	typedef itk::MultiplyImageFilter< AtlasImageType, AtlasImageType > MultiplyFilterType; 
	//	MultiplyFilterType::Pointer multiply = MultiplyFilterType::New(); 
	//	multiply->SetInput1(myProbMap); 
	//	multiply->SetInput2(labelFilter->GetOutput()); 
	//	multiply->Update(); 
	//
	///////////////// another /////////////////////////////////////
	//
	// convert to binary image 
	/*typedef itk::BinaryThresholdImageFilter< AtlasImageType, AtlasImageType > BinaryThresholdImageFilterType;
	BinaryThresholdImageFilterType::Pointer binaryFilter = BinaryThresholdImageFilterType::New();

	binaryFilter->SetInput(myProbMap);
	binaryFilter->SetLowerThreshold(0);
	binaryFilter->SetUpperThreshold(100);
	binaryFilter->SetInsideValue(0);
	binaryFilter->SetOutsideValue(1);
	binaryFilter->Update();

	this->addImageToNode(binaryFilter->GetOutput(), "bi", true); */
	//
	//	typedef itk::RescaleIntensityImageFilter< AtlasImageType, AtlasImageType > RescaleFilterType;
	//	RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	//	rescaleFilter->SetOutputMinimum(0);
	//	rescaleFilter->SetOutputMaximum(1);
	//	rescaleFilter->SetInput(binaryFilter->GetOutput());
	//	rescaleFilter->Update();
	//	
	//	this->writeAtlasImage(rescaleFilter->GetOutput(), this->m_fixedImgName + "_prob.nrrd");


	//// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	//	AtlasImageType::Pointer prob_leftMap = AtlasImageType::New();
	//	prob_leftMap->SetRegions(region);
	//	prob_leftMap->SetOrigin(inputImg->GetOrigin());
	//	prob_leftMap->SetSpacing(inputImg->GetSpacing());
	//	prob_leftMap->Allocate();
	//	prob_leftMap->FillBuffer(itk::NumericTraits< AtlasImageType::PixelType >::Zero);
	//
	//	// create a copy for right_kidney 
	//	for (int x = 0; x < 128; x++)
	//	{
	//		for (int y = 0; y < prob_leftMap->GetLargestPossibleRegion().GetSize()[1]; y++)
	//		{
	//			for (int z = 0; z < prob_leftMap->GetLargestPossibleRegion().GetSize()[2]; z++)
	//			{
	//				AtlasImageType::IndexType index = { x, y, z }; 
	//				if (binaryFilter->GetOutput()->GetPixel(index) > 0)
	//					prob_leftMap->SetPixel(index, 1); 
	//			}
	//		}
	//	}
	//
	//	// get the connected components 
	//	typedef itk::ConnectedComponentImageFilter< AtlasImageType, AtlasImageType> ConnectedComponentImageFilterType;
	//	ConnectedComponentImageFilterType::Pointer connectedFilterLeft = ConnectedComponentImageFilterType::New();
	//	connectedFilterLeft->SetInput(prob_leftMap);
	//	connectedFilterLeft->Update();
	//
	//	typedef itk::RelabelComponentImageFilter< AtlasImageType, AtlasImageType > RelabelFilterType;
	//	RelabelFilterType::Pointer relabelFilterLeft = RelabelFilterType::New();
	//	relabelFilterLeft->SetInput(connectedFilterLeft->GetOutput());
	//	relabelFilterLeft->Update();
	//
	//	typedef itk::LabelShapeKeepNObjectsImageFilter< AtlasImageType > LabelShapeKeepNObjectsImageFilterType;
	//	LabelShapeKeepNObjectsImageFilterType::Pointer labelFilter_left = LabelShapeKeepNObjectsImageFilterType::New();
	//	labelFilter_left->SetInput(relabelFilterLeft->GetOutput());
	//	labelFilter_left->SetBackgroundValue(0);
	//	labelFilter_left->SetNumberOfObjects(1);   // right 
	//	labelFilter_left->SetAttribute(LabelShapeKeepNObjectsImageFilterType::LabelObjectType::PHYSICAL_SIZE);
	//	labelFilter_left->Update();
	//
	//	typedef itk::RescaleIntensityImageFilter< AtlasImageType, AtlasImageType > RescaleFilterType;
	//	RescaleFilterType::Pointer rescaleFilter_left = RescaleFilterType::New();
	//	rescaleFilter_left->SetOutputMinimum(0);
	//	rescaleFilter_left->SetOutputMaximum(1);
	//	rescaleFilter_left->SetInput(labelFilter_left->GetOutput());
	//	rescaleFilter_left->Update();
	//
	//	this->writeAtlasImage(rescaleFilter_left->GetOutput(), this->m_fixedImgName + "_right.nrrd");
	//
	//
	////// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	//
	//	AtlasImageType::Pointer prob_rightMap = AtlasImageType::New();
	//	prob_rightMap->SetRegions(region);
	//	prob_rightMap->SetOrigin(inputImg->GetOrigin());
	//	prob_rightMap->SetSpacing(inputImg->GetSpacing());
	//	prob_rightMap->Allocate();
	//	prob_rightMap->FillBuffer(itk::NumericTraits< AtlasImageType::PixelType >::Zero);
	//
	//
	//	// create a copy for left_kidney 
	//	for (int x = 128; x < prob_rightMap->GetLargestPossibleRegion().GetSize()[0] ; x++)
	//	{
	//		for (int y = 0; y < prob_rightMap->GetLargestPossibleRegion().GetSize()[1]; y++)
	//		{
	//			for (int z = 0; z < prob_rightMap->GetLargestPossibleRegion().GetSize()[2] ; z++)
	//			{
	//				AtlasImageType::IndexType index = { x, y, z };
	//				if (binaryFilter->GetOutput()->GetPixel(index) > 0)
	//					prob_rightMap->SetPixel(index, 1);	
	//			}
	//		}
	//	}
	//
	//	// get the connected components 
	//	ConnectedComponentImageFilterType::Pointer connectedFilterRight = ConnectedComponentImageFilterType::New();
	//	connectedFilterRight->SetInput(prob_rightMap);
	//	connectedFilterRight->Update();
	//
	//	RelabelFilterType::Pointer relabelFilterRight = RelabelFilterType::New();
	//	relabelFilterRight->SetInput(connectedFilterRight->GetOutput());
	//	relabelFilterRight->Update();
	//
	//	LabelShapeKeepNObjectsImageFilterType::Pointer labelFilterRight = LabelShapeKeepNObjectsImageFilterType::New();
	//	labelFilterRight->SetInput(relabelFilterRight->GetOutput());
	//	labelFilterRight->SetBackgroundValue(0);
	//	labelFilterRight->SetNumberOfObjects(1);   // right 
	//	labelFilterRight->SetAttribute(LabelShapeKeepNObjectsImageFilterType::LabelObjectType::PHYSICAL_SIZE);
	//	labelFilterRight->Update();
	//
	//	RescaleFilterType::Pointer rescaleFilterRight = RescaleFilterType::New();
	//	rescaleFilterRight->SetOutputMinimum(0);
	//	rescaleFilterRight->SetOutputMaximum(1);
	//	rescaleFilterRight->SetInput(labelFilterRight->GetOutput());
	//	rescaleFilterRight->Update();
	//
	//	this->writeAtlasImage(rescaleFilterRight->GetOutput(), this->m_fixedImgName + "_left.nrrd");

	std::cout << "--------------- Probability Pancreas Map Creation Done --------------------------" << std::endl;
}


void SegmentationView::MergeKidneyMap()
{
	if (!this->m_fixedImageForReg)
	{
		std::cout << " Please select one image to crop ... " << this->m_fixedImgName << std::endl;
		return;
	}

	typedef itk::ConnectedComponentImageFilter< AtlasImageType, AtlasImageType> ConnectedComponentImageFilterType;
	typedef itk::RelabelComponentImageFilter< AtlasImageType, AtlasImageType > RelabelFilterType;
	typedef itk::LabelShapeKeepNObjectsImageFilter< AtlasImageType > LabelShapeKeepNObjectsImageFilterType;
	typedef itk::RescaleIntensityImageFilter< AtlasImageType, AtlasImageType > RescaleFilterType;

	AtlasImageType::Pointer inputImg = AtlasImageType::New();
	mitk::CastToItkImage(this->m_fixedImageForReg, inputImg);

	AtlasImageType::Pointer myProbMap = AtlasImageType::New();
	AtlasImageType::SizeType size = { 256, 256, 256 };
	AtlasImageType::IndexType index = { 0, 0, 0 };
	AtlasImageType::RegionType region(index, size);
	myProbMap->SetRegions(region);
	myProbMap->SetOrigin(inputImg->GetOrigin());
	myProbMap->SetSpacing(inputImg->GetSpacing());
	myProbMap->Allocate();
	myProbMap->FillBuffer(itk::NumericTraits< AtlasImageType::PixelType >::Zero);

	QFileDialog fd;
	fd.setDirectory(m_recentDir.c_str());
	fd.setFileMode(QFileDialog::Directory);
	QString selected;

	if (fd.exec() == QDialog::Accepted)
	{
		m_recentDir = fd.directory().absolutePath().toAscii();
		QStringList myfiles = fd.selectedFiles();
		if (!myfiles.isEmpty())
			selected = myfiles[0];
	}

	QDir* imageDir = new QDir(selected);
	QStringList::Iterator it;
	QStringList files = imageDir->entryList();
	it = files.begin();

	int Cnt = 0;

	while (it != files.end())
	{
		if (QFileInfo(*imageDir, *it).isFile() && ((*it).endsWith)(".pic") || (*it).endsWith(".nrrd") || (*it).endsWith(".mhd"))
		{
			mitk::DataNodeFactory::Pointer nodeReader = mitk::DataNodeFactory::New();

			QFileInfo& fileInfoIMG_ = QFileInfo(*imageDir, *it);
			std::string imgFileName = fileInfoIMG_.absoluteFilePath().toStdString();
			std::string imgAbsName = QFileInfo(*imageDir, *it).baseName().toStdString();

			AtlasImageType::Pointer itkImage = AtlasImageType::New();

			Cnt++;

			try{

				nodeReader->SetFileName(imgFileName);
				nodeReader->Update();

				mitk::DataNode::Pointer imageNode = nodeReader->GetOutput();
				mitk::Image::Pointer mitkImage = dynamic_cast<mitk::Image*>(imageNode->GetData());
				mitk::CastToItkImage(mitkImage, itkImage);
			}
			catch (...){}

			itk::ImageRegionIterator< AtlasImageType > it(itkImage, itkImage->GetRequestedRegion());
			itk::ImageRegionIterator< AtlasImageType > prob(myProbMap, myProbMap->GetRequestedRegion());
			for (it.GoToBegin(), prob.GoToBegin(); !it.IsAtEnd() && !prob.IsAtEnd(); ++it, ++prob)
			{
				prob.Set(prob.Get() + it.Get());
			}
		}

		++it;
	}

	delete imageDir;

	this->writeAtlasImage(myProbMap, this->m_fixedImgName + "_map.nrrd");

	//// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Kidney >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	AtlasImageType::IndexType minBound = { 0, 0, 0 };
	AtlasImageType::IndexType maxBound = { 254, 254, 254 };

	itk::ImageRegionIterator< AtlasImageType > initProb(myProbMap, myProbMap->GetRequestedRegion());
	for (initProb.GoToBegin(); !initProb.IsAtEnd(); ++initProb)
	{
		if (initProb.GetIndex()[0] >= minBound[0] && initProb.GetIndex()[0] <= maxBound[0] &&
			initProb.GetIndex()[1] >= minBound[1] && initProb.GetIndex()[1] <= maxBound[1] &&
			initProb.GetIndex()[2] >= minBound[2] && initProb.GetIndex()[2] <= maxBound[2] && 
			initProb.Get() > 100){
		}
		else
			initProb.Set(0);
	}
	
	//// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	AtlasImageType::Pointer prob_leftMap = AtlasImageType::New();
	prob_leftMap->SetRegions(region);
	prob_leftMap->SetOrigin(inputImg->GetOrigin());
	prob_leftMap->SetSpacing(inputImg->GetSpacing());
	prob_leftMap->Allocate();
	prob_leftMap->FillBuffer(itk::NumericTraits< AtlasImageType::PixelType >::Zero);

	// create a copy for right_kidney 
	for (int x = minBound[0]; x < 128; x++)
	{
		for (int y = minBound[1]; y < maxBound[1]; y++)
		{
			for (int z = minBound[2]; z < maxBound[2]; z++)
			{
				AtlasImageType::IndexType index = { x, y, z };
				prob_leftMap->SetPixel(index, myProbMap->GetPixel(index));
			}
		}
	}

	
	ConnectedComponentImageFilterType::Pointer connectedFilterLeft = ConnectedComponentImageFilterType::New();
	connectedFilterLeft->SetInput(prob_leftMap);
	connectedFilterLeft->Update();

	
	RelabelFilterType::Pointer relabelFilterLeft = RelabelFilterType::New();
	relabelFilterLeft->SetInput(connectedFilterLeft->GetOutput());
	relabelFilterLeft->Update();

	
	LabelShapeKeepNObjectsImageFilterType::Pointer labelFilter_left = LabelShapeKeepNObjectsImageFilterType::New();
	labelFilter_left->SetInput(relabelFilterLeft->GetOutput());
	labelFilter_left->SetBackgroundValue(0);
	labelFilter_left->SetNumberOfObjects(1);   // right 
	labelFilter_left->SetAttribute(LabelShapeKeepNObjectsImageFilterType::LabelObjectType::PHYSICAL_SIZE);
	labelFilter_left->Update();

	itk::ImageRegionIterator< AtlasImageType > leftProb(prob_leftMap, prob_leftMap->GetRequestedRegion()); 
	itk::ImageRegionIterator< AtlasImageType > leftMask(labelFilter_left->GetOutput(), labelFilter_left->GetOutput()->GetRequestedRegion()); 
	for (leftProb.GoToBegin(), leftMask.GoToBegin(); !leftProb.IsAtEnd() && !leftMask.IsAtEnd(); ++leftProb, ++leftMask)
	{
		if (leftMask.Get() == 0) leftProb.Set(0); 
	}

	this->writeAtlasImage(prob_leftMap, this->m_fixedImgName + "_right.nrrd");


	/*RescaleFilterType::Pointer rescaleFilter_left = RescaleFilterType::New();
	rescaleFilter_left->SetOutputMinimum(0);
	rescaleFilter_left->SetOutputMaximum(1);
	rescaleFilter_left->SetInput(labelFilter_left->GetOutput());
	rescaleFilter_left->Update();*/

	
	//// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	AtlasImageType::Pointer prob_rightMap = AtlasImageType::New();
	prob_rightMap->SetRegions(region);
	prob_rightMap->SetOrigin(inputImg->GetOrigin());
	prob_rightMap->SetSpacing(inputImg->GetSpacing());
	prob_rightMap->Allocate();
	prob_rightMap->FillBuffer(itk::NumericTraits< AtlasImageType::PixelType >::Zero);


	// create a copy for left_kidney 
	for (int x = 128; x < maxBound[0]; x++)
	{
		for (int y = minBound[1]; y < maxBound[1]; y++)
		{
			for (int z = minBound[2]; z < maxBound[2]; z++)
			{
				AtlasImageType::IndexType index = { x, y, z };
				prob_rightMap->SetPixel(index, myProbMap->GetPixel(index));
			}
		}
	}

	// get the connected components 
	ConnectedComponentImageFilterType::Pointer connectedFilterRight = ConnectedComponentImageFilterType::New();
	connectedFilterRight->SetInput(prob_rightMap);
	connectedFilterRight->Update();

	RelabelFilterType::Pointer relabelFilterRight = RelabelFilterType::New();
	relabelFilterRight->SetInput(connectedFilterRight->GetOutput());
	relabelFilterRight->Update();

	LabelShapeKeepNObjectsImageFilterType::Pointer labelFilter_right = LabelShapeKeepNObjectsImageFilterType::New();
	labelFilter_right->SetInput(relabelFilterRight->GetOutput());
	labelFilter_right->SetBackgroundValue(0);
	labelFilter_right->SetNumberOfObjects(1);   // right 
	labelFilter_right->SetAttribute(LabelShapeKeepNObjectsImageFilterType::LabelObjectType::PHYSICAL_SIZE);
	labelFilter_right->Update();

	itk::ImageRegionIterator< AtlasImageType > rightProb(prob_rightMap, prob_rightMap->GetRequestedRegion());
	itk::ImageRegionIterator< AtlasImageType > rightMask(labelFilter_right->GetOutput(), labelFilter_right->GetOutput()->GetRequestedRegion());
	for (rightProb.GoToBegin(), rightMask.GoToBegin(); !rightProb.IsAtEnd() && !rightMask.IsAtEnd(); ++rightProb, ++rightMask)
	{
		if (rightMask.Get() == 0) rightProb.Set(0);
	}

	this->writeAtlasImage(prob_rightMap, this->m_fixedImgName + "_left.nrrd");


	/*RescaleFilterType::Pointer rescaleFilterRight = RescaleFilterType::New();
	rescaleFilterRight->SetOutputMinimum(0);
	rescaleFilterRight->SetOutputMaximum(1);
	rescaleFilterRight->SetInput(labelFilterRight->GetOutput());
	rescaleFilterRight->Update();*/

	std::cout << "--------------- Probability Map Creation Done --------------------------" << std::endl;
}


/////////////////////////// --- Important ToolBox ---- //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void SegmentationView::FixImageSelectedForReg(const mitk::DataNode* _node)
{
	this->m_fixedImageForReg = mitk::Image::Pointer();

	if (_node){
		std::string name;
		if (_node->GetName(name)){
			this->m_fixedImgName = _node->GetName();

			mitk::BaseData* basedata = _node->GetData();
			if (basedata){
				this->m_fixedImageForReg = dynamic_cast<mitk::Image*>(basedata);
			}
		}
	}
}


void SegmentationView::LoadAndResampleImage()
{
	cout << __FUNCTION__ << " Please Load and Resampling Images ... " << endl; 

	// import masks 
	QFileDialog fd;
	fd.setDirectory(m_recentDir.c_str());
	fd.setFileMode(QFileDialog::Directory);
	QString selected;

	if (fd.exec() == QDialog::Accepted)
	{
		m_recentDir = fd.directory().absolutePath().toAscii();
		QStringList myfiles = fd.selectedFiles();
		if (!myfiles.isEmpty())
			selected = myfiles[0];
	}

	QDir* imageDir = new QDir(selected);
	QStringList::Iterator it;
	QStringList files = imageDir->entryList();

	it = files.begin();

	while (it != files.end())
	{
		if (QFileInfo(*imageDir, *it).isFile() && ((*it).endsWith)(".nii.gz") || (*it).endsWith(".nrrd") || (*it).endsWith(".mhd"))
		{
			mitk::DataNodeFactory::Pointer nodeReader = mitk::DataNodeFactory::New();

			QFileInfo& fileInfoIMG_ = QFileInfo(*imageDir, *it);
			std::string imgFileName = fileInfoIMG_.absoluteFilePath().toStdString();
			std::string imgAbsName = QFileInfo(*imageDir, *it).baseName().toStdString();

			AtlasImageType::Pointer inputImg = AtlasImageType::New();

			try{
				nodeReader->SetFileName(imgFileName);
				nodeReader->Update();

				mitk::DataNode::Pointer imageNode = nodeReader->GetOutput();
				mitk::Image::Pointer mitkImage = dynamic_cast<mitk::Image*>(imageNode->GetData());
				mitk::CastToItkImage(mitkImage, inputImg);
			}
			catch (...){}

			// resampling image 
			AtlasImageType::SizeType inputSize;
			inputSize[0] = inputImg->GetLargestPossibleRegion().GetSize()[0];
			inputSize[1] = inputImg->GetLargestPossibleRegion().GetSize()[1];
			inputSize[2] = inputImg->GetLargestPossibleRegion().GetSize()[2];

			std::cout << __FUNCTION__ << inputSize[0] << " , " << inputSize[1] << " , " << inputSize[2] << std::endl;

			AtlasImageType::SpacingType inputSpacing;
			inputSpacing[0] = inputImg->GetSpacing()[0];
			inputSpacing[1] = inputImg->GetSpacing()[1];
			inputSpacing[2] = inputImg->GetSpacing()[2];

			std::cout << __FUNCTION__ << inputSpacing[0] << " , " << inputSpacing[1] << " , " << inputSpacing[2] << std::endl;

			AtlasImageType::SizeType outputSize;
			outputSize.Fill(0);
			outputSize[0] = int(m_Controls.resample_x->value());
			outputSize[1] = int(m_Controls.resample_y->value());
			outputSize[2] = int(m_Controls.resample_z->value());

			std::cout << __FUNCTION__ << outputSize[0] << " , " << outputSize[1] << " , " << outputSize[2] << std::endl;

			AtlasImageType::SpacingType outputSpacing;
			outputSpacing.Fill(0.0);
			outputSpacing[0] = inputSpacing[0] * (static_cast<double>(inputSize[0]) / static_cast<double>(outputSize[0]));
			outputSpacing[1] = inputSpacing[1] * (static_cast<double>(inputSize[1]) / static_cast<double>(outputSize[1]));
			outputSpacing[2] = inputSpacing[2] * (static_cast<double>(inputSize[2]) / static_cast<double>(outputSize[2]));

			std::cout << __FUNCTION__ << outputSpacing[0] << " , " << outputSpacing[1] << " , " << outputSpacing[2] << std::endl;

			typedef itk::ResampleImageFilter< AtlasImageType, AtlasImageType > ResampleFilterType;
			ResampleFilterType::Pointer resampler = ResampleFilterType::New();

			typedef itk::TranslationTransform< double, 3 > TransformType;
			TransformType::Pointer transform = TransformType::New();
			resampler->SetInput(inputImg);
			resampler->SetSize(outputSize);
			resampler->SetOutputSpacing(outputSpacing);
			resampler->SetTransform(transform);
			resampler->SetOutputOrigin(inputImg->GetOrigin());
			resampler->SetOutputDirection(inputImg->GetDirection());
			resampler->UpdateLargestPossibleRegion();

			AtlasImageType::Pointer output = resampler->GetOutput();
			std::cout << " Output Size = " << output->GetLargestPossibleRegion().GetSize() << std::endl;
			std::cout << " [ " << __FUNCTION__ << " ] resampled size[2] = " << resampler->GetOutput()->GetLargestPossibleRegion().GetSize()[2] << std::endl;


			// write out resampled image 
			this->writeAtlasImage(output, imgAbsName + ".nrrd");

		}

		++it;
	}

	delete imageDir;

	/////////////////////////// Resampling ///////////////////////////////////////////////////////

	
	


	///////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////// Remove Boundary /////////////////////////////////////////////////////////
	//// remove the boundary pixel = 1 
	//bool found = false; 
	//for (int zDim = 250; zDim < itkImage->GetLargestPossibleRegion().GetSize()[2]; zDim++)
	//{
	//	for (int xDim = 0; xDim < itkImage->GetLargestPossibleRegion().GetSize()[0]; xDim++)
	//	{
	//		for (int yDim = 0; yDim < itkImage->GetLargestPossibleRegion().GetSize()[1]; yDim++)
	//		{
	//			Image3DType::IndexType index = { xDim, yDim, zDim }; 
	//			if (itkImage->GetPixel(index) > 0)
	//			{
	//				found = true; 
	//				itkImage->SetPixel(index, 0); 
	//			}
	//		}
	//	}
	//}
	//if (found)
	//{
	//	std::cout << " !!!!! found in higher range of image : " << imgAbsName << std::endl;
	//}
	//itkImage->Update();
	//typedef itk::ImageFileWriter< Image3DType > WriterType; 
	//WriterType::Pointer writer = WriterType::New(); 
	//writer->SetInput(itkImage); 
	//writer->SetFileName(imgAbsName + ".nrrd"); 
	//writer->Update(); 


}


void SegmentationView::ComputeGMM()
{
	if (!this->m_segImage)
	{
		MITK_WARN << " Please load an image ... " << endl; 
		return; 
	}

	if (!this->m_probMap)
	{
		MITK_WARN << " Please load a mask ... " << endl; 
		return; 
	}

	DuplicatorType::Pointer dup = DuplicatorType::New(); 
	dup->SetInputImage(this->m_segImage); 
	dup->Update(); 
	AtlasImageType::Pointer roiImg = dup->GetOutput(); 
	DuplicatorType::Pointer dup1 = DuplicatorType::New();
	dup1->SetInputImage(this->m_segImage);
	dup1->Update();
	AtlasImageType::Pointer noiImg = dup1->GetOutput();


	itk::ImageRegionIterator< AtlasImageType > roi(roiImg, roiImg->GetRequestedRegion());
	itk::ImageRegionIterator< AtlasImageType > noi(noiImg, noiImg->GetRequestedRegion()); 
	itk::ImageRegionIterator< AtlasImageType > mask(this->m_probMap, this->m_probMap->GetRequestedRegion());
	
	for (roi.GoToBegin(), noi.GoToBegin(), mask.GoToBegin(); !roi.IsAtEnd() && !noi.IsAtEnd() && !mask.IsAtEnd(); ++roi, ++noi, ++mask)
	{
		if (mask.Value() > 0.0) noi.Set(-2000);
		
		else roi.Set(-2000);
	}

	cout << " ================= ROI Statistics =========================== " << endl; 
	statisticFilterType::Pointer statFilter = statisticFilterType::New();
	statFilter->SetInput(roiImg);
	statFilter->Update();

	cout << " ================= NOI Statistics =========================== " << endl;
	statFilter->SetInput(noiImg); 
	statFilter->Update(); 

	cout << " ============================================================ " << endl;


}


void SegmentationView::GenerateCorrespondingMeshFromBinaryImage()
{
	int numVertices = m_Controls.numVertices->value();
	unsigned int radiusOpening = 0;
	unsigned int radiusClosing = 0; //5;

	typedef itk::BinaryThresholdImageFilter< Image3DType, Image3DType> BinaryFilterType; 

	//== Load Binary Image Folder
	QFileDialog fd;
	fd.setDirectory(m_recentDir.c_str());
	fd.setFileMode(QFileDialog::Directory);
	QString selected;

	if (fd.exec() == QDialog::Accepted)
	{
		m_recentDir = fd.directory().absolutePath().toAscii();
		QStringList myfiles = fd.selectedFiles();
		if (!myfiles.isEmpty())
			selected = myfiles[0];
	}

	QDir* imageDir = new QDir(selected);
	QStringList::Iterator it;
	QStringList files = imageDir->entryList();

	it = files.begin();
	int Cnt = 0;

	while (it != files.end())
	{
		if (QFileInfo(*imageDir, *it).isFile() && ((*it).endsWith)(".pic") || (*it).endsWith(".nrrd") || (*it).endsWith(".mhd") || (*it).endsWith(".nii.gz"))
		{
			Cnt++;

			mitk::DataNodeFactory::Pointer nodeReader = mitk::DataNodeFactory::New();

			QFileInfo& fileInfoIMG_ = QFileInfo(*imageDir, *it);
			std::string imgFileName = fileInfoIMG_.absoluteFilePath().toStdString();
			std::string imgAbsName = QFileInfo(*imageDir, *it).baseName().toStdString();

			Image3DType::Pointer maskImg = Image3DType::New();

			try{
				nodeReader->SetFileName(imgFileName);
				nodeReader->Update();

				mitk::DataNode::Pointer imageNode = nodeReader->GetOutput();
				mitk::Image::Pointer mitkImage = dynamic_cast<mitk::Image*>(imageNode->GetData());
				mitk::CastToItkImage(mitkImage, maskImg);
			}
			catch (itk::ExceptionObject & err)
			{
				std::cerr << "Exception Object Caught in reading nrrd image !" << err << std::endl;
				return;
			}

			// cast to binary char image 
			BinaryFilterType::Pointer binary = BinaryFilterType::New();
			binary->SetOutsideValue(0); 
			binary->SetInsideValue(1); 
			binary->SetLowerThreshold(1); 
			binary->SetUpperThreshold(1000); 
			binary->SetInput(maskImg); 
			binary->Update(); 


			//>>> perform opening to the mask image 
			typedef itk::BinaryBallStructuringElement<Image3DType::PixelType, Image3DType::ImageDimension> StructuringElementType;
			StructuringElementType openingStructuringElement;
			openingStructuringElement.SetRadius(radiusOpening);
			openingStructuringElement.CreateStructuringElement();

			typedef itk::BinaryMorphologicalOpeningImageFilter< Image3DType, Image3DType, StructuringElementType> BinaryMorphologicalOpeningImageFilterType;
			BinaryMorphologicalOpeningImageFilterType::Pointer openingFilter = BinaryMorphologicalOpeningImageFilterType::New();
			openingFilter->SetInput(binary->GetOutput());
			openingFilter->SetForegroundValue(1.0);
			openingFilter->SetKernel(openingStructuringElement);
			openingFilter->Update();


			//>>> perform closing to the mask image 
			StructuringElementType closingStructuringElement;
			closingStructuringElement.SetRadius(radiusClosing);
			closingStructuringElement.CreateStructuringElement();

			typedef itk::BinaryMorphologicalClosingImageFilter< Image3DType, Image3DType, StructuringElementType > BinaryMorphologicalClosingImageFilterType;
			BinaryMorphologicalClosingImageFilterType::Pointer closingFilter = BinaryMorphologicalClosingImageFilterType::New();
			closingFilter->SetInput(openingFilter->GetOutput());
			closingFilter->SetForegroundValue(1);
			closingFilter->SetKernel(closingStructuringElement);
			closingFilter->Update();

			maskImg = closingFilter->GetOutput();

			this->createMeshFromBinaryImage(maskImg, imgAbsName, numVertices);

		}
		++it;
	}

	delete imageDir;

	std::cout << " Processing meshes done! " << std::endl;
}


void SegmentationView::FilterRegionOfInterest()
{
	std::cout << "  start filtering image - " << this->m_fixedImgName << std::endl;

	if (!this->m_fixedImageForReg)
	{
		std::cout << " Please select one image to crop ... " << std::endl;
		return;
	}

	AtlasImageType::Pointer filteredImg = AtlasImageType::New();
	mitk::CastToItkImage(this->m_fixedImageForReg, filteredImg);

	// Set Default Pixel 
	AtlasImageType::PixelType LK = m_Controls.HU_leftKidney->value();
	AtlasImageType::PixelType RK = m_Controls.HU_rightKidney->value();

	std::cout << "Pixels = " << LK << " , " << RK << std::endl;

	if (LK == RK)
	{
		std::cout << "           left == right! " << std::endl;
		return;
	}

	DuplicatorType::Pointer duplicator = DuplicatorType::New();
	duplicator->SetInputImage(filteredImg);
	duplicator->Update();
	AtlasImageType::Pointer outputImg_left = duplicator->GetOutput();

	//AtlasImageType::IndexType startIndex = { 600, 600, 600 };
	//AtlasImageType::IndexType endIndex = { 0, 0, 0 };

	itk::ImageRegionIterator<AtlasImageType> it_org(filteredImg, filteredImg->GetRequestedRegion());
	itk::ImageRegionIterator<AtlasImageType> it_out(outputImg_left, outputImg_left->GetRequestedRegion());

	for (it_org.GoToBegin(), it_out.GoToBegin(); !it_org.IsAtEnd() && !it_out.IsAtEnd(); ++it_org, ++it_out)
	{
		if (it_org.Value() == LK)
		{
			it_out.Set(it_org.Value());

			//AtlasImageType::IndexType index = it_org.GetIndex();
			//startIndex[0] = std::min(startIndex[0], index[0]);
			//startIndex[1] = std::min(startIndex[1], index[1]);
			//startIndex[2] = std::min(startIndex[2], index[2]);
			//endIndex[0] = std::max(endIndex[0], index[0]);
			//endIndex[1] = std::max(endIndex[1], index[1]);
			//endIndex[2] = std::max(endIndex[2], index[2]);
		}
		else
			it_out.Set(0);
	}

	// save index to write name 
	std::string writeName = this->m_fixedImgName + "_mask_left.nrrd";
	WriterFilterType::Pointer writer = WriterFilterType::New();
	try
	{
		writer->SetInput(outputImg_left);
		writer->SetFileName(writeName);
		writer->Update();
	}
	catch (itk::ExceptionObject& err)
	{
		std::cout << "Caught in writing registration file ! " << err.GetDescription() << std::endl;
		return;
	}

	std::cout << "--------------- Left Done --------------------------" << std::endl;

	/////////////////////////////////////////////////////////////////////////
	duplicator->SetInputImage(filteredImg);
	duplicator->Update();
	AtlasImageType::Pointer outputImg_right = duplicator->GetOutput();

	itk::ImageRegionIterator<AtlasImageType> it_org_right(filteredImg, filteredImg->GetRequestedRegion());
	itk::ImageRegionIterator<AtlasImageType> it_out_right(outputImg_right, outputImg_left->GetRequestedRegion());

	for (it_org_right.GoToBegin(), it_out_right.GoToBegin(); !it_org_right.IsAtEnd() && !it_out_right.IsAtEnd(); ++it_org_right, ++it_out_right)
	{
		if (it_org_right.Value() == RK)
		{
			it_out_right.Set(it_org_right.Value());
		}
		else
			it_out_right.Set(0);
	}

	// save index to write name 
	std::string writeName_right = this->m_fixedImgName + "_mask_right.nrrd";
	WriterFilterType::Pointer writer_right = WriterFilterType::New();

	try
	{
		writer_right->SetInput(outputImg_right);
		writer_right->SetFileName(writeName_right);
		writer_right->Update();
	}
	catch (itk::ExceptionObject& err)
	{
		std::cout << "Caught in writing registration file ! " << err.GetDescription() << std::endl;
		return;
	}

	std::cout << "--------------- Right Done --------------------------" << std::endl;

}


void SegmentationView::CropImageWithDefinedAxes()
{
	std::cout << "  start cropping image - " << this->m_fixedImgName << std::endl;

	if (!this->m_fixedImageForReg)
	{
		std::cout << " Please select one image to crop ... " << std::endl;
		return;
	}

	int lowerBounds[3] = { 0 };
	int upperBounds[3] = { 0 };

	lowerBounds[0] = m_Controls.cropX_lower->value();
	lowerBounds[1] = m_Controls.cropY_lower->value();
	lowerBounds[2] = m_Controls.cropZ_lower->value();

	upperBounds[0] = m_Controls.cropX_upper->value();
	upperBounds[1] = m_Controls.cropY_upper->value();
	upperBounds[2] = m_Controls.cropZ_upper->value();

	AtlasImageType::Pointer croppedImg = AtlasImageType::New();
	mitk::CastToItkImage(this->m_fixedImageForReg, croppedImg);

	// get name from Image Tree Node 
	int xDim = croppedImg->GetLargestPossibleRegion().GetSize()[0];
	int yDim = croppedImg->GetLargestPossibleRegion().GetSize()[1];
	int zDim = croppedImg->GetLargestPossibleRegion().GetSize()[2];

	if (upperBounds[0] >= xDim || upperBounds[0] < 1)
	{
		std::cout << " upper_0 exceed or forgot setting !" << std::endl;
		upperBounds[0] = xDim - 1;
	}

	if (upperBounds[1] >= yDim || upperBounds[1] < 1)
	{
		std::cout << " upper_1 exceed or forgot setting !" << std::endl;
		upperBounds[1] = yDim - 1;
	}

	if (upperBounds[2] >= zDim || upperBounds[2] < 1)
	{
		std::cout << " upper_2 exceed or forgot setting !" << std::endl;
		upperBounds[2] = zDim - 1;
	}

	AtlasImageType::IndexType startCropIndex = { lowerBounds[0], lowerBounds[1], lowerBounds[2] };
	AtlasImageType::SizeType cropSize;
	cropSize[0] = upperBounds[0] - lowerBounds[0] + 1;
	cropSize[1] = upperBounds[1] - lowerBounds[1] + 1;
	cropSize[2] = upperBounds[2] - lowerBounds[2] + 1;

	AtlasImageType::RegionType cropRegion(startCropIndex, cropSize);

	ExtractFilterType::Pointer extractor = ExtractFilterType::New();
	extractor->SetInput(croppedImg);
	extractor->SetExtractionRegion(cropRegion);
	extractor->SetDirectionCollapseToIdentity();
	extractor->Update();

	std::string writeName = "cp_" + this->m_fixedImgName + ".nrrd";

	this->addImageToNode(extractor->GetOutput(), writeName, true);


	WriterFilterType::Pointer writer = WriterFilterType::New();
	try{
		writer->SetInput(extractor->GetOutput());
		writer->SetFileName(writeName);
		writer->Update();
	}
	catch (itk::ExceptionObject& err)
	{
		std::cout << "Caught in writing registration file ! " << err.GetDescription() << std::endl;
		return;
	}

}


void SegmentationView::PolyDataReflectionX()
{
	// read poly meshes 
	QFileDialog fd;
	fd.setDirectory(m_recentDir.c_str());
	fd.setFileMode(QFileDialog::Directory);
	QString selected;

	if (fd.exec() == QDialog::Accepted)
	{
		m_recentDir = fd.directory().absolutePath().toAscii();
		QStringList myfiles = fd.selectedFiles();
		if (!myfiles.isEmpty())
			selected = myfiles[0];
	}

	QDir* meshDir = new QDir(selected);
	QStringList::Iterator it;
	QStringList files = meshDir->entryList();
	it = files.begin();

	int Cnt = 0;

	while (it != files.end())
	{
		if (QFileInfo(*meshDir, *it).isFile() && ((*it).endsWith)(".vtk"))
		{
			Cnt++;

			vtkSmartPointer< vtkPolyDataReader > pReader = vtkSmartPointer<vtkPolyDataReader>::New();

			pReader->SetFileName(QFileInfo(*meshDir, *it).absoluteFilePath().toAscii());
			std::string fileName = QFileInfo(*meshDir, *it).baseName().toStdString();
			pReader->Update();

			std::cout << " -- reflection the " << Cnt << "-th mesh : " << fileName << " : " << std::endl;

			vtkSmartPointer< vtkPolyData > loadedPoly = pReader->GetOutput();
			vtkSmartPointer< vtkPolyData > mirroredPoly = vtkSmartPointer< vtkPolyData >::New();
			mirroredPoly->DeepCopy(loadedPoly);

			vtkSmartPointer< vtkReflectionFilter > reflectionFilter = vtkSmartPointer< vtkReflectionFilter >::New();
			reflectionFilter->SetInputData(loadedPoly);
			reflectionFilter->SetPlaneToX();
			reflectionFilter->CopyInputOff();
			reflectionFilter->Update();

			mirroredPoly->DeepCopy(reflectionFilter->GetOutput());

			std::string outputName = fileName + "_mirrored_X.vtk";

			// save output
			vtkSmartPointer< vtkPolyDataWriter > writer = vtkSmartPointer< vtkPolyDataWriter >::New();
			writer->SetInputData(mirroredPoly);
			writer->SetFileName(outputName.c_str());
			writer->Update();

		}
		++it;
	}

	delete meshDir;

	std::cout << " ---- reflection Done ---- " << std::endl;
}


void SegmentationView::PolyDataReflectionY()
{
	// read poly meshes 
	QFileDialog fd;
	fd.setDirectory(m_recentDir.c_str());
	fd.setFileMode(QFileDialog::Directory);
	QString selected;

	if (fd.exec() == QDialog::Accepted)
	{
		m_recentDir = fd.directory().absolutePath().toAscii();
		QStringList myfiles = fd.selectedFiles();
		if (!myfiles.isEmpty())
			selected = myfiles[0];
	}

	QDir* meshDir = new QDir(selected);
	QStringList::Iterator it;
	QStringList files = meshDir->entryList();
	it = files.begin();

	int Cnt = 0;

	while (it != files.end())
	{
		if (QFileInfo(*meshDir, *it).isFile() && ((*it).endsWith)(".vtk"))
		{
			Cnt++;

			vtkSmartPointer< vtkPolyDataReader > pReader = vtkSmartPointer<vtkPolyDataReader>::New();

			pReader->SetFileName(QFileInfo(*meshDir, *it).absoluteFilePath().toAscii());
			std::string fileName = QFileInfo(*meshDir, *it).baseName().toStdString();
			pReader->Update();

			std::cout << " -- reflection the " << Cnt << "-th mesh : " << fileName << " : " << std::endl;

			vtkSmartPointer< vtkPolyData > loadedPoly = pReader->GetOutput();
			vtkSmartPointer< vtkPolyData > mirroredPoly = vtkSmartPointer< vtkPolyData >::New();
			mirroredPoly->DeepCopy(loadedPoly);

			vtkSmartPointer< vtkReflectionFilter > reflectionFilter = vtkSmartPointer< vtkReflectionFilter >::New();
			reflectionFilter->SetInputData(loadedPoly);
//			reflectionFilter->SetPlaneToY();
			reflectionFilter->SetPlaneToYMin();
			reflectionFilter->CopyInputOff();
			reflectionFilter->Update();

			mirroredPoly->DeepCopy(reflectionFilter->GetOutput());

			std::string outputName = fileName + "_mirrored_Y.vtk";

			// save output
			vtkSmartPointer< vtkPolyDataWriter > writer = vtkSmartPointer< vtkPolyDataWriter >::New();
			writer->SetInputData(mirroredPoly);
			writer->SetFileName(outputName.c_str());
			writer->Update();

		}
		++it;
	}

	delete meshDir;

	std::cout << " ---- reflection Done ---- " << std::endl;
}


void SegmentationView::PolyDataReflectionZ()
{
	std::vector< vtkSmartPointer<vtkPolyData> > datasets; 

	// read poly meshes 
	QFileDialog fd;
	fd.setDirectory(m_recentDir.c_str());
	fd.setFileMode(QFileDialog::Directory);
	QString selected;

	if (fd.exec() == QDialog::Accepted)
	{
		m_recentDir = fd.directory().absolutePath().toAscii();
		QStringList myfiles = fd.selectedFiles();
		if (!myfiles.isEmpty())
			selected = myfiles[0];
	}

	QDir* meshDir = new QDir(selected);
	QStringList::Iterator it;
	QStringList files = meshDir->entryList();
	it = files.begin();

	int Cnt = 0;

	while (it != files.end())
	{
		if (QFileInfo(*meshDir, *it).isFile() && ((*it).endsWith)(".vtk"))
		{
			Cnt++;

			vtkSmartPointer< vtkPolyDataReader > pReader = vtkSmartPointer<vtkPolyDataReader>::New();

			pReader->SetFileName(QFileInfo(*meshDir, *it).absoluteFilePath().toAscii());
			std::string fileName = QFileInfo(*meshDir, *it).baseName().toStdString();
			pReader->Update();

			vtkSmartPointer< vtkPolyData > readPoly = pReader->GetOutput(); 

			datasets.push_back(readPoly); 

			//std::cout << " -- reflection the " << Cnt << "-th mesh : " << fileName << " : " << std::endl;
			//vtkSmartPointer< vtkPolyData > loadedPoly = pReader->GetOutput();
			//vtkSmartPointer< vtkPolyData > mirroredPoly = vtkSmartPointer< vtkPolyData >::New();
			//mirroredPoly->DeepCopy(loadedPoly);
			//vtkSmartPointer< vtkReflectionFilter > reflectionFilter = vtkSmartPointer< vtkReflectionFilter >::New();
			//reflectionFilter->SetInputData(loadedPoly);
			//reflectionFilter->SetPlaneToZ();
			//reflectionFilter->CopyInputOff();
			//reflectionFilter->Update();
			//mirroredPoly->DeepCopy(reflectionFilter->GetOutput());
			//std::string outputName = fileName + "_mirrored_Z.vtk";
			//// save output
			//vtkSmartPointer< vtkPolyDataWriter > writer = vtkSmartPointer< vtkPolyDataWriter >::New();
			//writer->SetInputData(mirroredPoly);
			//writer->SetFileName(outputName.c_str());
			//writer->Update();

		}
		++it;
	}

	delete meshDir;

	this->CenterOfMassToOrigin(datasets, true); 

	this->visualizeDatasets(datasets, "center", false); 

	// landmark transformation 
	vtkSmartPointer< vtkPointLocator > locator = vtkPointLocator::New(); 
	locator->SetDataSet(datasets.at(0)); 
	locator->BuildLocator(); 

	vtkPoints* points = datasets.at(1)->GetPoints(); 

	for (int ptID = 0; ptID < points->GetNumberOfPoints(); ptID++)
	{
		double* pnt = points->GetPoint(ptID); 
		vtkIdType id = locator->FindClosestPoint(pnt);


	}




	std::cout << " ---- reflection Done ---- " << std::endl;
}


void SegmentationView::SeparateMasks()
{
	// import masks 
	QFileDialog fd;
	fd.setDirectory(m_recentDir.c_str());
	fd.setFileMode(QFileDialog::Directory);
	QString selected;

	if (fd.exec() == QDialog::Accepted)
	{
		m_recentDir = fd.directory().absolutePath().toAscii();
		QStringList myfiles = fd.selectedFiles();
		if (!myfiles.isEmpty())
			selected = myfiles[0];
	}

	QDir* imageDir = new QDir(selected);
	QStringList::Iterator it;
	QStringList files = imageDir->entryList();

	it = files.begin();

	typedef itk::ImageDuplicator< AtlasImageType > ClonerType; 

	while (it != files.end())
	{
		if (QFileInfo(*imageDir, *it).isFile() && ((*it).endsWith)(".pic") || (*it).endsWith(".nrrd") || (*it).endsWith(".mhd"))
		{
			mitk::DataNodeFactory::Pointer nodeReader = mitk::DataNodeFactory::New();

			QFileInfo& fileInfoIMG_ = QFileInfo(*imageDir, *it);
			std::string imgFileName = fileInfoIMG_.absoluteFilePath().toStdString();
			std::string imgAbsName = QFileInfo(*imageDir, *it).baseName().toStdString();

			AtlasImageType::Pointer mask = AtlasImageType::New();

			try{
				nodeReader->SetFileName(imgFileName);
				nodeReader->Update();

				mitk::DataNode::Pointer imageNode = nodeReader->GetOutput();
				mitk::Image::Pointer mitkImage = dynamic_cast<mitk::Image*>(imageNode->GetData());
				mitk::CastToItkImage(mitkImage, mask);
			}
			catch (...){}

			int xDim = mask->GetLargestPossibleRegion().GetSize()[0];
			int yDim = mask->GetLargestPossibleRegion().GetSize()[1];
			int zDim = mask->GetLargestPossibleRegion().GetSize()[2];

			AtlasImageType::SizeType size = { xDim, yDim, zDim }; 
			AtlasImageType::IndexType index = { 0, 0, 0 }; 
			AtlasImageType::RegionType region(index, size); 

			AtlasImageType::Pointer left = AtlasImageType::New(); 
			left->SetRegions(region); 
			left->SetOrigin(mask->GetOrigin()); 
			left->SetSpacing(mask->GetSpacing()); 
			left->Allocate(); 
			left->FillBuffer(itk::NumericTraits< AtlasImageType::PixelType >::Zero);


			AtlasImageType::Pointer right = AtlasImageType::New();
			right->SetRegions(region);
			right->SetOrigin(mask->GetOrigin());
			right->SetSpacing(mask->GetSpacing());
			right->Allocate();
			right->FillBuffer(itk::NumericTraits< AtlasImageType::PixelType >::Zero);
			

			for (int x = 0; x < xDim / 2; x++)
			{
				for (int y = 0; y < yDim; y++)
				{
					for (int z = 0; z < zDim; z++)
					{
						AtlasImageType::IndexType index = { x, y, z }; 
						if (mask->GetPixel(index) > 0)
							right->SetPixel(index, 1); 
					}
				}
			}

			for (int x = xDim / 2; x < xDim ; x++)
			{
				for (int y = 0; y < yDim; y++)
				{
					for (int z = 0; z < zDim; z++)
					{
						AtlasImageType::IndexType index = { x, y, z };
						if (mask->GetPixel(index) > 0)
							left->SetPixel(index, 1);
					}
				}
			}

			this->addImageToNode(left, imgAbsName + "_left", false); 
			this->addImageToNode(right, imgAbsName + "_right", false); 

			this->writeAtlasImage(left, imgAbsName + "_left.nrrd"); 
			this->writeAtlasImage(right, imgAbsName + "_right.nrrd"); 
		
		}

		++it;
	}

	delete imageDir;

	std::cout << "--------------- Separation Done --------------------------" << std::endl;
}


void SegmentationView::NormalizeMask()
{
	typedef itk::ChangeInformationImageFilter< AtlasImageType > FilterType; 

	QFileDialog fd;
	fd.setDirectory(m_recentDir.c_str());
	fd.setFileMode(QFileDialog::Directory);
	QString selected;

	if (fd.exec() == QDialog::Accepted)
	{
		m_recentDir = fd.directory().absolutePath().toAscii();

		QStringList myfiles = fd.selectedFiles();
		if (!myfiles.isEmpty())
			selected = myfiles[0];
	}

	QDir* imageDir = new QDir(selected);

	QStringList::Iterator it;
	QStringList files = imageDir->entryList();
	it = files.begin();

	while (it != files.end())
	{
		if (QFileInfo(*imageDir, *it).isFile() && ((*it).endsWith)(".pic") || (*it).endsWith(".nrrd") || (*it).endsWith(".mhd"))
		{
			mitk::DataNodeFactory::Pointer nodeReader = mitk::DataNodeFactory::New();

			QFileInfo& fileInfoIMG_ = QFileInfo(*imageDir, *it);
			std::string imgFileName = fileInfoIMG_.absoluteFilePath().toStdString();
			std::string imgAbsName = QFileInfo(*imageDir, *it).baseName().toStdString();

			AtlasImageType::Pointer itkImage = AtlasImageType::New();

			try{
				nodeReader->SetFileName(imgFileName);
				nodeReader->Update();

				mitk::DataNode::Pointer imageNode = nodeReader->GetOutput();
				mitk::Image::Pointer mitkImage = dynamic_cast<mitk::Image*>(imageNode->GetData());
				mitk::CastToItkImage(mitkImage, itkImage);
			}
			catch (...){}	

			itk::ImageRegionIterator< AtlasImageType > it(itkImage, itkImage->GetRequestedRegion()); 
			AtlasImageType::IndexType minIndex = { 512, 512, 512 }; 
			AtlasImageType::IndexType maxIndex = { 0, 0, 0 }; 

			for (it.GoToBegin(); !it.IsAtEnd(); ++it)
			{
				if (it.Value() > 0)
				{
					AtlasImageType::IndexType currIndex = it.GetIndex(); 

					if (currIndex[0] < minIndex[0]) minIndex[0] = currIndex[0]; 
					if (currIndex[0] > maxIndex[0]) maxIndex[0] = currIndex[0]; 

					if (currIndex[1] < minIndex[1]) minIndex[1] = currIndex[1]; 
					if (currIndex[1] > maxIndex[1]) maxIndex[1] = currIndex[1];

					if (currIndex[2] < minIndex[2]) minIndex[2] = currIndex[2]; 
					if (currIndex[2] > maxIndex[2]) maxIndex[2] = currIndex[2]; 
				}
			}

			std::cout << " current Image " << imgAbsName << " in the region : [" << minIndex[0] << " - " << maxIndex[0] << "] , ["
				<< minIndex[1] << " - " << maxIndex[1] << "] , ["
				<< minIndex[2] << " - " << maxIndex[2] << "] " << std::endl;


			// comment - 1 Image Origin Correlation 
			// comment - 2 Mask Normalization 
		/*	AtlasImageType::Pointer itkCorrectImg = AtlasImageType::New(); 
			AtlasImageType::IndexType startIndex = { 0, 0, 0 }; 
			AtlasImageType::SizeType size; 
			size[0] = itkImage->GetLargestPossibleRegion().GetSize()[0]; 
			size[1] = itkImage->GetLargestPossibleRegion().GetSize()[1]; 
			size[2] = itkImage->GetLargestPossibleRegion().GetSize()[2]; 
			AtlasImageType::RegionType region(startIndex, size); 
			itkCorrectImg->SetRegions(region); 
			itkCorrectImg->SetSpacing(itkImage->GetSpacing()); 
			itkCorrectImg->Allocate(); 
			itkCorrectImg->FillBuffer(0); 

			for (int x = 0; x < size[0]; x++)
			{
				for (int y = 0; y < size[1]; y++)
				{
					for (int z = 0; z < size[2]; z++)
					{
						AtlasImageType::IndexType origIndex = { x, size[1] - 1 - y, size[2] - 1 - z };
						AtlasImageType::IndexType currIndex = { x, y, z }; 

						AtlasImageType::PixelType value = itkImage->GetPixel(origIndex); 
						itkCorrectImg->SetPixel(currIndex, value); 
					}
				}
			}

			itkCorrectImg->Modified(); */
			/*FilterType::Pointer filter = FilterType::New(); 
			filter->SetInput(itkImage); 
			filter->ChangeSpacingOff(); 
			filter->ChangeDirectionOff(); 

			AtlasImageType::PointType origin; 
			origin[0] = 0; 
			origin[1] = 0; 
			origin[2] = 0; 

			filter->SetOutputOrigin(origin); 
			filter->ChangeOriginOn(); 
			filter->CenterImageOn(); 

			filter->UpdateOutputInformation(); 
			filter->UpdateLargestPossibleRegion();
			filter->Update(); */

			// this->writeAtlasImage(itkCorrectImg, imgAbsName + ".nrrd");
		}

		++it;
	}

	delete imageDir;

	std::cout << "--------------- Normalization Done --------------------------" << std::endl;
}


void SegmentationView::CheckMaskBoundary()
{
	std::cout << __FUNCTION__ << " : checking boundary size and index ... " << std::endl; 

	// import masks 
	QFileDialog fd;
	fd.setDirectory(m_recentDir.c_str());
	fd.setFileMode(QFileDialog::Directory);
	QString selected;

	if (fd.exec() == QDialog::Accepted)
	{
		m_recentDir = fd.directory().absolutePath().toAscii();
		QStringList myfiles = fd.selectedFiles();
		if (!myfiles.isEmpty())
			selected = myfiles[0];
	}

	QDir* imageDir = new QDir(selected);
	QStringList::Iterator it;
	QStringList files = imageDir->entryList();

	it = files.begin();

	while (it != files.end())
	{
		if (QFileInfo(*imageDir, *it).isFile() && ((*it).endsWith)(".nii.gz") || (*it).endsWith(".nrrd") || (*it).endsWith(".mhd"))
		{
			mitk::DataNodeFactory::Pointer nodeReader = mitk::DataNodeFactory::New();

			QFileInfo& fileInfoIMG_ = QFileInfo(*imageDir, *it);
			std::string imgFileName = fileInfoIMG_.absoluteFilePath().toStdString();
			std::string imgAbsName = QFileInfo(*imageDir, *it).baseName().toStdString();

			AtlasImageType::Pointer mask = AtlasImageType::New();

			try{
				nodeReader->SetFileName(imgFileName);
				nodeReader->Update();

				mitk::DataNode::Pointer imageNode = nodeReader->GetOutput();
				mitk::Image::Pointer mitkImage = dynamic_cast<mitk::Image*>(imageNode->GetData());
				mitk::CastToItkImage(mitkImage, mask);
			}
			catch (...){}

			int xDim = mask->GetLargestPossibleRegion().GetSize()[0];
			int yDim = mask->GetLargestPossibleRegion().GetSize()[1];
			int zDim = mask->GetLargestPossibleRegion().GetSize()[2];

			AtlasImageType::IndexType maxIndex = { 0, 0, 0 }; 
			AtlasImageType::IndexType minIndex = { xDim, yDim, zDim }; 
			AtlasImageType::IndexType currIndex = { 0, 0, 0 }; 

			itk::ImageRegionIterator< AtlasImageType > iter(mask, mask->GetRequestedRegion()); 
			for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
			{
				currIndex = iter.GetIndex(); 

				if (iter.Get() > 0)
				{
					if (currIndex[0] < minIndex[0]) minIndex[0] = currIndex[0];
					if (currIndex[1] < minIndex[1]) minIndex[1] = currIndex[1]; 
					if (currIndex[2] < minIndex[2]) minIndex[2] = currIndex[2]; 

					if (currIndex[0] > maxIndex[0]) maxIndex[0] = currIndex[0]; 
					if (currIndex[1] > maxIndex[1]) maxIndex[1] = currIndex[1]; 
					if (currIndex[2] > maxIndex[2]) maxIndex[2] = currIndex[2]; 
				}
			}

			cout << __FUNCTION__ << " : " << imgAbsName << " [ " << xDim << " - " << yDim << " - " << zDim << " ] " ;
			cout << " ROI = X [ " << minIndex[0] << " - " << maxIndex[0] << " ] = " << maxIndex[0] - minIndex[0]; 
			cout << " Y [ " << minIndex[1] << " - " << maxIndex[1] << " ] = " << maxIndex[1] - minIndex[1]; 
			cout << " Z [ " << minIndex[2] << " - " << maxIndex[2] << " ] = " << maxIndex[2] - minIndex[2] << endl; 

		}

		++it;
	}

	delete imageDir;

	std::cout << "--------------- Done --------------------------" << std::endl;
}


void SegmentationView::ComputeAndCropROI()
{
	if (!this->m_segImage)
	{
		MITK_WARN << " Please Select an Image to crop ... " << endl; 
		return; 
	}

	if (!this->m_probMap)
	{
		MITK_WARN << " Please Select an Mask to crop ... " << endl; 
	}

	int padding = m_Controls.closeRadius->value(); 
	std::cout << __FUNCTION__ << " cropped padding = " << padding ;

	AtlasImageType::IndexType minIndex = { this->m_segImage->GetLargestPossibleRegion().GetSize()[0], this->m_segImage->GetLargestPossibleRegion().GetSize()[1], this->m_segImage->GetLargestPossibleRegion().GetSize()[2] };
	AtlasImageType::IndexType maxIndex = { 0, 0, 0 }; 
	
	itk::ImageRegionIterator< AtlasImageType > iter(this->m_probMap, this->m_probMap->GetRequestedRegion()); 
	for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
	{
		if (iter.Get() > 0)
		{
			AtlasImageType::IndexType currIndex = iter.GetIndex();

			if (currIndex[0] < minIndex[0]) minIndex[0] = currIndex[0]; 
			if (currIndex[1] < minIndex[1]) minIndex[1] = currIndex[1]; 
			if (currIndex[2] < minIndex[2]) minIndex[2] = currIndex[2]; 

			if (currIndex[0] > maxIndex[0]) maxIndex[0] = currIndex[0]; 
			if (currIndex[1] > maxIndex[1]) maxIndex[1] = currIndex[1]; 
			if (currIndex[2] > maxIndex[2]) maxIndex[2] = currIndex[2]; 
		}
	}

	cout << " cropped size = " << maxIndex[0] - minIndex[0] << " - " << maxIndex[1] - minIndex[1] << " - " << maxIndex[2] - minIndex[2] << endl; 

	//////////////////////////////////////////////////////////////////////////////

	AtlasImageType::IndexType startIndex; 
	for (int i = 0; i < 3; i++)
		startIndex[i] = std::max(0, (int)(minIndex[i] - padding)); 

	AtlasImageType::SizeType cropSize; 
	for (int i = 0; i < 3; i++)
	{
		if (maxIndex[i] + padding - startIndex[i] >= this->m_segImage->GetLargestPossibleRegion().GetSize()[i] - startIndex[i])
			cropSize[i] = this->m_segImage->GetLargestPossibleRegion().GetSize()[i] - startIndex[i]; 
		else
			cropSize[i] = std::min((int)(maxIndex[i] + padding - startIndex[i] + 1), (int)this->m_segImage->GetLargestPossibleRegion().GetSize()[i]);
	}
		

	AtlasImageType::RegionType cropRegion(startIndex, cropSize);

	ExtractFilterType::Pointer extractor = ExtractFilterType::New();
	extractor->SetInput(this->m_segImage);
	extractor->SetExtractionRegion(cropRegion);
	extractor->SetDirectionCollapseToIdentity();
	extractor->Update();

	std::string writeName = this->m_segImageName +  "_cp.nrrd";

	WriterFilterType::Pointer writer = WriterFilterType::New();
	try{
		writer->SetInput(extractor->GetOutput());
		writer->SetFileName(writeName);
		writer->Update();
	}
	catch (itk::ExceptionObject& err)
	{
		std::cout << "Caught in writing registration file ! " << err.GetDescription() << std::endl;
		return;
	}

	///////////////////////////////////////////////////////////

	extractor->SetInput(this->m_probMap);
	extractor->SetExtractionRegion(cropRegion);
	extractor->SetDirectionCollapseToIdentity();
	extractor->Update();

	writeName = this->m_segImageName + "_cp_mask.nrrd";

	try{
		writer->SetInput(extractor->GetOutput());
		writer->SetFileName(writeName);
		writer->Update();
	}
	catch (itk::ExceptionObject& err)
	{
		std::cout << "Caught in writing registration file ! " << err.GetDescription() << std::endl;
		return;
	}

	cout << __FUNCTION__ << " Done cropping ....." << endl; 
}


/////////////////////////// --- Deep Learning Utilities ---- //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void SegmentationView::generateSlicesWithMode(std::string fileName, std::string sliceMode, Image3DType::Pointer itkImage)
{	
	// just write all jpg images to the folder and move manually ! 
	int xDim = itkImage->GetLargestPossibleRegion().GetSize()[0]; 
	int yDim = itkImage->GetLargestPossibleRegion().GetSize()[1]; 
	int zDim = itkImage->GetLargestPossibleRegion().GetSize()[2]; 

	if (m_Controls.selectSliceMode->currentText() == "")
	{
		sliceMode = "SC"; // be default 
	}

	typedef itk::ExtractImageFilter<Image3DType, Image2DType> filter2dType;
	typedef itk::ImageFileWriter<Image2DType> write2DImageType;

	filter2dType::Pointer filter2D = filter2dType::New();
	write2DImageType::Pointer writer2D = write2DImageType::New();

	// remove the "_mask"
	std::string subString = "_mask";
	std::string::size_type foundPos = fileName.find(subString);

	while (foundPos != std::string::npos)
	{
		fileName.erase(foundPos, subString.length());
		foundPos = fileName.find(subString, foundPos);
	}

	std::string baseName = fileName + sliceMode; 
	
	if (m_Controls.selectSliceMode->currentText() == "SC") // slice from A (z-axes)
	{
		Image3DType::RegionType desiredRegion;
		Image3DType::SizeType desiredSize;
		Image3DType::IndexType desiredIndex;
		desiredIndex[0] = 0;
		desiredIndex[1] = 0;
		desiredSize[0] = xDim;
		desiredSize[1] = yDim;
		desiredSize[2] = 0; 

		for (int z = 0; z < zDim; z++)
		{
			std::stringstream ss;
			ss << std::setfill('0') << setw(3) << z + 1; 
			std::string fileName = baseName + ss.str(); 
			if (m_Controls.check_mask->isChecked())
			{
				fileName = fileName + "_mask.jpg"; 
			}
			else
			{
				fileName = fileName + ".jpg";
			}

			desiredIndex[2] = z;
			desiredRegion.SetIndex(desiredIndex);
			desiredRegion.SetSize(desiredSize);


			filter2D->InPlaceOn();
			filter2D->SetDirectionCollapseToIdentity();
			filter2D->SetInput(itkImage);
			filter2D->SetExtractionRegion(desiredRegion);
			
			try{
				writer2D->SetInput(filter2D->GetOutput());
				writer2D->SetFileName(fileName);
				writer2D->Update();
			}
			catch (itk::ExceptionObject &err)
			{
				MITK_WARN << err.GetDescription();
				return;
			}		
		}
	}
	else if (m_Controls.selectSliceMode->currentText() == "SA") // slice from C (y-axes)
	{
		Image3DType::RegionType desiredRegion;
		Image3DType::SizeType desiredSize;
		Image3DType::IndexType desiredIndex;
		desiredIndex[0] = 0;
		desiredIndex[2] = 0;
		desiredSize[0] = xDim;
		desiredSize[1] = 0;
		desiredSize[2] = zDim;

		for (int y = 0; y < yDim; y++)
		{
			std::stringstream ss;
			ss << std::setfill('0') << setw(3) << y + 1;
			std::string fileName = baseName + ss.str();
			if (m_Controls.check_mask->isChecked() )
			{
				fileName = fileName + "_mask.jpg";
			}
			else fileName = fileName + ".jpg";

			desiredIndex[1] = y;
			desiredRegion.SetIndex(desiredIndex);
			desiredRegion.SetSize(desiredSize);

			filter2D->InPlaceOn();
			filter2D->SetDirectionCollapseToIdentity();
			filter2D->SetInput(itkImage);
			filter2D->SetExtractionRegion(desiredRegion);

			writer2D->SetInput(filter2D->GetOutput());
			writer2D->SetFileName(fileName);
			writer2D->Update();
		}
	}
	else if (m_Controls.selectSliceMode->currentText() == "CA") // slice from S (x-axes)
	{
		Image3DType::RegionType desiredRegion;
		Image3DType::SizeType desiredSize;
		Image3DType::IndexType desiredIndex;
		desiredIndex[1] = 0;
		desiredIndex[2] = 0;
		desiredSize[0] = 0;
		desiredSize[1] = yDim;
		desiredSize[2] = zDim;

		for (int x = 0; x < xDim; x++)
		{
			std::stringstream ss;
			ss << std::setfill('0') << setw(3) << x + 1;
			std::string fileName = baseName + ss.str();
			if (m_Controls.check_mask->isChecked())
			{
				fileName = fileName + "_mask.jpg";
			}
			else fileName = fileName + ".jpg";

			desiredIndex[0] = x;
			desiredRegion.SetIndex(desiredIndex);
			desiredRegion.SetSize(desiredSize);

			filter2D->InPlaceOn();
			filter2D->SetDirectionCollapseToIdentity();
			filter2D->SetInput(itkImage);
			filter2D->SetExtractionRegion(desiredRegion);

			writer2D->SetInput(filter2D->GetOutput());
			writer2D->SetFileName(fileName);
			writer2D->Update();
		}
	}
	else
	{
		MITK_WARN << " Please select a slice mode !" << std::endl;
		return;
	}

	
}


/////////////////////////// --- Segmentation Utilities ---- ///////////////////////////////////////////////////////////////////////


void SegmentationView::loadProbShapeForLocalization()
{
	if (!this->m_probShape) this->m_probShape = vtkPolyData::New(); 

	std::string probShapeName = this->m_probMapName + ".vtk";

	QString qFilename = "J:/IGDMedApp_VS12/build/bin/Release/Prob Vtk/" + QString::fromStdString(probShapeName); 
	QByteArray byteArr = qFilename.toAscii();
	const char * filename = byteArr.data();  //QByteArray( qFilename ).data();
	
	mitk::DataNodeFactory::Pointer nodeReader = mitk::DataNodeFactory::New();
	try
	{
		nodeReader->SetFileName(filename);
		nodeReader->Update();
		mitk::DataNode::Pointer loadedSurfaceNode = nodeReader->GetOutput(); // loaded model from button "Load Organ Model"

		if (loadedSurfaceNode)
		{
			mitk::Surface::Pointer surface = dynamic_cast<mitk::Surface*>(loadedSurfaceNode->GetData());
			if (surface)
			{
				std::cout << __FUNCTION__ << " probability map is loaded : " << probShapeName << std::endl;
				this->m_probShape->DeepCopy(surface->GetVtkPolyData());
			}
		}
	}
	catch (...)
	{
		fprintf(stderr, "Could not open file %s \n\n", filename);
		return;
	}

	this->m_probShape->Modified(); 

}


void SegmentationView::loadDeformShapeFromLocal()
{
	std::cout << __FUNCTION__ ;

	if (!this->m_deformedShape) m_deformedShape = vtkSmartPointer< vtkPolyData >::New();

	mitk::DataNodeFactory::Pointer nodeReader = mitk::DataNodeFactory::New();

	QString qFilename = "J:/IGDMedApp_VS12/build/bin/Release/Models/" + m_Controls.m_modelNameSelect->currentText() + ".vtk";
	QByteArray byteArr = qFilename.toAscii();
	const char * filename = byteArr.data();  //QByteArray( qFilename ).data();
	try
	{
		nodeReader->SetFileName(filename);
		nodeReader->Update();
		mitk::DataNode::Pointer loadedSurfaceNode = nodeReader->GetOutput(); // loaded model from button "Load Organ Model"

		if (loadedSurfaceNode)
		{
			mitk::Surface::Pointer surface = dynamic_cast<mitk::Surface*>(loadedSurfaceNode->GetData());
			if (surface)
			{
				std::cout << " deformShape is loaded from " << filename << std::endl;
				m_deformedShape->DeepCopy(surface->GetVtkPolyData());
			}
		}
	}
	catch (...)
	{
		fprintf(stderr, "Could not open file %s \n\n", filename);
		return;
	}

	this->m_deformedShape->Modified(); 

}


void SegmentationView::loadGroundTruthShapeForPreTraining()
{
	if (!this->m_preTrainShape) this->m_preTrainShape = vtkPolyData::New();

	std::string GTShapeName = this->m_probMapName + ".vtk";

	QString qFilename = "J:/IGDMedApp_VS12/build/bin/Release/Train Vtk/" + QString::fromStdString(GTShapeName);
	QByteArray byteArr = qFilename.toAscii();
	const char * filename = byteArr.data();  //QByteArray( qFilename ).data();

	mitk::DataNodeFactory::Pointer nodeReader = mitk::DataNodeFactory::New();
	try
	{
		nodeReader->SetFileName(filename);
		nodeReader->Update();
		mitk::DataNode::Pointer loadedSurfaceNode = nodeReader->GetOutput(); // loaded model from button "Load Organ Model"

		if (loadedSurfaceNode)
		{
			mitk::Surface::Pointer surface = dynamic_cast<mitk::Surface*>(loadedSurfaceNode->GetData());
			if (surface)
			{
				this->m_preTrainShape->DeepCopy(surface->GetVtkPolyData());
			}
		}
	}
	catch (...)
	{
		fprintf(stderr, "Could not open file %s \n\n", filename);
		return;
	}

	this->m_preTrainShape->Modified();

}


void SegmentationView::alignDeformShapeWithProbShape(vtkPolyData* _sourcePoly, vtkPolyData* _targetPoly)
{
	cout << __FUNCTION__ << endl; 

	vtkSmartPointer< vtkIterativeClosestPointTransform > icp = vtkSmartPointer< vtkIterativeClosestPointTransform >::New();
	icp->SetSource(_sourcePoly);
	icp->SetTarget(_targetPoly);  // m_probShape
//	icp->GetLandmarkTransform()->SetModeToRigidBody();
	icp->GetLandmarkTransform()->SetModeToSimilarity();  
	icp->GetLandmarkTransform()->Modified();
	icp->StartByMatchingCentroidsOn();
	icp->SetMaximumNumberOfIterations(50);
	icp->Modified();
	icp->Update();

	vtkSmartPointer<vtkTransformPolyDataFilter> icpTransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	icpTransformFilter->SetInputData(_sourcePoly);
	icpTransformFilter->SetTransform(icp);
	icpTransformFilter->Update();

	_sourcePoly->DeepCopy(icpTransformFilter->GetOutput());
	_sourcePoly->Modified();
}


void SegmentationView::transformDeformShapeToProbMap(vtkPolyData* _inputPoly)
{	
	cout << __FUNCTION__ << endl; 

	vtkPolyData* probCopy = vtkPolyData::New(); 
	probCopy->DeepCopy(this->m_probShape); 

	this->reComputeNormals(_inputPoly);
//	this->reComputeNormals(probCopy); 

	vtkSmartPointer<vtkPolyDataNormals> normalsGen = vtkSmartPointer<vtkPolyDataNormals>::New();

	normalsGen->SplittingOff();
	normalsGen->ConsistencyOff();
	normalsGen->AutoOrientNormalsOff();
	normalsGen->ComputePointNormalsOn(); 
	normalsGen->ComputeCellNormalsOn(); 
	normalsGen->SetFeatureAngle(180);
	normalsGen->SetInputData(probCopy); 
	normalsGen->Update(); 
	probCopy->DeepCopy(normalsGen->GetOutput()); 
	probCopy->Modified(); 

	vtkSmartPointer< vtkDataArray > weights = vtkSmartPointer< vtkDoubleArray >::New();
	weights = _inputPoly->GetPointData()->GetArray("weights"); 


	vtkDataArray* normalDeform = _inputPoly->GetPointData()->GetNormals();
	vtkDataArray* normalProb = probCopy->GetPointData()->GetNormals(); 
	vtkDataArray* normalCellProb = probCopy->GetCellData()->GetNormals(); 

	vtkSmartPointer< vtkCellLocator > cellLocator = vtkSmartPointer< vtkCellLocator >::New(); 
	cellLocator->SetDataSet(probCopy);
	cellLocator->SetNumberOfCellsPerBucket(1); 
	cellLocator->BuildLocator(); 
	cellLocator->Update(); 

	vtkPoints* outputPoints = _inputPoly->GetPoints();

	for (int ptId = 0; ptId < _inputPoly->GetNumberOfPoints(); ptId++)
	{

		double point[3]; 
		_inputPoly->GetPoint(ptId, point);

//>>>>>> Find Closest Point >>>>>>>>>>>>>>
		
		double closestPoint[3]; 
		double closestDist; 
		int subId;
		vtkIdType cellId; 
		vtkIdType pointId; 

 		cellLocator->FindClosestPoint(point, closestPoint, cellId, subId, closestDist); 

//>>>>>> Condition Judge >>>>>>>>>>>>>>>>>>

		double normalCell[3];
		normalCellProb->GetTuple(cellId, normalCell); 

		double normalPoint[3];
		normalDeform->GetTuple(ptId, normalPoint);

		int Count = 0; 
		for (int j = 0; j < 3; j++)
		{
			if (normalCell[j] * normalPoint[j] >= 0)
				Count++; 
		}

		if (Count >= 2)
		{
			outputPoints->SetPoint(ptId, closestPoint[0], closestPoint[1], closestPoint[2]); 
			weights->SetTuple1(ptId, 1.0); 
		}
		else
		{
			outputPoints->SetPoint(ptId, point[0], point[1], point[2]);
			weights->SetTuple1(ptId, 0.0);
		}
		
	}

	weights->Modified(); 
	outputPoints->Modified(); 
	_inputPoly->GetPointData()->SetScalars(weights);
	_inputPoly->Modified();

	vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
	smoothFilter->SetInputData(_inputPoly); 
	smoothFilter->SetNumberOfIterations(15);
	smoothFilter->SetRelaxationFactor(0.1);
	smoothFilter->FeatureEdgeSmoothingOff();
	smoothFilter->BoundarySmoothingOn();
	smoothFilter->Update();

	_inputPoly->DeepCopy(smoothFilter->GetOutput()); 
	_inputPoly->Modified(); 

	this->reComputeNormals(_inputPoly); 

}


void SegmentationView::reComputeNormals(vtkPolyData* _poly)
{
	vtkSmartPointer<vtkPolyDataNormals> normalsGen = vtkSmartPointer<vtkPolyDataNormals>::New();

	normalsGen->SplittingOff();
	normalsGen->ConsistencyOff();
//	normalsGen->AutoOrientNormalsOff();
	normalsGen->AutoOrientNormalsOn(); 
	normalsGen->ComputePointNormalsOn();
	normalsGen->SetFeatureAngle(180);
	
	//	normalsGen->SplittingOn();  // preserve sharp edges
//	normalsGen->ConsistencyOn(); // polygon order 
////	normalsGen->ConsistencyOff();
//	normalsGen->AutoOrientNormalsOn(); 
////	normalsGen->AutoOrientNormalsOff();
//	normalsGen->ComputePointNormalsOn();
//	normalsGen->ComputeCellNormalsOff(); 
//	normalsGen->SetFeatureAngle(180);
////	normalsGen->SetSplitting(0); 

	normalsGen->SetInputData(_poly);
	normalsGen->Update();

	_poly->DeepCopy(normalsGen->GetOutput());
	_poly->Modified();
}


mitk::Surface::Pointer SegmentationView::GetTransformedGeometrySurface(mitk::Surface::Pointer _surface)
{

	vtkSmartPointer<vtkPolyData> polyData = _surface->GetVtkPolyData();

	vtkSmartPointer<vtkTransformPolyDataFilter> filter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();

	vtkSmartPointer<vtkLinearTransform> transform = _surface->GetTimeGeometry()->GetGeometryForTimeStep(0)->GetVtkTransform();

	transform->Update();

	filter->SetInputData(polyData);
	filter->SetTransform(transform);
	filter->Update();

	mitk::Surface::Pointer result = mitk::Surface::New();
	result->SetVtkPolyData(filter->GetOutput());
	result->CalculateBoundingBox();

	return result;
}


void SegmentationView::createMeshFromBinaryImage(Image3DType::Pointer _mask, std::string imagename, int numvertices)
{
	std::cout << " creating mesh from binary image : " << imagename << " with " << numvertices << " landmarks ... " << std::endl;

	// convert to binary image 
	typedef itk::BinaryThresholdImageFilter< Image3DType, Image3DType > binarythresholdimagefiltertype;
	binarythresholdimagefiltertype::Pointer binaryfilter = binarythresholdimagefiltertype::New();

	binaryfilter->SetInput(_mask);
	binaryfilter->SetLowerThreshold(1);
	binaryfilter->SetUpperThreshold(500);
	binaryfilter->SetInsideValue(1);
	binaryfilter->SetOutsideValue(0);
	binaryfilter->Update();

	_mask = binaryfilter->GetOutput();

	mitk::Image::Pointer mitkimage = mitk::Image::New();
	mitk::CastToMitkImage(_mask, mitkimage);

	mitk::ManualSegmentationToSurfaceFilter::Pointer surfacefilter = mitk::ManualSegmentationToSurfaceFilter::New();
	surfacefilter->SetInput(mitkimage);

	/*surfacefilter->SetMedianKernelSize(3, 3, 3); 
	surfacefilter->SetTargetReduction(1.0); 
	surfacefilter->SetMedianFilter3D(true); 
	surfacefilter->SetSmooth(true); 
	surfacefilter->UpdateLargestPossibleRegion(); */

	surfacefilter->MedianFilter3DOn(); 
	surfacefilter->SetGaussianStandardDeviation(1.5); 
	surfacefilter->InterpolationOn(); 
	surfacefilter->UseGaussianImageSmoothOn(); 
	surfacefilter->SetThreshold(0.5); 
	surfacefilter->SetDecimate(mitk::ImageToSurfaceFilter::DecimatePro); 
	surfacefilter->SetTargetReduction(0.05f); 
	surfacefilter->SmoothOn(); 
	surfacefilter->Update(); 


	vtkSmartPointer<vtkPolyData> _poly = surfacefilter->GetOutput()->GetVtkPolyData();

	vtkPolyData* remeshedpoly = vtkPolyData::New();
	remeshedpoly->DeepCopy(_poly);
	remeshedpoly = this->remeshSurface(_poly, numvertices);


	vtkSmartPointer< vtkPolyDataConnectivityFilter > connectivityfilter = vtkSmartPointer< vtkPolyDataConnectivityFilter >::New();
	connectivityfilter->SetInputData(remeshedpoly); 
	connectivityfilter->SetExtractionModeToLargestRegion();
	connectivityfilter->Update();


	// write output 
	vtkPolyDataWriter* writer = vtkPolyDataWriter::New();
	writer->SetInputData(connectivityfilter->GetOutput());
	std::string writename = imagename + ".vtk";
	writer->SetFileName(writename.c_str());
	writer->Update();

}


vtkPolyData* SegmentationView::remeshSurface(vtkPolyData* surfacePoly, int numVertices)
{
	// set default parameters 
	double gradation = 1.0;
	int subsampling = 10;
	double edgeSplitting = 0.0;
	int optimizationLevel = 1;
	bool forceManifold = false;
	bool boundaryFixing = false;

	if (surfacePoly == nullptr)
		throw std::exception("PolyData of input surface at time step 1 is NULL!");

	if (surfacePoly->GetNumberOfPolys() == 0)
		throw std::exception("Input surface has no polygons at time step 1!");

	vtkSmartPointer<vtkPolyData> surfacePolyData = vtkSmartPointer<vtkPolyData>::New();
	surfacePolyData->DeepCopy(surfacePoly);

	vtkSmartPointer<vtkSurface> mesh = vtkSmartPointer<vtkSurface>::New();

	mesh->CreateFromPolyData(surfacePolyData);
	mesh->GetCellData()->Initialize();
	mesh->GetPointData()->Initialize();

	mesh->DisplayMeshProperties();

	if (numVertices == 0)
		numVertices = surfacePolyData->GetNumberOfPoints();

	if (edgeSplitting != 0.0)
		mesh->SplitLongEdges(edgeSplitting);

	vtkSmartPointer<vtkIsotropicDiscreteRemeshing> remesher = vtkSmartPointer<vtkIsotropicDiscreteRemeshing>::New();

	remesher->GetMetric()->SetGradation(gradation);
	remesher->SetBoundaryFixing(boundaryFixing);
	remesher->SetConsoleOutput(1);
	remesher->SetForceManifold(forceManifold);
	remesher->SetInput(mesh);
	remesher->SetNumberOfClusters(numVertices);
	//remesher->SetNumberOfThreads(vtkMultiThreader::GetGlobalDefaultNumberOfThreads());
	remesher->SetSubsamplingThreshold(subsampling);

	remesher->Remesh();

	// Optimization: Minimize distance between input surface and remeshed surface
	if (optimizationLevel != 0)
	{
		ClustersQuadrics clustersQuadrics(numVertices);

		vtkSmartPointer<vtkIdList> faceList = vtkSmartPointer<vtkIdList>::New();
		vtkSmartPointer<vtkIntArray> clustering = remesher->GetClustering();
		vtkSmartPointer<vtkSurface> remesherInput = remesher->GetInput();
		int clusteringType = remesher->GetClusteringType();
		int numItems = remesher->GetNumberOfItems();
		int numMisclassifiedItems = 0;

		for (int i = 0; i < numItems; ++i)
		{
			int cluster = clustering->GetValue(i);

			if (cluster >= 0 && cluster < numVertices)
			{
				if (clusteringType != 0)
				{
					remesherInput->GetVertexNeighbourFaces(i, faceList);
					int numIds = static_cast<int>(faceList->GetNumberOfIds());

					for (int j = 0; j < numIds; ++j)
						vtkQuadricTools::AddTriangleQuadric(clustersQuadrics.Elements[cluster], remesherInput, faceList->GetId(j), false);
				}
				else
				{
					vtkQuadricTools::AddTriangleQuadric(clustersQuadrics.Elements[cluster], remesherInput, i, false);
				}
			}
			else
			{
				++numMisclassifiedItems;
			}
		}

		if (numMisclassifiedItems != 0)
			std::cout << numMisclassifiedItems << " items with wrong cluster association" << std::endl;

		vtkSmartPointer<vtkSurface> remesherOutput = remesher->GetOutput();
		double point[3];

		for (int i = 0; i < numVertices; ++i)
		{
			remesherOutput->GetPoint(i, point);
			vtkQuadricTools::ComputeRepresentativePoint(clustersQuadrics.Elements[i], point, optimizationLevel);
			remesherOutput->SetPointCoordinates(i, point);
		}

		std::cout << "After quadrics post-processing:" << std::endl;
		remesherOutput->DisplayMeshProperties();
	}

	vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();

	normals->SetInputData(remesher->GetOutput());
	normals->AutoOrientNormalsOn();
	normals->ComputeCellNormalsOff();
	normals->ComputePointNormalsOn();
	normals->ConsistencyOff();
	normals->FlipNormalsOff();
	normals->NonManifoldTraversalOff();
	normals->SplittingOff();

	normals->Update();

	//5. ensure that all filters will be destroyed
	vtkSmartPointer<vtkPolyData> polyData = normals->GetOutput();
	polyData->Register(nullptr);
	//polyData->SetSource(nullptr);  

	return polyData;

}


void SegmentationView::addImageToNode(AtlasImageType::Pointer _img, std::string imgName, bool visable)
{
	mitk::Image::Pointer outputMITK = mitk::Image::New();
	mitk::CastToMitkImage(_img, outputMITK);

	mitk::DataNode::Pointer outputNode = mitk::DataNode::New();
	outputNode->SetData(outputMITK);
	outputNode->SetName(imgName);
	outputNode->SetVisibility(visable);
	GetDataStorage()->Add(outputNode);
}


void SegmentationView::addCharImageToNode(Image3DType::Pointer _img, std::string imgName, bool visable)
{
	mitk::Image::Pointer outputMITK = mitk::Image::New();
	mitk::CastToMitkImage(_img, outputMITK);

	mitk::DataNode::Pointer outputNode = mitk::DataNode::New();
	outputNode->SetData(outputMITK);
	outputNode->SetName(imgName);
	outputNode->SetVisibility(visable);
	GetDataStorage()->Add(outputNode);
}


void SegmentationView::addPolyToNode(vtkPolyData* _poly, std::string meshName, bool visable)
{
	mitk::Surface::Pointer surface = mitk::Surface::New();
	surface->SetVtkPolyData(_poly);

	mitk::DataNode::Pointer node = mitk::DataNode::New();
	node->SetData(surface);
	node->SetName(meshName);
	node->SetVisibility(visable);
	GetDataStorage()->Add(node);
}


void SegmentationView::CenterOfMassToOrigin(std::vector< vtkSmartPointer< vtkPolyData >> &m_meshes, bool weighting)
{
	cout << " [ " << __FUNCTION__ << " ] Center meshes to Origin .. " << std::endl;

	for (unsigned int ns = 0; ns < m_meshes.size(); ns++)
	{
		vtkSmartPointer< vtkCenterOfMass > centerOfMassFilter = vtkSmartPointer< vtkCenterOfMass >::New();

		double center[3];

		if (weighting == false)
		{
			centerOfMassFilter->SetInputData(m_meshes.at(ns));
			centerOfMassFilter->SetUseScalarsAsWeights(false);
			centerOfMassFilter->Update();
			centerOfMassFilter->GetCenter(center);
		}
		else
		{
			centerOfMassFilter->ComputeCenterOfMass(m_meshes.at(ns)->GetPoints(), m_meshes.at(ns)->GetPointData()->GetArray("Weighting"), center);
		}

		// translate
		for (unsigned int np = 0; np < m_meshes.at(ns)->GetNumberOfPoints(); np++)
		{
			double point_old[3];
			m_meshes.at(ns)->GetPoint(np, point_old);

			double point_new[3];
			point_new[0] = point_old[0] - center[0];
			point_new[1] = point_old[1] - center[1];
			point_new[2] = point_old[2] - center[2];
			m_meshes.at(ns)->GetPoints()->SetPoint(np, point_new);
		}
	}
}


void SegmentationView::visualizeDatasets(std::vector< vtkSmartPointer< vtkPolyData > > &m_datasets, std::string name0, bool visible)
{
	std::cout << " [ " << __FUNCTION__ << " ] Visualize " << m_datasets.size() << " datasets starting with " << name0 << std::endl;

	for (unsigned int i = 0; i < m_datasets.size(); i++)
	{
		vtkSmartPointer< vtkPolyData > poly = m_datasets.at(i);
		mitk::DataNode::Pointer node = mitk::DataNode::New();
		mitk::Surface::Pointer surface = mitk::Surface::New();
		surface->SetVtkPolyData(poly);
		node->SetData(surface);

		std::stringstream ss;
		ss << i + 1;
		std::string countNum = ss.str();

		if (i + 1 < 10) countNum = "00" + countNum;
		else countNum = "0" + countNum;

		std::string fileName = name0 + countNum;
		node->SetName(fileName.c_str());
		node->SetVisibility(visible);
		GetDataStorage()->Add(node);
	}

}


void SegmentationView::writeAtlasImage(AtlasImageType::Pointer img, std::string writeName)
{
	WriterFilterType::Pointer imgWriter = WriterFilterType::New();
	imgWriter->SetInput(img);
	imgWriter->SetFileName(writeName.c_str());
	imgWriter->Update();
}


void SegmentationView::writeChar3DImage(Image3DType::Pointer img, std::string writeName)
{
	typedef itk::ImageFileWriter< Image3DType > WriterType; 
	WriterType::Pointer writer = WriterType::New(); 
	writer->SetInput(img); 
	writer->SetFileName(writeName.c_str()); 
	writer->Update(); 
}


void SegmentationView::WriteOutPolyData(vtkSmartPointer< vtkPolyData > poly, std::string fileName)
{
	vtkSmartPointer< vtkPolyDataWriter > writer = vtkSmartPointer< vtkPolyDataWriter >::New();

	writer->SetInputData(poly);
	writer->SetFileName(fileName.c_str());
	writer->Update();
}


void SegmentationView::ResetPolyData(vtkSmartPointer< vtkPolyData > poly)
{
	poly = NULL;
	poly = vtkSmartPointer< vtkPolyData >::New();
}


