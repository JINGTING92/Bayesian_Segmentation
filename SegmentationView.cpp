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

		m_Controls.setSigma->setEnabled(true); 
		m_Controls.setSigma->setMinimum(0.001); 
		m_Controls.setSigma->setMaximum(10000); 
		m_Controls.setSigma->setValue(100); 

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
		}
	}


/// >>>>>>>>>>>>>>>>>>>>> create m_imageAnalyser && Localization >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	mitk::DataNode::Pointer nodeImg = m_Controls.m_segImageSelect->GetSelectedNode();
	mitk::BaseData* data = nodeImg->GetData();
	mitk::Image* image = dynamic_cast<mitk::Image*>(data);

	m_imageAnalyser = new ImageAnalyser(nodeImg, this->m_probMap);

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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
