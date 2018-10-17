
// local includes
#include "ImageAnalyser.h"


ImageAnalyser::ImageAnalyser(mitk::DataNode::Pointer imageNode, AtlasImageType::Pointer probMapImage)
{	
	m_probMapImage = AtlasImageType::New(); 
	ImageDuplicatorType::Pointer duplicator = ImageDuplicatorType::New(); 
	duplicator->SetInputImage(probMapImage); 
	duplicator->Update(); 
	m_probMapImage = duplicator->GetOutput(); 

	m_itkImage = AtlasImageType::New();
	if (!imageNode) return;
	m_imageNode = imageNode;
	m_mitkImage = dynamic_cast<mitk::Image*>(imageNode->GetData());
	if(!m_mitkImage) return;
	mitk::CastToItkImage(m_mitkImage, m_itkImage);

	m_ROIImage = AtlasImageType::New(); 
	m_NOIImage = AtlasImageType::New(); 

	RescaleIntensityProbMap();  // get the largest component by the way 

	this->PriorProb_ROI = 0.0; 
	this->PriorProb_NOI = 0.0; 

	this->PriorProb_ROI_STD = 0.0; 
	this->PriorProb_NOI_STD = 0.0; 

	this->ROI_MEAN = 90;  // by default 
	this->ROI_STD = 30; 
	this->NOI_MEAN = 30; 
	this->NOI_STD = 100;  

}


ImageAnalyser::~ImageAnalyser()
{

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void ImageAnalyser::landmarkDeformation(vtkPolyData* _inputPoly, vtkPolyData* outputPoly, double _lambda, int _radius, int & _onBoundary)
{
	if (!_inputPoly->GetPointData()->GetArray("weights"))
	{
		MITK_WARN << __FUNCTION__ << " no weights property found ... ";
		return;
	}

	vtkSmartPointer< vtkDataArray > weights = vtkSmartPointer< vtkDoubleArray >::New();
	weights = _inputPoly->GetPointData()->GetArray("weights");

	vtkPoints* outputPoints = outputPoly->GetPoints();
	vtkDataArray* normals = _inputPoly->GetPointData()->GetNormals();

	_onBoundary = 0;

	for (long ptId = 0; ptId < _inputPoly->GetNumberOfPoints(); ptId++)
	{
		itk::Point< double, 3 > basePoint;
		basePoint[0] = _inputPoly->GetPoint(ptId)[0];
		basePoint[1] = _inputPoly->GetPoint(ptId)[1];
		basePoint[2] = _inputPoly->GetPoint(ptId)[2];

		itk::Point< double, 3 > baseNormal;
		baseNormal[0] = normals->GetTuple(ptId)[0];
		baseNormal[1] = normals->GetTuple(ptId)[1];
		baseNormal[2] = normals->GetTuple(ptId)[2];

		itk::Point< double, 3 > result;
		result[0] = basePoint[0];
		result[1] = basePoint[1];
		result[2] = basePoint[2];

		itk::Point< double, 3 > interiorPoint; 
		interiorPoint[0] = basePoint[0] - baseNormal[0]; 
		interiorPoint[1] = basePoint[1] - baseNormal[1]; 
		interiorPoint[2] = basePoint[2] - baseNormal[2]; 

		itk::Point< double, 3 > outeriorPoint; 
		outeriorPoint[0] = basePoint[0] + baseNormal[0]; 
		outeriorPoint[1] = basePoint[1] + baseNormal[1]; 
		outeriorPoint[2] = basePoint[2] + baseNormal[2]; 

		AtlasImageType::IndexType baseIndex;
		if (!this->m_itkImage->TransformPhysicalPointToIndex(basePoint, baseIndex))
		{
			MITK_WARN << __FUNCTION__ << " exceed boundary " << endl;

			outputPoints->SetPoint(ptId, basePoint[0], basePoint[1], basePoint[2]);

			continue;  // not consider this point with current radius 
		}

		bool inside = false;
		bool outside = false;

		double InROI = 0.0;
		double OutNOI = 0.0;

		this->calculateProbInside(baseIndex, _lambda, _radius, InROI);
		this->calculateProbOutside(baseIndex, _lambda, _radius, OutNOI);

		// >>> check whether inside organ satisfies 
		if (this->checkBoundary(basePoint, _lambda, _radius))
		{
//			cout << " ptId = " << ptId << " on boundary : inROI = " << InROI << ", outNOI = " << OutNOI << endl;

			_onBoundary++;

			outputPoints->SetPoint(ptId, result[0], result[1], result[2]);

			weights->SetTuple1(ptId, 1.0);

			continue;
		}

		// Inside not in range but Outside in in range
		if (InROI < this->PriorProb_ROI - this->PriorProb_ROI_STD && OutNOI >= this->PriorProb_NOI - this->PriorProb_NOI_STD)
		{
			// move inside
//			cout << " ptId = " << ptId << " inside false (move inside) : inROI = " << InROI << ", outNOI = " << OutNOI << endl;

			result[0] = interiorPoint[0];
			result[1] = interiorPoint[1];
			result[2] = interiorPoint[2];

			outputPoints->SetPoint(ptId, result[0], result[1], result[2]);

			weights->SetTuple1(ptId, 0.8);
		}

		// Outside not in range but Inside is in range
		else if (InROI >= this->PriorProb_ROI - this->PriorProb_ROI_STD && OutNOI < this->PriorProb_NOI - this->PriorProb_NOI_STD)
		{
//			cout << " ptId = " << ptId << " outside false (move outside) : inROI = " << InROI << ", outNOI = " << OutNOI << endl;

			result[0] = outeriorPoint[0];
			result[1] = outeriorPoint[1];
			result[2] = outeriorPoint[2];

			outputPoints->SetPoint(ptId, result[0], result[1], result[2]);

			weights->SetTuple1(ptId, 0.8);
		}

		// Find the next Point 
		else
		{		
			AtlasImageType::IndexType interiorIndex; 
			if (!m_itkImage->TransformPhysicalPointToIndex(interiorPoint, interiorIndex))
			{
//				cout << " ptId = " << ptId << " exceed (move outside) : inROI = " << InROI << ", outNOI = " << OutNOI << endl;

				result[0] = outeriorPoint[0]; 
				result[1] = outeriorPoint[1]; 
				result[2] = outeriorPoint[2]; 
				
				outputPoints->SetPoint(ptId, result[0], result[1], result[2]);

				weights->SetTuple1(ptId, 0.6);

				continue; 
			}

			AtlasImageType::IndexType outeriorIndex; 
			if (!m_itkImage->TransformPhysicalPointToIndex(outeriorPoint, outeriorIndex))
			{
//				cout << " ptId = " << ptId << " exceed (move inside) : inROI = " << InROI << ", outNOI = " << OutNOI << endl;

				result[0] = interiorPoint[0]; 
				result[1] = interiorPoint[1];
				result[2] = interiorPoint[2]; 
				
				outputPoints->SetPoint(ptId, result[0], result[1], result[2]);

				weights->SetTuple1(ptId, 0.6);

				continue;
			}

			double interiorROI = 0.0;
			double interiorNOI = 0.0;

			double outeriorROI = 0.0;
			double outeriorNOI = 0.0;

			this->calculateProbOfNeighborhood(interiorIndex, _lambda, 2, interiorROI, interiorNOI); 
			this->calculateProbOfNeighborhood(outeriorIndex, _lambda, 2, outeriorROI, outeriorNOI); 

			if (interiorROI >= InROI && interiorROI > outeriorROI)
			{
				// move inside 
//				cout << " ptId = " << ptId << " interiorROI > outeriorROI (move inside) : inROI = " << InROI << ", outNOI = " << OutNOI << endl;

				result[0] = (basePoint[0] + interiorPoint[0]) / 2; 
				result[1] = (basePoint[1] + interiorPoint[1]) / 2; 
				result[2] = (basePoint[2] + interiorPoint[2]) / 2; 
				
				outputPoints->SetPoint(ptId, result[0], result[1], result[2]);

				weights->SetTuple1(ptId, 0.2);
			}
			else if (outeriorROI >= interiorROI && outeriorROI > InROI)
			{
				// move outside 
//				cout << " ptId = " << ptId << " outeriorNOI > interiorNOI (move outside) : inROI = " << InROI << ", outNOI = " << OutNOI << endl;

				result[0] = (basePoint[0] + outeriorPoint[0]) / 2; 
				result[1] = (basePoint[1] + outeriorPoint[1]) / 2;
				result[2] = (basePoint[2] + outeriorPoint[2]) / 2;
				
				outputPoints->SetPoint(ptId, result[0], result[1], result[2]);

				weights->SetTuple1(ptId, 0.2);
			}

			else
			{
				outputPoints->SetPoint(ptId, basePoint[0], basePoint[1], basePoint[2]);

				weights->SetTuple1(ptId, 0.0);
			}
			
		}

		if (this->checkBoundary(result, _lambda, _radius))
		{
			weights->SetTuple1(ptId, 1.0);

			_onBoundary++;
		}

	}

	weights->Modified();
	outputPoly->GetPointData()->SetScalars(weights);
	outputPoints->Modified();
	outputPoly->Modified();

	cout << __FUNCTION__ << " boundary number = " << _onBoundary << " .. end deformation ... " << endl;

//	this->trainPriorGMM(outputPoly);    // update GMM 

//	this->computePriorProbs(outputPoly, 0.0, _radius);   // update m_imgROI && m_imgNOI 

}


void ImageAnalyser::calculateProbInside(AtlasImageType::IndexType _baseIndex, double _lambdaNN, int _radius, double& aveProbROI)
{
	if (!m_ROIImage) return; 

	// extract region with radius 
	AtlasImageType::IndexType startIndex;
	startIndex[0] = _baseIndex[0] - _radius;
	startIndex[1] = _baseIndex[1] - _radius;
	startIndex[2] = _baseIndex[2] - _radius;

	AtlasImageType::SizeType size;
	size[0] = 2 * _radius + 1;
	size[1] = 2 * _radius + 1;
	size[2] = 2 * _radius + 1;

	for (int i = 0; i < 3; i++)  // important ! 
	{
		if (startIndex[i] < m_ROIImage->GetLargestPossibleRegion().GetIndex()[i] || startIndex[i] + size[i] > m_ROIImage->GetLargestPossibleRegion().GetSize()[i])
		{
			MITK_WARN << __FUNCTION__ << " image boundary exceed !";
			return;
		}
	}

	AtlasImageType::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(startIndex);

	itk::ImageRegionConstIterator< AtlasImageType > it(m_ROIImage, desiredRegion);

	int TotalPixels = 0; 
	aveProbROI = 0.0; 
	
	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		InternalPixelType gValue = it.Get();

		if (gValue < -1024) continue; 

		AtlasImageType::IndexType index = it.GetIndex();

		double currProbROI = this->getGaussianProbOfROI(gValue); 	// double currProbROI = this->getProbOfROI(gValue); 
		double currProbInNOI = this->getGaussianProbOfNOI(gValue); 
		double currProbNN = this->getPriorProbBasedOnNNMap(index);

		//>>> weighted GMM <<<// 
		
		if (_lambdaNN == 0) currProbNN = 0.5; 

		double ProbTotal = currProbROI * currProbNN + currProbInNOI * (1 - currProbNN);

		if (ProbTotal == 0) continue; 

		aveProbROI += currProbROI * currProbNN / ProbTotal;

		//>>>>>>>>>> <<<<<<<<<<// 

		/*TotalProbROI += currProbROI + _lambdaNN * currProbNN;
		TotalInNOI += currProbInNOI + _lambdaNN * (1 - currProbNN);*/

		TotalPixels++; 
	}

	if (TotalPixels == 0) aveProbROI = 0.0;
	
	else aveProbROI /= TotalPixels;

}


void ImageAnalyser::calculateProbOutside(AtlasImageType::IndexType _baseIndex, double _lambdaNN, int _radius, double& aveProbNOI)
{
	if (!m_NOIImage) return; 

	// extract region with radius 
	AtlasImageType::IndexType startIndex;
	startIndex[0] = _baseIndex[0] - _radius;
	startIndex[1] = _baseIndex[1] - _radius;
	startIndex[2] = _baseIndex[2] - _radius;

	AtlasImageType::SizeType size;
	size[0] = 2 * _radius + 1;
	size[1] = 2 * _radius + 1;
	size[2] = 2 * _radius + 1;

	for (int i = 0; i < 3; i++)  // important ! 
	{
		if (startIndex[i] < m_NOIImage->GetLargestPossibleRegion().GetIndex()[i] || startIndex[i] + size[i] > m_NOIImage->GetLargestPossibleRegion().GetSize()[i])
		{
			MITK_WARN << __FUNCTION__ << " image boundary exceed !";
			return;
		}
	}

	AtlasImageType::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(startIndex);

	itk::ImageRegionConstIterator< AtlasImageType > it(m_NOIImage, desiredRegion);

	int TotalPixels = 0;
	aveProbNOI = 0.0; 

	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		InternalPixelType gValue = it.Get();

		if (gValue < -1024) continue; 

		AtlasImageType::IndexType index = it.GetIndex();

		double currOutROI = this->getGaussianProbOfROI(gValue);  
		double currProbNOI = this->getGaussianProbOfNOI(gValue);  
		double currProbNN = this->getPriorProbBasedOnNNMap(index);

		//>>>>> Weighted GMM <<<<<//

		if (_lambdaNN == 0) currProbNN = 0.5;

		double ProbTotal = currOutROI * currProbNN + currProbNOI * (1 - currProbNN);

		if (ProbTotal == 0) continue;

		aveProbNOI += currProbNOI * (1 - currProbNN) / ProbTotal;

		//>>>>>>>>>>> <<<<<<<<<<<//
		
		/*TotalProbNOI += currProbNOI + _lambdaNN * (1 - currProbNN);
		TotalOutROI += currOutROI + _lambdaNN * currProbNN;*/

		TotalPixels++;
	}

	if (TotalPixels == 0) aveProbNOI = 0.0;
	
	else aveProbNOI /= TotalPixels; 

}


void ImageAnalyser::calculateProbOfNeighborhood(AtlasImageType::IndexType _baseIndex, double _lambdaNN, int _radius, double & ProbROI, double & ProbNOI)
{
	ProbROI = 0.0; 
	ProbNOI = 0.0; 

	// extract region with radius 
	AtlasImageType::IndexType startIndex;
	startIndex[0] = _baseIndex[0] - _radius;
	startIndex[1] = _baseIndex[1] - _radius;
	startIndex[2] = _baseIndex[2] - _radius;

	AtlasImageType::SizeType size;
	size[0] = 2 * _radius + 1;
	size[1] = 2 * _radius + 1;
	size[2] = 2 * _radius + 1;

	for (int i = 0; i < 3; i++)  // important ! 
	{
		if (startIndex[i] < m_itkImage->GetLargestPossibleRegion().GetIndex()[i] || startIndex[i] + size[i] > m_itkImage->GetLargestPossibleRegion().GetSize()[i])
		{
			MITK_WARN << __FUNCTION__ << " image boundary exceed !";
			return;
		}
	}

	AtlasImageType::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(startIndex);

	itk::ImageRegionConstIterator< AtlasImageType > it(m_itkImage, desiredRegion);

	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		InternalPixelType gValue = it.Get();

		AtlasImageType::IndexType index = it.GetIndex();

		double currProbROI = this->getGaussianProbOfROI(gValue); 	
		double currProbInNOI = this->getGaussianProbOfNOI(gValue);
		double currProbNN = this->getPriorProbBasedOnNNMap(index);

		//>>> weighted GMM <<<// 

		if (_lambdaNN == 0) currProbNN = 0.5;

		double ProbTotal = currProbROI * currProbNN + currProbInNOI * (1 - currProbNN);

		ProbROI += currProbROI * currProbNN / ProbTotal; 
		ProbNOI += currProbInNOI * (1 - currProbNN) / ProbTotal; 
	}

	ProbROI /= size[0] * size[1] * size[2]; 
	ProbNOI /= size[0] * size[1] * size[2]; 

}


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Pre Training >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


void ImageAnalyser::trainProbsROIandNOI(vtkPolyData* _inputPoly, double _lambda, int _radius)
{
	mitk::Surface::Pointer surface = mitk::Surface::New();
	surface->SetVtkPolyData(_inputPoly);
	surface->Modified();

	mitk::SurfaceToImageFilter::Pointer surfaceToImageFilter = mitk::SurfaceToImageFilter::New();
	surfaceToImageFilter->SetInput(surface);
	surfaceToImageFilter->SetImage(m_mitkImage);  
	surfaceToImageFilter->Update();

	mitk::ImageTimeSelector::Pointer timeSelector = mitk::ImageTimeSelector::New();
	timeSelector->SetInput(surfaceToImageFilter->GetOutput());
	timeSelector->SetTimeNr(0);
	timeSelector->UpdateLargestPossibleRegion();
	timeSelector->Update();

	/** get m_ROIImage **/
	mitk::CastToItkImage(timeSelector->GetOutput(), m_ROIImage);

	/** get m_NOIImage **/
	ImageDuplicatorType::Pointer duplicator = ImageDuplicatorType::New();
	duplicator->SetInputImage(m_ROIImage);
	duplicator->Update();
	m_NOIImage = duplicator->GetOutput();

	itk::ImageRegionIterator< AtlasImageType > imgITK(m_itkImage, m_itkImage->GetRequestedRegion());
	itk::ImageRegionIterator< AtlasImageType > imgNOI(m_NOIImage, m_NOIImage->GetRequestedRegion());
	for (imgITK.GoToBegin(), imgNOI.GoToBegin(); !imgITK.IsAtEnd() && !imgNOI.IsAtEnd(); ++imgITK, ++imgNOI)
	{
		AtlasImageType::IndexType currIndex = imgNOI.GetIndex();

		if (imgNOI.Get() > -10000) imgNOI.Set(-10000);
		else imgNOI.Set(imgITK.Get());
	}

	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	vtkPoints* inputPoints = _inputPoly->GetPoints();
	vtkDataArray* normals = _inputPoly->GetPointData()->GetNormals();

	vnl_vector< double > InROIs(_inputPoly->GetNumberOfPoints(), 0.0);
	vnl_vector< double > OutNOIs(_inputPoly->GetNumberOfPoints(), 0.0); 
		
	for (long ptId = 0; ptId < _inputPoly->GetNumberOfPoints(); ptId++)
	{
		itk::Point< double, 3 > basePoint;
		basePoint[0] = inputPoints->GetPoint(ptId)[0];
		basePoint[1] = inputPoints->GetPoint(ptId)[1];
		basePoint[2] = inputPoints->GetPoint(ptId)[2];

		itk::Point< double, 3 > baseNormal;
		baseNormal[0] = normals->GetTuple(ptId)[0];
		baseNormal[1] = normals->GetTuple(ptId)[1];
		baseNormal[2] = normals->GetTuple(ptId)[2];

		itk::Point< double, 3 > result;
		result[0] = basePoint[0]; result[1] = basePoint[1]; result[2]; basePoint[2];

		AtlasImageType::IndexType baseIndex;
		if (!this->m_itkImage->TransformPhysicalPointToIndex(basePoint, baseIndex))
		{
			MITK_WARN << __FUNCTION__ << " exceed boundary " << endl;
			InROIs[ptId] = 0.5; 
			OutNOIs[ptId] = 0.5;   
			continue;
		}

		this->calculateProbInside(baseIndex, _lambda, _radius, InROIs[ptId]);
		this->calculateProbOutside(baseIndex, _lambda, _radius, OutNOIs[ptId]);

		cout << " ptID : " << ptId << " InROI = " << InROIs[ptId] << " , OutNOIs = " << OutNOIs[ptId] << endl; 

	}

	this->PriorProb_ROI = InROIs.mean(); 
	this->PriorProb_NOI = OutNOIs.mean(); 

	double sumInROI = 0.0; 
	double sumOutNOI = 0.0;
	for (int i = 0; i < _inputPoly->GetNumberOfPoints(); i++)
	{
		sumInROI += pow(InROIs[i] - this->PriorProb_ROI, 2);
		sumOutNOI += pow(OutNOIs[i] - this->PriorProb_NOI, 2);
	}

	this->PriorProb_ROI_STD = sqrt(sumInROI / _inputPoly->GetNumberOfPoints()); 
	this->PriorProb_NOI_STD = sqrt(sumOutNOI / _inputPoly->GetNumberOfPoints());   
	
	cout << __FUNCTION__ << " Prior_InROI(std) = " << std::setprecision(4) << this->PriorProb_ROI << " +- " << this->PriorProb_ROI_STD <<
		" , Prior_OutNOI(std) = " << this->PriorProb_NOI << " +- " << this->PriorProb_NOI_STD << endl;
	cout << endl; 
}


void ImageAnalyser::computePriorProbs(vtkPolyData* _inputPoly, double _lambda, int _radius)
{
	cout << __FUNCTION__ << " Updating prior probs : ";

	mitk::Surface::Pointer surface = mitk::Surface::New();
	surface->SetVtkPolyData(_inputPoly);
	surface->Modified();

	mitk::SurfaceToImageFilter::Pointer surfaceToImageFilter = mitk::SurfaceToImageFilter::New();
	surfaceToImageFilter->SetInput(surface);
	surfaceToImageFilter->SetImage(m_mitkImage);
	surfaceToImageFilter->Update();

	mitk::ImageTimeSelector::Pointer timeSelector = mitk::ImageTimeSelector::New();
	timeSelector->SetInput(surfaceToImageFilter->GetOutput());
	timeSelector->SetTimeNr(0);
	timeSelector->UpdateLargestPossibleRegion();
	timeSelector->Update();

	/** get m_ROIImage **/
	mitk::CastToItkImage(timeSelector->GetOutput(), m_ROIImage);

	std::vector< double > ProbROIs;
	this->PriorProb_ROI = 0.0;

	itk::ImageRegionIterator< AtlasImageType > imgROI(m_ROIImage, m_ROIImage->GetRequestedRegion());
	for (imgROI.GoToBegin(); !imgROI.IsAtEnd(); ++imgROI)
	{
		if (imgROI.Get() < -1024) continue;

		InternalPixelType gValue = imgROI.Get();

		double currInROI = this->getGaussianProbOfROI(gValue);
		double currInNOI = this->getGaussianProbOfNOI(gValue);
		double currProbNN = this->getPriorProbBasedOnNNMap(imgROI.GetIndex());

		if (_lambda == 0) currProbNN = 0.5;

		double probTotal = currInROI * currProbNN + currInNOI * (1 - currProbNN);

		if (probTotal == 0) continue;

		ProbROIs.push_back(currInROI * currProbNN / probTotal);
		this->PriorProb_ROI += currInROI * currProbNN / probTotal;
	}

	if (ProbROIs.size() == 0)
	{
		this->PriorProb_ROI = 0.0;
		this->PriorProb_ROI_STD = 0.0;
	}
	else
	{
		this->PriorProb_ROI /= ProbROIs.size();

		double sumROI = 0.0;
		for (int i = 0; i < ProbROIs.size(); i++)
			sumROI += pow(ProbROIs.at(i) - this->PriorProb_ROI, 2);

		this->PriorProb_ROI_STD = sumROI / ProbROIs.size();
		this->PriorProb_ROI_STD = sqrt(this->PriorProb_ROI_STD);
	}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	/** get m_NOIImage **/
	ImageDuplicatorType::Pointer duplicator = ImageDuplicatorType::New();
	duplicator->SetInputImage(m_ROIImage);
	duplicator->Update();
	m_NOIImage = duplicator->GetOutput();

	itk::ImageRegionIterator< AtlasImageType > imgITK(m_itkImage, m_itkImage->GetRequestedRegion());
	itk::ImageRegionIterator< AtlasImageType > img_NOI(m_NOIImage, m_NOIImage->GetRequestedRegion());
	for (imgITK.GoToBegin(), img_NOI.GoToBegin(); !imgITK.IsAtEnd() && !img_NOI.IsAtEnd(); ++imgITK, ++img_NOI)
	{
		AtlasImageType::IndexType currIndex = img_NOI.GetIndex();

		if (img_NOI.Get() > -10000) img_NOI.Set(-10000);
		else img_NOI.Set(imgITK.Get());
	}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	vtkDataArray* normals = _inputPoly->GetPointData()->GetNormals();

	vnl_vector< double > ProbNOIs(_inputPoly->GetNumberOfPoints(), 0.0); 

	int radiusBackup = _radius; 

	for (int ptId = 0; ptId < _inputPoly->GetNumberOfPoints(); ptId++)
	{
		itk::Point< double, 3 > basePoint;
		basePoint[0] = _inputPoly->GetPoint(ptId)[0];
		basePoint[1] = _inputPoly->GetPoint(ptId)[1];
		basePoint[2] = _inputPoly->GetPoint(ptId)[2];

		AtlasImageType::IndexType baseIndex;
		if (!this->m_itkImage->TransformPhysicalPointToIndex(basePoint, baseIndex))
		{
			MITK_WARN << __FUNCTION__ << " exceed boundary outside " << endl;
			ProbNOIs[ptId] = 0.5; 
			continue;
		}

		this->calculateProbOutside(baseIndex, _lambda, _radius, ProbNOIs[ptId]);

	}

	this->PriorProb_NOI = ProbNOIs.mean(); 

	double sumNOI = 0.0; 
	for (int i = 0; i < ProbNOIs.size(); i++)
		sumNOI += pow(ProbNOIs[i] - this->PriorProb_NOI, 2);

	this->PriorProb_NOI_STD = sumNOI / ProbNOIs.size(); 
	this->PriorProb_NOI_STD = sqrt(this->PriorProb_NOI_STD); 

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	cout << " InROI (std) = " << this->PriorProb_ROI << " ( " << this->PriorProb_ROI_STD << " ) , OutNOI (std) = " << this->PriorProb_NOI << " ( " << this->PriorProb_NOI_STD << " ) " << endl; 
	
}


void ImageAnalyser::trainPriorGMM(vtkPolyData* _inputPoly)
{
	cout << __FUNCTION__ << " with input surface : "; 

	double* bounds = _inputPoly->GetBounds(); 
	itk::Point< double, 3 > leftPoint;  leftPoint[0] = bounds[0]; leftPoint[1] = bounds[2]; leftPoint[2] = bounds[4]; 
	itk::Point< double, 3 > rightPoint; rightPoint[0] = bounds[1]; rightPoint[1] = bounds[3]; rightPoint[2] = bounds[5]; 
	
	AtlasImageType::IndexType leftCorner;  this->m_itkImage->TransformPhysicalPointToIndex(leftPoint, leftCorner); 
	AtlasImageType::IndexType rightCorner; this->m_itkImage->TransformPhysicalPointToIndex(rightPoint, rightCorner); 

	leftCorner[0] -= 10; leftCorner[1] -= 10; leftCorner[2] -= 10; 
	rightCorner[0] += 10; rightCorner[1] += 10; rightCorner[2] += 10; 

	for (int i = 0; i < 3; i++)
	{
		leftCorner[i] = max(m_itkImage->GetLargestPossibleRegion().GetIndex()[i], leftCorner[i]); 
		rightCorner[i] = min((AtlasImageType::IndexValueType)m_itkImage->GetLargestPossibleRegion().GetSize()[i], rightCorner[i]);
	}

	mitk::Surface::Pointer surface = mitk::Surface::New();
	surface->SetVtkPolyData(_inputPoly);
	surface->Modified();

	mitk::SurfaceToImageFilter::Pointer surfaceToImageFilter = mitk::SurfaceToImageFilter::New();
	surfaceToImageFilter->SetInput(surface);
	surfaceToImageFilter->SetImage(m_mitkImage);
	surfaceToImageFilter->Update();

	mitk::ImageTimeSelector::Pointer timeSelector = mitk::ImageTimeSelector::New();
	timeSelector->SetInput(surfaceToImageFilter->GetOutput());
	timeSelector->SetTimeNr(0);
	timeSelector->UpdateLargestPossibleRegion();
	timeSelector->Update();

	/** get m_ROIImage **/
	mitk::CastToItkImage(timeSelector->GetOutput(), m_ROIImage);

	statisticFilterType::Pointer statFilterROI = statisticFilterType::New();
	statFilterROI->SetInput(m_ROIImage);
	statFilterROI->Update();

	this->ROI_MEAN = statFilterROI->GetMeanVoxelValue();
	this->ROI_STD = statFilterROI->GetStandardDeviation();


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


	ImageDuplicatorType::Pointer duplicator = ImageDuplicatorType::New(); 
	duplicator->SetInputImage(m_ROIImage); 
	duplicator->Update(); 
	m_NOIImage = duplicator->GetOutput(); 

	/** get m_NOIImage **/
	itk::ImageRegionIterator< AtlasImageType > imgITK(m_itkImage, m_itkImage->GetRequestedRegion()); 
	itk::ImageRegionIterator< AtlasImageType > imgNOI(m_NOIImage, m_NOIImage->GetRequestedRegion()); 
	for (imgITK.GoToBegin(), imgNOI.GoToBegin(); !imgITK.IsAtEnd() && !imgNOI.IsAtEnd(); ++imgITK, ++imgNOI)
	{
		AtlasImageType::IndexType currIndex = imgNOI.GetIndex(); 

		if (currIndex[0] >= leftCorner[0] && currIndex[0] <= rightCorner[0]
			&& currIndex[1] >= leftCorner[1] && currIndex[1] <= rightCorner[1]
			&& currIndex[2] >= leftCorner[2] && currIndex[2] <= rightCorner[2])
		{
			if (imgNOI.Get() > -10000) imgNOI.Set(-10000);
			else imgNOI.Set(imgITK.Get());
		}
		
	}

	statisticFilterType::Pointer statFilterNOI = statisticFilterType::New(); 
	statFilterNOI->SetInput(m_NOIImage); 
	statFilterNOI->Update(); 

	this->NOI_MEAN = statFilterNOI->GetMeanVoxelValue(); 
	this->NOI_STD = statFilterNOI->GetStandardDeviation(); 

	cout << " ROI_MEAN (ROI_STD) = " << this->ROI_MEAN << " ( " << this->ROI_STD << " ) NOI_MEAN (NOI_STD) = " << this->NOI_MEAN << " ( " << this->NOI_STD << " ) " << endl; 

}


bool ImageAnalyser::checkBoundary(itk::Point< double, 3 > _point, double _lambda, int _radius)
{

	AtlasImageType::IndexType baseIndex;
	if (!this->m_itkImage->TransformPhysicalPointToIndex(_point, baseIndex))
	{
		MITK_WARN << __FUNCTION__ << " exceed boundary " << endl;
		return false;
	}

	double InROI = 0.0; 
	double OutNOI = 0.0; 

	this->calculateProbInside(baseIndex, _lambda, _radius, InROI);
	this->calculateProbOutside(baseIndex, _lambda, _radius, OutNOI);

	if (InROI >= 0.5 && OutNOI >= 0.5) return true; 
	else return false; 


}


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


double ImageAnalyser::getGaussianProbOfROI(InternalPixelType _pIntens)
{
	double probROI = 0.0;

	if (_pIntens <= this->ROI_MEAN - 3 * this->ROI_STD || _pIntens >= this->ROI_MEAN + 3 * this->ROI_STD) return 0.0; 

	// Gaussian model of intensity 
	probROI = exp(- (double)std::pow(_pIntens - this->ROI_MEAN, 2) / (double)(2 * std::pow(this->ROI_STD, 2)));

	probROI /= this->ROI_STD * sqrt(2 * M_PI);

	return probROI;
}


double ImageAnalyser::getGaussianProbOfNOI(InternalPixelType _pIntens)
{
	double probNOI = 0.0;

	if (_pIntens <= this->NOI_MEAN - 3 * this->NOI_STD || _pIntens >= this->NOI_MEAN + 3 * this->NOI_STD) return 0.0; 

	probNOI = exp(-(double)std::pow(_pIntens - this->NOI_MEAN, 2) / (double)(2 * std::pow(this->NOI_STD, 2)));

	probNOI /= this->NOI_STD * sqrt(2 * M_PI);

	return probNOI;
}


double ImageAnalyser::getPriorProbBasedOnNNMap(AtlasImageType::IndexType _inputIndex)
{
	InternalPixelType pixelValue = 0.0f; 
	pixelValue = this->m_probMapImage->GetPixel(_inputIndex);

	if (pixelValue == 255) return 1.0; 
	else return 1 / (1 + exp(- 3 * pixelValue / 255));
	
}


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


void ImageAnalyser::RescaleIntensityProbMap()
{
	if (!this->m_probMapImage)
	{
		MITK_WARN << __FUNCTION__ << " please load probability map ; "; 
		return; 
	}

	itk::ImageRegionIterator< AtlasImageType > it(m_probMapImage, m_probMapImage->GetRequestedRegion()); 
	for (it.GoToBegin(); !it.IsAtEnd(); it++)
	{
		if (it.Get() >= 600) it.Set(255); 
		else if (it.Get() <= 200) it.Set(0); 
		else it.Set((it.Get() - 200) * 255 / 400 ); 
	}

	m_probMapImage->Modified(); 
	m_probMapImage->Update(); 
}
