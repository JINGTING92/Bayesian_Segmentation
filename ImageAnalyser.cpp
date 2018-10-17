// ==================================================================== //
//																		//
// Filename: ImageAnalyser.cpp					        			    //
//																		//
// Author:	Fraunhofer Institute for Computer Graphics (IGD)			//
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
// Modified         : 2017-2018     Pancreas Segmentation is added
// ===================================================================	//

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

	//this->PriorProb_InROI = 0.0; 
	//this->PriorProb_InNOI = 0.0; 
	//this->PriorProb_OutROI = 0.0; 
	//this->PriorProb_OutNOI = 0.0; 

	this->PriorProb_ROI = 0.0; 
	this->PriorProb_NOI = 0.0; 

	this->PriorProb_ROI_STD = 0.0; 
	this->PriorProb_NOI_STD = 0.0; 

	this->ROI_MEAN = 90;  // by default 
	this->ROI_STD = 30; 
	this->NOI_MEAN = 30; 
	this->NOI_STD = 100;  

	////////////////////////////

//	EnhanceIntensityWithProbMap(); 

//	ComputeEdgeImage(m_itkImage);  // use m_maskImage -> edge value is too large : over 1,000 

	cout << __FUNCTION__ << " PreprocessImage Done .. " << endl;
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

//	double* bounds = _inputPoly->GetBounds();
//	itk::Point< double, 3 > leftPoint;  leftPoint[0] = bounds[0]; leftPoint[1] = bounds[2]; leftPoint[2] = bounds[4];
//	itk::Point< double, 3 > rightPoint; rightPoint[0] = bounds[1]; rightPoint[1] = bounds[3]; rightPoint[2] = bounds[5];
//
//	AtlasImageType::IndexType leftCorner;  this->m_itkImage->TransformPhysicalPointToIndex(leftPoint, leftCorner);
//	AtlasImageType::IndexType rightCorner; this->m_itkImage->TransformPhysicalPointToIndex(rightPoint, rightCorner);
//
//	leftCorner[0] -= 3; leftCorner[1] -= 3; leftCorner[2] -= 3;
//	rightCorner[0] += 3; rightCorner[1] += 3; rightCorner[2] += 3;
//
//	for (int i = 0; i < 3; i++)
//	{
//		leftCorner[i] = max(m_itkImage->GetLargestPossibleRegion().GetIndex()[i], leftCorner[i]);
//		rightCorner[i] = min((AtlasImageType::IndexValueType)m_itkImage->GetLargestPossibleRegion().GetSize()[i], rightCorner[i]);
//	}
//
	itk::ImageRegionIterator< AtlasImageType > imgITK(m_itkImage, m_itkImage->GetRequestedRegion());
	itk::ImageRegionIterator< AtlasImageType > img_NOI(m_NOIImage, m_NOIImage->GetRequestedRegion());
	for (imgITK.GoToBegin(), img_NOI.GoToBegin(); !imgITK.IsAtEnd() && !img_NOI.IsAtEnd(); ++imgITK, ++img_NOI)
	{
		AtlasImageType::IndexType currIndex = img_NOI.GetIndex();

		//if (currIndex[0] >= leftCorner[0] && currIndex[0] <= rightCorner[0]
		//	&& currIndex[1] >= leftCorner[1] && currIndex[1] <= rightCorner[1]
		//	&& currIndex[2] >= leftCorner[2] && currIndex[2] <= rightCorner[2])
		//{
		//	if (img_NOI.Get() > -10000) img_NOI.Set(-10000);
		//	else img_NOI.Set(imgITK.Get());
		//}

		if (img_NOI.Get() > -10000) img_NOI.Set(-10000);
		else img_NOI.Set(imgITK.Get());
	}
//
//>>>>>>>>>>>>>>>>>>>>
//
//	std::vector< double > ProbNOIs;
//	this->PriorProb_NOI = 0.0;
//
//	itk::ImageRegionIterator< AtlasImageType > imgNOI(m_NOIImage, m_NOIImage->GetRequestedRegion());
//	for (imgNOI.GoToBegin(); !imgNOI.IsAtEnd(); ++imgNOI)
//	{
//		if (imgNOI.Get() < -1024) continue;
//
//		InternalPixelType gValue = imgNOI.Get();
//
//		double currOutROI = this->getGaussianProbOfROI(gValue);
//		double currOutNOI = this->getGaussianProbOfNOI(gValue);
//		double currProbNN = this->getPriorProbBasedOnNNMap(imgNOI.GetIndex());
//
//		if (_lambda == 0) currProbNN = 0.5;
//
//		double probTotal = currOutROI * currProbNN + currOutNOI * (1 - currProbNN);
//
//		if (probTotal == 0) continue;
//
//		ProbNOIs.push_back(currOutNOI * (1 - currProbNN) / probTotal);
//		this->PriorProb_NOI += currOutNOI * (1 - currProbNN) / probTotal;
//	}
//
//	if (ProbNOIs.size() == 0)
//	{
//		this->PriorProb_NOI = 0.0;
//		this->PriorProb_NOI_STD = 0.0;
//	}
//	else
//	{
//		this->PriorProb_NOI /= ProbNOIs.size();
//
//		double sumNOI = 0.0;
//		for (int i = 0; i < ProbNOIs.size(); i++)
//			sumNOI += pow(ProbNOIs.at(i) - this->PriorProb_NOI, 2);
//
//		this->PriorProb_NOI_STD = sumNOI / ProbNOIs.size();
//		this->PriorProb_NOI_STD = sqrt(this->PriorProb_NOI_STD);
//	}
//

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

	
	//if (InROI >= this->PriorProb_ROI - this->PriorProb_ROI_STD // && InROI <= this->PriorProb_InROI + PROB_INROI_STD
	//	&& InNOI <= this->PriorProb_NOI + this->PriorProb_NOI_STD // && InNOI >= this->PriorProb_InNOI - PROB_INNOI_STD && 
	//	&& InROI - OutROI >= this->PriorProb_ROI_STD && OutNOI - InNOI >= this->PriorProb_NOI_STD)
	//	inside = true;

	//if (OutNOI >= this->PriorProb_NOI - this->PriorProb_NOI_STD // && OutNOI <= this->PriorProb_OutNOI + PROB_OUTNOI_STD
	//	&& OutROI <= this->PriorProb_ROI + this->PriorProb_ROI_STD // && OutROI >= this->PriorProb_OutROI - PROB_OUTROI_STD 
	//	&& OutNOI - InNOI >= this->PriorProb_NOI_STD && InROI - OutROI >= this->PriorProb_ROI_STD)
	//	outside = true;

	//if (inside == true && outside == true) return true; 
	//else return false; 

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


void ImageAnalyser::ComputeEdgeImage(AtlasImageType::Pointer _img)
{
	m_edgeImage = AtlasImageType::New();
	AtlasImageType::RegionType region;
	region.SetIndex(_img->GetLargestPossibleRegion().GetIndex());
	region.SetSize(_img->GetLargestPossibleRegion().GetSize());
	m_edgeImage->SetRegions(region);
	m_edgeImage->Allocate();
	m_edgeImage->SetSpacing(_img->GetSpacing());
	m_edgeImage->SetOrigin(_img->GetOrigin());
	m_edgeImage->FillBuffer(itk::NumericTraits<AtlasImageType::PixelType>::Zero);

	InitEdgeMatrix();

	typedef itk::ConstNeighborhoodIterator< AtlasImageType > ConstNeghbIteratorType;
	typedef itk::ImageRegionIterator< AtlasImageType> IteratorType;

	ConstNeghbIteratorType::RadiusType radius;
	radius.Fill(s / 2);

	ConstNeghbIteratorType it(radius, _img, _img->GetRequestedRegion());
	IteratorType out(m_edgeImage, m_edgeImage->GetRequestedRegion());

	it.GoToBegin();
	out.GoToBegin();

	InternalPixelType data;
	while (!it.IsAtEnd())
	{
		float xEdge = 0;
		float yEdge = 0;
		float zEdge = 0;

		for (int z = 0; z < s; z++)
		{
			for (int y = 0; y < s; y++)
			{
				for (int x = 0; x < s; x++)
				{
					int index = x + it.GetStride(1)*y + it.GetStride(2) *z;
					data = it.GetPixel(index);

					xEdge += data * mX[x][y][z];
					yEdge += data * mY[x][y][z];
					zEdge += data * mZ[x][y][z];
				}
			}
		}

		float intens = sqrt(xEdge*xEdge + yEdge*yEdge + zEdge*zEdge);

		++it;
		if (intens >= 10)
		{
			out.Set(intens);
		}
		else
			out.Set(0);
		//		out.Set(intens);
		++out;
	}

	m_edgeImage->Modified();
}


void ImageAnalyser::InitEdgeMatrix()
{
	float m[s][s][s];

	float sigma = 1;
	float sum = 0;
	for (int z = 0; z < s; z++)
		for (int y = 0; y < s; y++)
			for (int x = 0; x < s; x++){
				int X = x - s / 2; int Y = y - s / 2; int Z = z - s / 2;

				float gauss_derivative = (-X / pow(2 * M_PI, 3.0f / 2.0f)) *  exp(-(X*X + Y*Y + Z*Z) / 2.0f);
				m[x][y][z] = gauss_derivative;
			}

	for (int z = 0; z < s; z++)
		for (int y = 0; y < s; y++)
			for (int x = 0; x < s; x++)	{

				int x_rot_y = z; int y_rot_y = y; int z_rot_y = s - 1 - x;
				int x_rot_z = s - 1 - y; int y_rot_z = x; int z_rot_z = z;

				mX[x][y][z] = m[x][y][z];
				mY[x][y][z] = m[x_rot_z][y_rot_z][z_rot_z];
				mZ[x][y][z] = m[x_rot_y][y_rot_y][z_rot_y];
			}
}


void ImageAnalyser::EnhanceIntensityWithProbMap()
{
	m_maskImage = AtlasImageType::New();
	m_edgeImage = AtlasImageType::New();

	// >>>>>>>> morphology probability image >>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	typedef itk::CastImageFilter< AtlasImageType, CharImageType> CastFilterType;
	CastFilterType::Pointer caster = CastFilterType::New();
	caster->SetInput(this->m_probMapImage);
	caster->Update();

	typedef itk::MaskImageFilter< AtlasImageType, CharImageType, AtlasImageType > MaskFilterType;
	MaskFilterType::Pointer maskFilter = MaskFilterType::New();
	maskFilter->SetInput(m_itkImage);

	AtlasImageType::IndexType minBoundIdx = { 512, 512, 512 };
	AtlasImageType::IndexType maxBoundIdx = { 0, 0, 0 };

	if (this->m_closeRadius > 1)
	{
		typedef itk::BinaryBallStructuringElement< CharImageType::PixelType, CharImageType::ImageDimension > StructuringElement;
		StructuringElement ballClose;
		ballClose.SetRadius(this->m_closeRadius);
		ballClose.CreateStructuringElement();

		typedef itk::BinaryMorphologicalClosingImageFilter< CharImageType, CharImageType, StructuringElement > BinaryMorphologicalClosingImageFilterType;
		BinaryMorphologicalClosingImageFilterType::Pointer closer = BinaryMorphologicalClosingImageFilterType::New();
		closer->SetKernel(ballClose);
		closer->SetInput(caster->GetOutput());
		closer->SetForegroundValue(1);           // necessary 
		closer->UpdateLargestPossibleRegion();   // necessary 
		closer->Update();                        // necessary 

		/*typedef itk::ImageFileWriter< CharImageType > CharWriterType;
		CharWriterType::Pointer charWriter = CharWriterType::New();
		charWriter->SetInput(closer->GetOutput());
		charWriter->SetFileName(m_nameSegImg + "_closer.nrrd");
		charWriter->Update();*/

		itk::ImageRegionIterator< CharImageType > it(closer->GetOutput(), closer->GetOutput()->GetRequestedRegion());
		for (it.GoToBegin(); !it.IsAtEnd(); ++it)
		{
			if (it.Value() > 0)
			{
				CharImageType::IndexType index = it.GetIndex();
				if (index[0] <= minBoundIdx[0]) minBoundIdx[0] = (int)index[0];
				if (index[1] <= minBoundIdx[1]) minBoundIdx[1] = (int)index[1];
				if (index[2] <= minBoundIdx[2]) minBoundIdx[2] = (int)index[2];

				if (index[0] >= maxBoundIdx[0]) maxBoundIdx[0] = (int)index[0];
				if (index[1] >= maxBoundIdx[1]) maxBoundIdx[1] = (int)index[1];
				if (index[2] >= maxBoundIdx[2]) maxBoundIdx[2] = (int)index[2];
			}
		}

		maskFilter->SetMaskImage(closer->GetOutput());
		maskFilter->SetOutsideValue(-10000);
		maskFilter->Update();
	}
	else
	{
		cout << " [ " << __FUNCTION__ << " ] NO Closing , crop the pure mask ... " << endl;

		itk::ImageRegionIterator< CharImageType > it(caster->GetOutput(), caster->GetOutput()->GetRequestedRegion());
		for (it.GoToBegin(); !it.IsAtEnd(); ++it)
		{
			if (it.Value() > 0)
			{
				CharImageType::IndexType index = it.GetIndex();
				if (index[0] <= minBoundIdx[0]) minBoundIdx[0] = (int)index[0];
				if (index[1] <= minBoundIdx[1]) minBoundIdx[1] = (int)index[1];
				if (index[2] <= minBoundIdx[2]) minBoundIdx[2] = (int)index[2];

				if (index[0] >= maxBoundIdx[0]) maxBoundIdx[0] = (int)index[0];
				if (index[1] >= maxBoundIdx[1]) maxBoundIdx[1] = (int)index[1];
				if (index[2] >= maxBoundIdx[2]) maxBoundIdx[2] = (int)index[2];
			}
		}

		maskFilter->SetMaskImage(caster->GetOutput());
		maskFilter->SetOutsideValue(-10000);
		maskFilter->Update();
	}

	// >>>>>>>> extract region from tempImg >>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	AtlasImageType::IndexType extractIndex = { minBoundIdx[0], minBoundIdx[1], minBoundIdx[2] };
	AtlasImageType::SizeType extractSize;
	extractSize[0] = maxBoundIdx[0] - minBoundIdx[0] + 1;
	extractSize[1] = maxBoundIdx[1] - minBoundIdx[1] + 1;
	extractSize[2] = maxBoundIdx[2] - minBoundIdx[2] + 1;

	AtlasImageType::PointType pointOrigin;
	pointOrigin[0] = m_itkImage->GetOrigin()[0];
	pointOrigin[1] = m_itkImage->GetOrigin()[1];
	pointOrigin[2] = m_itkImage->GetOrigin()[2];

	itk::Point< double, 3 > origin;
	origin[0] = minBoundIdx[0] * m_itkImage->GetSpacing()[0] + m_itkImage->GetOrigin()[0];
	origin[1] = minBoundIdx[1] * m_itkImage->GetSpacing()[1] + m_itkImage->GetOrigin()[1];
	origin[2] = minBoundIdx[2] * m_itkImage->GetSpacing()[2] + m_itkImage->GetOrigin()[2];

	AtlasImageType::IndexType startIndex = { 0, 0, 0 };
	//	AtlasImageType::SizeType sizeFull = m_itkImage->GetLargestPossibleRegion().GetSize();
	AtlasImageType::RegionType newRegion(startIndex, extractSize);
	m_maskImage->SetRegions(newRegion);
	m_maskImage->SetOrigin(origin);
	m_maskImage->SetSpacing(m_itkImage->GetSpacing());
	m_maskImage->Allocate();
	m_maskImage->FillBuffer(-10000);


	m_edgeImage->SetRegions(newRegion);
	m_edgeImage->SetOrigin(origin);
	m_edgeImage->SetSpacing(m_itkImage->GetSpacing());
	m_edgeImage->Allocate();
	m_edgeImage->FillBuffer(0);


	AtlasImageType::RegionType extractRegion(extractIndex, extractSize);
	typedef itk::PasteImageFilter< AtlasImageType, AtlasImageType > PasteFilterType;
	PasteFilterType::Pointer paster = PasteFilterType::New();
	paster->SetSourceImage(maskFilter->GetOutput());
	paster->SetSourceRegion(extractRegion);
	paster->SetDestinationImage(m_maskImage);
	paster->SetDestinationIndex(startIndex);
	paster->Update();

	/////////////////////////////////////////////// Added ///////////////////////////////////////////

	/*AtlasImageType::Pointer croppedImg = AtlasImageType::New();
	croppedImg->SetRegions(newRegion);
	croppedImg->SetOrigin(origin);
	croppedImg->SetSpacing(m_itkImage->GetSpacing());
	croppedImg->Allocate();
	croppedImg->FillBuffer(0);

	PasteFilterType::Pointer pasterNew = PasteFilterType::New();
	pasterNew->SetSourceImage(m_itkImage);
	pasterNew->SetSourceRegion(extractRegion);
	pasterNew->SetDestinationImage(croppedImg);
	pasterNew->SetDestinationIndex(startIndex);
	pasterNew->Update();
	croppedImg = pasterNew->GetOutput();

	AtlasImageWriterType::Pointer writerNew = AtlasImageWriterType::New();
	writerNew->SetInput(croppedImg);
	writerNew->SetFileName(this->m_nameSegImg + "_ROI.nrrd");
	writerNew->Update();*/

	//////////////////////////////////////////////// Added //////////////////////////////////////////////////

	m_maskImage = paster->GetOutput();

	//AtlasImageWriterType::Pointer writer = AtlasImageWriterType::New();
	//writer->SetInput(m_maskImage);
	//writer->SetFileName(this->m_nameSegImg + "_ROI_L.nrrd");
	//writer->Update();

	ComputeEdgeImage(m_maskImage);

	// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


	/*AtlasImageWriterType::Pointer writer = AtlasImageWriterType::New();
	writer->SetInput(m_maskImage);
	writer->SetFileName(this->m_nameSegImg + "mask.nrrd");
	writer->Update();*/


	/*AtlasImageWriterType::Pointer writer0 = AtlasImageWriterType::New();
	writer0->SetInput(m_edgeImage);
	writer0->SetFileName(this->m_nameSegImg + "edge.nrrd");
	writer0->Update();*/

	// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	//typedef itk::Function::CustomColormapFunction<CharImageType::PixelType, ColorImageType::PixelType> ColormapType;
	//ColormapType::ChannelType redChannel;
	//ColormapType::ChannelType greenChannel;
	//ColormapType::ChannelType blueChannel;
	//redChannel.push_back(static_cast<ColormapType::RealType>(0.231373));
	//greenChannel.push_back(static_cast<ColormapType::RealType>(0.298039));
	//blueChannel.push_back(static_cast<ColormapType::RealType>(0.752941));
	//redChannel.push_back(static_cast<ColormapType::RealType>(1));
	//greenChannel.push_back(static_cast<ColormapType::RealType>(1));
	//blueChannel.push_back(static_cast<ColormapType::RealType>(1));
	//redChannel.push_back(static_cast<ColormapType::RealType>(0.705882));
	//greenChannel.push_back(static_cast<ColormapType::RealType>(0.0156863));
	//blueChannel.push_back(static_cast<ColormapType::RealType>(0.14902));
	//ColormapType::Pointer colormap = ColormapType::New();
	//colormap->SetRedChannel(redChannel);
	//colormap->SetGreenChannel(greenChannel);
	//colormap->SetBlueChannel(blueChannel);
	//// Setup conversion
	//typedef itk::ScalarToRGBColormapImageFilter<CharImageType, ColorImageType> RGBFilterType;
	//RGBFilterType::Pointer rgbfilter = RGBFilterType::New();
	//rgbfilter->SetInput(m_maskImage);
	//rgbfilter->SetColormap(colormap);
	//typedef itk::ImageFileWriter< ColorImageType > WriterType; 
	//WriterType::Pointer rgbWriter = WriterType::New(); 
	//rgbWriter->SetInput(rgbfilter->GetOutput());
	//rgbWriter->SetFileName("rgbNew.nrrd"); 
	//rgbWriter->Update(); 

}



//>>>>>>>>>>>>>>>>>>> Utilities not in use >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


void ImageAnalyser::SearchPointsOfInterestKidney(vtkPolyData* modelPoly, vtkPolyData* structurePoly, vtkFloatArray* alphaWeights, vtkFloatArray* internalWeights)
{
	vtkSmartPointer< vtkDoubleArray> scalarArray = vtkSmartPointer< vtkDoubleArray >::New();
	scalarArray->SetNumberOfValues(modelPoly->GetNumberOfPoints());

	float edgeThresholdMin;
	float edgeThresholdMax;
	float energyIntern;   // constraint data[0] 
	float energyExtern;   // constraint data[7] 
	float HUmin = -50;    // constraint data[2] - 25% quantile 
	float HUmax = 270;    // constraint data[3] - 75% quantile 
	float gradientRange;  // constraint data[4] - data[5] 
	float HURange;

	double* tensorData;

	vtkPointData* modelPointData = modelPoly->GetPointData();
	vtkDataArray* constraints = modelPointData->GetTensors();  /**> constrains for every point **/

	vtkPoints* modelPoints = modelPoly->GetPoints();
	vtkPoints* structurePoints = structurePoly->GetPoints();
	vtkDataArray* normals = modelPoly->GetPointData()->GetNormals();

	float highestIntens = 0.0;       // highest intensity of m_edgeImage 

	for (long ptId = 0; ptId < modelPoints->GetNumberOfPoints(); ptId++)
	{
		tensorData = constraints->GetTuple(ptId); // set constraints is important 

		gradientRange = tensorData[0];
		HURange = tensorData[1];  // HU range 
		HUmin = tensorData[2];
		HUmax = tensorData[3];
		edgeThresholdMin = tensorData[4];
		edgeThresholdMax = tensorData[5];
		energyIntern = tensorData[0] / 5.0;
		energyExtern = tensorData[7];

		if (HURange > 70)
		{
			edgeThresholdMin = 5;
		}

		double point[3];
		modelPoints->GetPoint(ptId, point);

		double normal[3];
		normals->GetTuple(ptId, normal);

		double result[3];  /**> the points moving to **/
		result[0] = point[0]; result[1] = point[1]; result[2] = point[2];

		highestIntens = 0.0;       // highest intensity of m_edgeImage 
		float predecessorIntens = 0.0;   // previous 
		float successorIntens = 0.0;     // later 

		double predecessorPoint[3];
		predecessorPoint[0] = point[0]; predecessorPoint[1] = point[1]; predecessorPoint[2] = point[2];

		double successorPoint[3];
		successorPoint[0] = point[0]; successorPoint[1] = point[1]; successorPoint[2] = point[2];


		bool meshPointIsInside;
		bool breakOuterSearch = false;
		bool breakInnerSearch = false;


		double currentPointIntensity;
		getPixelIntensity(point, currentPointIntensity, m_maskImage);/*m_itkImage*/


		if (pointIsInside(point, HUmin, HUmax) && currentPointIntensity > 10)
		{
			meshPointIsInside = true;
		}
		else
		{
			meshPointIsInside = false;
		}

		double intens = 0.0f;
		int maxStep = 5;

		bool fixed = false;

		for (int i = 0; i < maxStep; i++)   // search along the normal of the current mesh point
		{

			// >>>>>>>>>>>>>>>>>>>>>>> SuccessorPoint Search Start >>>>>>>>>>>>>>>>>>>>>>>>>>

			getPixelIntensity(successorPoint, intens, m_edgeImage); // get the edge intensity 
			successorIntens = intens;

			if (pointIsOnBoundary(successorPoint))
			{
				result[0] = successorPoint[0]; result[1] = successorPoint[1]; result[2] = successorPoint[2];
				energyExtern = 0.0;
				breakOuterSearch = true;
				break;
			}

			if (meshPointIsInside && !breakOuterSearch)
			{
				if ((successorIntens > edgeThresholdMin) && (successorIntens < edgeThresholdMax) && (successorIntens > highestIntens) && pointIsInside(successorPoint, HUmin, HUmax))
				{
					result[0] = successorPoint[0]; result[1] = successorPoint[1]; result[2] = successorPoint[2];
					highestIntens = successorIntens;	// update highest edge intensity found so far
				}
				else if (pointIsInside(successorPoint, HUmin, HUmax) && highestIntens == 0)   // successorIntens <= 0
				{
					if (energyExtern != 0 && HURange < 70) //  no valid edge was found so far 
					{
						result[0] = successorPoint[0]; result[1] = successorPoint[1]; result[2] = successorPoint[2];
					}
				}
				else if (highestIntens == 0) // point does not change 
				{
					result[0] = point[0]; result[1] = point[1]; result[2] = point[2];
					energyExtern = 0.0;
					breakOuterSearch = true;  // then never search outer surface 
				}
			}

			successorPoint[0] += normal[0];		successorPoint[1] += normal[1];		successorPoint[2] += normal[2];

			// >>>>>>>>>>>>>>>>>>>>>>> SuccessorPoint Search End >>>>>>>>>>>>>>>>>>>>>>>>>>>

			predecessorPoint[0] -= normal[0];		predecessorPoint[1] -= normal[1];		predecessorPoint[2] -= normal[2];

			getPixelIntensity(predecessorPoint, intens, m_edgeImage);

			predecessorIntens = intens;

			if (!breakInnerSearch)
			{
				if (meshPointIsInside)
				{
					if (!pointIsInside(predecessorPoint, HUmin, HUmax)) // this inner searched point is not inside 
					{
						energyExtern = 0.0; // do not search outer ? 
					}
				}
				else
				{
					if ((predecessorIntens > edgeThresholdMin) && (predecessorIntens < edgeThresholdMax) && (predecessorIntens > highestIntens) && pointIsInside(predecessorPoint, HUmin, HUmax))
					{
						result[0] = predecessorPoint[0]; result[1] = predecessorPoint[1]; result[2] = predecessorPoint[2];
						highestIntens = predecessorIntens;
					}
					else if (highestIntens == 0) // no valid edge has been found
					{
						result[0] = predecessorPoint[0]; result[1] = predecessorPoint[1]; result[2] = predecessorPoint[2];
					}
				}

			} // previous break inner search 

		} // deformation done  

		energyIntern /= 2.0; // normalize internal weight

		internalWeights->InsertComponent(ptId, 0, energyIntern);
		alphaWeights->InsertComponent(ptId, 0, energyExtern);

		structurePoints->SetPoint(ptId, result[0], result[1], result[2]);

		//		scalarArray->SetValue(ptId, (double)(ptId * 100 / modelPoly->GetNumberOfPoints())); 

		if (pointIsOnBoundary(result))
		{
			scalarArray->SetValue(ptId, 100);
		}
		else if (pointIsInside(result, HUmin, HUmax))
		{
			scalarArray->SetValue(ptId, 20.0);
		}
		else
		{
			scalarArray->SetValue(ptId, 80.0);
		}

	} // All point iteration done 

	modelPoly->GetPointData()->SetScalars(scalarArray);
	structurePoly->GetPointData()->SetScalars(scalarArray);

	structurePoints->Modified();
	structurePoly->Modified();
	alphaWeights->Modified();
	internalWeights->Modified();


	cout << __FUNCTION__ << " number of points inside = " << endl;

}


void ImageAnalyser::AutoSetKidneyModelConstraints(mitk::Surface::Pointer surface)
{

	if (!m_mitkImage || !surface) return;	// make sure an image and model have been loaded

	vtkSmartPointer<vtkPolyData> modelPoly = surface->GetVtkPolyData();  // get polydata ..
	vtkSmartPointer<vtkPoints> modelPoints = modelPoly->GetPoints();     // .. and points

	
	statisticFilterType::Pointer statFilter = statisticFilterType::New();


	mitk::SurfaceToImageFilter::Pointer surfaceToImageFilter = mitk::SurfaceToImageFilter::New();
	surfaceToImageFilter->SetInput(surface);
	surfaceToImageFilter->SetImage(m_mitkImage);
	surfaceToImageFilter->Update();

	mitk::Image::Pointer outputImage = surfaceToImageFilter->GetOutput();  // only the segmented part is used 

	mitk::ImageTimeSelector::Pointer timeSelector = mitk::ImageTimeSelector::New();
	timeSelector->SetInput(outputImage);
	timeSelector->SetTimeNr(0);
	timeSelector->UpdateLargestPossibleRegion();
	timeSelector->Update();

	AtlasImageType::Pointer timeSelectorOutputItkImage = AtlasImageType::New();
	mitk::CastToItkImage(timeSelector->GetOutput(), timeSelectorOutputItkImage);
	statFilter->SetInput(timeSelectorOutputItkImage);
	statFilter->Update();


	int quantile25 = statFilter->Get25PercentQuantile();
	int quantile75 = statFilter->Get75PercentQuantile();

	long numPts = modelPoints->GetNumberOfPoints();		// number of model points
	long ptId;

	vtkSmartPointer<vtkPointData> pointData = modelPoly->GetPointData();
	vtkSmartPointer<vtkDataArray> constraints = pointData->GetTensors();

	int constraintsCount = 9;

	if (constraints == 0)  // if no constraints found, add constraints data to poly
	{
		cout << " [ " << __FUNCTION__ << " ] Set default constraints (constraints = 0)." << endl;
		constraints = vtkSmartPointer<vtkFloatArray>::New();
		constraints->SetNumberOfComponents(constraintsCount);

		for (int i = 0; i < numPts; i++)
		{
			double data[9];
			for (int componentNr = 0; componentNr < constraintsCount; componentNr++)
			{
				data[componentNr] = 0.0;
			}

			// set some defaults
			data[0] = 5.0;		// energy intern
			data[4] = 5.0;		// edgeMin
			data[5] = 2200.0;	// edgeMax
			data[7] = 2.0;		// energy extern

			constraints->InsertTuple(i, data);
		}
		constraints->Modified();
		pointData->SetTensors(constraints);
	}

	// define min / max settings for the desired organ
	double fixedMin = DBL_MIN;
	double fixedMax = DBL_MAX;
	double currMin = (double)statFilter->Get5PercentQuantile();
	double currMax = (double)statFilter->GetMaximum(); //(double)quantile75;

	fixedMin = (double)statFilter->Get5PercentQuantile();
	fixedMax = (double)statFilter->GetMaximum();
	currMax = fixedMax;

	// iterate over all points of the surface and set the min/max HU values
	for (ptId = 0; ptId < numPts; ptId++)
	{
		double* tdata = constraints->GetTuple(ptId);

		//// set HU min and max values to the quantiles calculated from the model position
		//if (currMin > fixedMin && currMin < fixedMax)
		//{
		//	tdata[2] = currMin; // get the larger minHU
		//}
		//else
		//{
		//	tdata[2] = fixedMin;
		//}
		//if (currMax > fixedMin && currMax < fixedMax)
		//{
		//	tdata[3] = currMax;   // get the smaller maxHU 
		//}
		//else
		//{
		//	tdata[3] = fixedMax;
		//}

		tdata[2] = std::min(currMin, -(double)20.0);
		tdata[3] = currMax;

		constraints->SetTuple(ptId, tdata);
	}

	constraints->SetName("auto kidney constraints");

	cout << " [ " << __FUNCTION__ << " ] constraint HU - [ " << min(currMin, -40.0) << " , " << currMax << " ] " << endl;
}


void ImageAnalyser::AutoSetPancreasModelConstraints(mitk::Surface::Pointer surface)
{
	if (!m_mitkImage || !surface) return;	// make sure an image and model have been loaded

	vtkSmartPointer<vtkPolyData> modelPoly = surface->GetVtkPolyData();  // get polydata ..
	vtkSmartPointer<vtkPoints> modelPoints = modelPoly->GetPoints();     // .. and points

	typedef itk::VolumeGreyLevelStatisticsFilter< AtlasImageType, AtlasImageType> statisticFilterType;
	statisticFilterType::Pointer statFilter = statisticFilterType::New();

	mitk::SurfaceToImageFilter::Pointer surfaceToImageFilter = mitk::SurfaceToImageFilter::New();
	surfaceToImageFilter->SetInput(surface);
	surfaceToImageFilter->SetImage(m_mitkImage);
	surfaceToImageFilter->Update();

	mitk::Image::Pointer outputImage = surfaceToImageFilter->GetOutput();  // only the segmented part is used 

	mitk::ImageTimeSelector::Pointer timeSelector = mitk::ImageTimeSelector::New();
	timeSelector->SetInput(outputImage);
	timeSelector->SetTimeNr(0);
	timeSelector->UpdateLargestPossibleRegion();
	timeSelector->Update();

	AtlasImageType::Pointer timeSelectorOutputItkImage = AtlasImageType::New();
	mitk::CastToItkImage(timeSelector->GetOutput(), timeSelectorOutputItkImage);
	statFilter->SetInput(timeSelectorOutputItkImage);
	statFilter->Update();


	int quantile25 = statFilter->Get25PercentQuantile();
	int quantile75 = statFilter->Get75PercentQuantile();

	long numPts = modelPoints->GetNumberOfPoints();		// number of model points
	long ptId;

	vtkSmartPointer<vtkPointData> pointData = modelPoly->GetPointData();
	vtkSmartPointer<vtkDataArray> constraints = pointData->GetTensors();

	int constraintsCount = 9;

	if (constraints == 0)  // if no constraints found, add constraints data to poly
	{
		cout << " [ " << __FUNCTION__ << " ] Set default constraints (constraints = 0)." << endl;
		constraints = vtkSmartPointer<vtkFloatArray>::New();
		constraints->SetNumberOfComponents(constraintsCount);

		for (int i = 0; i < numPts; i++)
		{
			double data[9];
			for (int componentNr = 0; componentNr < constraintsCount; componentNr++)
			{
				data[componentNr] = 0.0;
			}

			// set some defaults
			data[0] = 5.0;		// energy intern
			data[4] = 10.0;		// edgeMin
			data[5] = 40.0;	// edgeMax
			data[7] = 2.0;		// energy extern

			constraints->InsertTuple(i, data);
		}
		constraints->Modified();
		pointData->SetTensors(constraints);
	}

	// define min / max settings for the desired organ
	double fixedMin = DBL_MIN;
	double fixedMax = DBL_MAX;
	double currMin = (double)quantile25;
	double currMax = (double)quantile75;

	fixedMin = (double)quantile25;
	fixedMax = (double)quantile75;
	currMax = fixedMax;

	// iterate over all points of the surface and set the min/max HU values
	for (ptId = 0; ptId < numPts; ptId++)
	{
		double* tdata = constraints->GetTuple(ptId);

		tdata[2] = 50; //std::min(currMin, (double)70.0);
		tdata[3] = 120; //std::min(currMax, (double)120.0);

		constraints->SetTuple(ptId, tdata);
	}

	constraints->SetName("auto kidney constraints");

}


void ImageAnalyser::SearchPointsOfInterestPancreas(vtkPolyData* modelPoly, vtkPolyData* structurePoly, vtkFloatArray* alphaWeights, vtkFloatArray* internalWeights)
{
	cout << __FUNCTION__ << endl;

	vtkSmartPointer< vtkDataArray > weights = vtkSmartPointer< vtkDoubleArray >::New();
	if (!modelPoly->GetPointData()->GetArray("weights"))
	{
		MITK_WARN << __FUNCTION__ << " no weights property found ... ";
		return;
	}

	weights = modelPoly->GetPointData()->GetArray("weights");


	float edgeThresholdMin;
	float edgeThresholdMax;
	float energyIntern;   // constraint data[0] 
	float energyExtern;   // constraint data[7] 
	float HUmin = 50;    // constraint data[2] - 25% quantile 
	float HUmax = 120;    // constraint data[3] - 75% quantile 
	float gradientRange;  // constraint data[4] - data[5] 
	float HURange;

	double* tensorData;

	vtkPointData* modelPointData = modelPoly->GetPointData();
	vtkDataArray* constraints = modelPointData->GetTensors();  /**> constrains for every point **/

	vtkPoints* modelPoints = modelPoly->GetPoints();
	vtkPoints* structurePoints = structurePoly->GetPoints();
	vtkDataArray* normals = modelPoly->GetPointData()->GetNormals();

	float highestIntens = 0.0;       // highest intensity of m_edgeImage 

	for (long ptId = 0; ptId < modelPoints->GetNumberOfPoints(); ptId++)
	{
		tensorData = constraints->GetTuple(ptId); // set constraints is important 

		gradientRange = tensorData[0];
		HURange = tensorData[1];  // HU range 
		HUmin = tensorData[2];
		HUmax = tensorData[3];
		edgeThresholdMin = tensorData[4];
		edgeThresholdMax = tensorData[5];
		energyIntern = tensorData[0] / 5.0;
		energyExtern = tensorData[7];

		if (HURange > 40)
		{
			edgeThresholdMin = 5;
		}

		double point[3];
		modelPoints->GetPoint(ptId, point);

		double normal[3];
		normals->GetTuple(ptId, normal);

		double result[3];  /**> the points moving to **/
		result[0] = point[0]; result[1] = point[1]; result[2] = point[2];

		highestIntens = 0.0;       // highest intensity of m_edgeImage 
		float predecessorIntens = 0.0;   // previous 
		float successorIntens = 0.0;     // later 

		double predecessorPoint[3];
		predecessorPoint[0] = point[0]; predecessorPoint[1] = point[1]; predecessorPoint[2] = point[2];

		double successorPoint[3];
		successorPoint[0] = point[0]; successorPoint[1] = point[1]; successorPoint[2] = point[2];


		bool meshPointIsInside;
		bool breakOuterSearch = false;
		bool breakInnerSearch = false;


		double currentPointIntensity;
		getPixelIntensity(point, currentPointIntensity, m_itkImage);


		if (pointIsInside(point, HUmin, HUmax) && currentPointIntensity > 10)
		{
			meshPointIsInside = true;
		}
		else
		{
			meshPointIsInside = false;
		}

		double intens = 0.0f;
		int maxStep = 5;
		double baseStep = 1.0;

		bool fixed = false;

		for (int i = 0; i < maxStep; i++)   // search along the normal of the current mesh point
		{

			// >>>>>>>>>>>>>>>>>>>>>>> SuccessorPoint Search Start >>>>>>>>>>>>>>>>>>>>>>>>>>

			getPixelIntensity(successorPoint, intens, m_edgeImage); // get the edge intensity 
			successorIntens = intens;

			if (pointIsOnBoundary(successorPoint))
			{
				result[0] = successorPoint[0]; result[1] = successorPoint[1]; result[2] = successorPoint[2];
				energyExtern = 0.0;
				breakOuterSearch = true;
				break;
			}

			if (meshPointIsInside && !breakOuterSearch)
			{
				if ((successorIntens > edgeThresholdMin) && (successorIntens < edgeThresholdMax) && (successorIntens > highestIntens) && pointIsInside(successorPoint, HUmin, HUmax))
				{
					result[0] = successorPoint[0]; result[1] = successorPoint[1]; result[2] = successorPoint[2];
					highestIntens = successorIntens;	// update highest edge intensity found so far
				}
				else if (pointIsInside(successorPoint, HUmin, HUmax) && highestIntens == 0)   // successorIntens <= 0
				{
					if (energyExtern != 0 && HURange < 70) //  no valid edge was found so far 
					{
						result[0] = successorPoint[0]; result[1] = successorPoint[1]; result[2] = successorPoint[2];
					}
				}
				else if (highestIntens == 0) // point does not change 
				{
					result[0] = point[0]; result[1] = point[1]; result[2] = point[2];
					energyExtern = 0.0;
					breakOuterSearch = true;  // then never search outer surface 
				}
			}

			successorPoint[0] += normal[0];		successorPoint[1] += normal[1];		successorPoint[2] += normal[2];

			// >>>>>>>>>>>>>>>>>>>>>>> SuccessorPoint Search End >>>>>>>>>>>>>>>>>>>>>>>>>>>

			predecessorPoint[0] -= normal[0];		predecessorPoint[1] -= normal[1];		predecessorPoint[2] -= normal[2];

			getPixelIntensity(predecessorPoint, intens, m_edgeImage);

			predecessorIntens = intens;

			if (!breakInnerSearch)
			{
				if (meshPointIsInside)
				{
					if (!pointIsInside(predecessorPoint, HUmin, HUmax)) // this inner searched point is not inside 
					{
						energyExtern = 0.0; // do not search outer ? 
					}
				}
				else
				{
					if ((predecessorIntens > edgeThresholdMin) && (predecessorIntens < edgeThresholdMax) && (predecessorIntens > highestIntens) && pointIsInside(predecessorPoint, HUmin, HUmax))
					{
						result[0] = predecessorPoint[0]; result[1] = predecessorPoint[1]; result[2] = predecessorPoint[2];
						highestIntens = predecessorIntens;
					}
					else if (highestIntens == 0) // no valid edge has been found
					{
						result[0] = predecessorPoint[0]; result[1] = predecessorPoint[1]; result[2] = predecessorPoint[2];
					}
				}

			} // previous break inner search 

		} // deformation done  

		energyIntern /= 3.0; // normalize internal weight

		internalWeights->InsertComponent(ptId, 0, energyIntern);
		alphaWeights->InsertComponent(ptId, 0, energyExtern);

		structurePoints->SetPoint(ptId, result[0], result[1], result[2]);

		//		scalarArray->SetValue(ptId, (double)(ptId * 100 / modelPoly->GetNumberOfPoints())); 

		if (pointIsInside(result, HUmin, HUmax))
		{
			weights->SetTuple1(ptId, 1.0);
		}
		else
		{
			weights->SetTuple1(ptId, 0.0);
		}

	} // All point iteration done 

	weights->Modified(); 
	modelPoly->GetPointData()->SetScalars(weights);
	structurePoly->GetPointData()->SetScalars(weights);

	structurePoints->Modified();
	structurePoly->Modified();
	alphaWeights->Modified();
	internalWeights->Modified();


	cout << __FUNCTION__ << endl;

}


inline bool ImageAnalyser::pointIsInside(double point[3], float& HUmin, float& HUmax)
{
	typedef itk::ImageRegionConstIterator< AtlasImageType > ConstIteratorType;

	AtlasImageType::IndexType start;
	itk::Point<double, 3> itkPoint;

	itkPoint[0] = point[0]; itkPoint[1] = point[1]; itkPoint[2] = point[2];

	if (!m_itkImage->TransformPhysicalPointToIndex(itkPoint, start)) return false;


	int s_neighbourhood = 3;

	start[0] -= s_neighbourhood / 2;
	start[1] -= s_neighbourhood / 2;
	start[2] -= s_neighbourhood / 2;

	AtlasImageType::SizeType size;
	size[0] = s_neighbourhood;
	size[1] = s_neighbourhood;
	size[2] = s_neighbourhood;

	for (int i = 0; i < 3; i++)  // important ! 
	{
		if (start[i] < m_itkImage->GetLargestPossibleRegion().GetIndex()[i] || start[i] + size[i] > m_itkImage->GetLargestPossibleRegion().GetSize()[i])
		{
			return false;
		}
	}


	AtlasImageType::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);

	ConstIteratorType it(m_itkImage, desiredRegion);
	it.GoToBegin();

	InternalPixelType gValue = 0.0f;

	int count = 0;

	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		gValue = it.Get();
		if ((gValue <= HUmax) && (gValue >= HUmin))
		{
			count++;
		}
	}

	if (count >= ((s_neighbourhood*s_neighbourhood*s_neighbourhood) / 4))
		return true;
	else
		return false;
}


inline bool ImageAnalyser::pointIsOnBoundary(double point[3])
{
	typedef itk::ImageRegionConstIterator< AtlasImageType > ConstIteratorType;

	AtlasImageType::IndexType start;
	itk::Point<double, 3> itkPoint;

	itkPoint[0] = point[0]; itkPoint[1] = point[1]; itkPoint[2] = point[2];

	if (!m_itkImage->TransformPhysicalPointToIndex(itkPoint, start)) return false;

	double intens = 0.0;
	getPixelIntensity(point, intens, m_edgeImage);

	int s_neighbourhood = 5;

	start[0] -= s_neighbourhood / 2;
	start[1] -= s_neighbourhood / 2;
	start[2] -= s_neighbourhood / 2;

	AtlasImageType::SizeType size;
	size[0] = s_neighbourhood;
	size[1] = s_neighbourhood;
	size[2] = s_neighbourhood;

	for (int i = 0; i < 3; i++)
	{
		if (start[i] < m_edgeImage->GetLargestPossibleRegion().GetIndex()[i] || start[i] + size[i] > m_edgeImage->GetLargestPossibleRegion().GetSize()[i])
		{
			return false;
		}
	}


	AtlasImageType::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);

	ConstIteratorType it(m_itkImage, desiredRegion);
	it.GoToBegin();

	InternalPixelType gValue = 0.0f;

	int count = 0;

	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		gValue = it.Get();
		if (gValue > 10)
		{
			count++;
		}
	}

	if ((count >= (s_neighbourhood*s_neighbourhood*s_neighbourhood) * 0.45) && (count <= (s_neighbourhood*s_neighbourhood*s_neighbourhood) * 0.55))
		return true;
	else
		return false;
}


inline void ImageAnalyser::getPixelIntensity(double point[3], double& intensity, AtlasImageType::Pointer image)
{
	AtlasImageType::IndexType pointIndex;

	itk::Point<double, 3> itkPoint;

	itkPoint[0] = point[0]; itkPoint[1] = point[1]; itkPoint[2] = point[2];
	if (image->TransformPhysicalPointToIndex(itkPoint, pointIndex))
	{
		intensity = image->GetPixel(pointIndex);
	}
	else
	{
		intensity = 0.0;
	}
}
