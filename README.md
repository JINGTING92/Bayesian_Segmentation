A Bayesian Model Incorporating Deep Neural Network with Statistical Shape Model for Segmentation

The project is built on the public source MITK, a unified framework consisting of ITK, VTK with QT interface. 

SegmentationView.cpp: the entrance of the program, which connects with a QT interface. 
ImageAnalyser.cpp: the main codes are here.

############################################################################################################################

<SegmentationView.cpp>

Pre-training functions:

A. statistical shape modeling (optional): 
   SegmentationView::LoadDatasetsForModeling() 
   SegmentationView::SelectModelType()
   SegmentationView::TrainShapeModel() 
   
B. load test image to be segmented: 
   SegmentationView::ImageSelected() : denote as m_segImage 
   
############################################################################################################################
   
Step 1: load the probability score map and initialize the shape model 

run SegmentationView::Deform() 

(1) loadProbShapeForLocalization: extract a triangle mesh from the probability socre map (to initialize the shape model)
(2) loadDeformShapeFromLocal(): if the shape model does not exist, please load a mean shape instead as a template shape
(3) initialize the weights: as we assign each landmark a weight to represent its probability to be ROI and NOI
(4) alignDeformShapeWithProbShape: this step is used to initialize the shape model (the template shape)
    transformDeformShapeToProbMap: in order to fit the template to the probability score map  
    
(5) create an ImageAanlyser to perform the Bayesian model 

    ImageAnalyser::trainPriorGMM : initialize the GMM 
    ImageAnalyser::computePriorProbs : initialize the probabilities 
    
(6) re-initialize the window of MITK (optional)

##########################################################################################################################

Step 2: Perform the Bayesian model for segmentation 

SegmentationView::AutoOrientationWithoutModel() 

(1) convert the current surface from the window as the shape to be computed 
(2) AutoSearchBasedOnIntensModel (just for test, not in use in practice, as you can define the number of iterations as well as other variables here)

(3) AutoOrientationWithModel() (optional): this is used for model back projection 
(4) AutoSearchBasedOnIntensModel() : this is the main part 

    In this function, it starts with the transformation from the current surface to the shape to be focused on 
    
    ImageAnalyser::landmarkDeformation() is the key function 
    
(5) Segment() : save the segment 

