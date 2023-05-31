SpatialStegoDectect:  A CNN for Detecting Images with Secret Data in Spatial Domain

Introduction

SpatialStegoDetect is an innovative software proposed to significantly improve the accuracy of detecting concealed information within digital images. The software operates through three distinct phases: pre-processing, feature extraction, and classification. In the pre-processing phase, spatial rich model filters are employed to enhance image quality by reducing noise artifacts. Feature extraction is performed using two-dimensional convolutions, enabling the software to capture intricate details from the images. Classification is accomplished through the utilization of fully connected layers and a softmax function.
	It is worth noting that we use Syndrome Trellis Codes with Matlab to embed the secret data in the cover image. I the folder “Matlab files”, you will find the code to embed the secret data. (See: GitHub Link)
	 It also worth noting that here we avail eight files classified in three categories (these files can also be found in ………  in a folder called “Python_codes”): One file for the SpatialStegoDectect model (SpatialStegoDectect.ipynb), six files for the datasets (X_training.npy, y_training.npy, X_testing.npy, y_testing.npy, X_validating.npy, y_validating.npy), and one file for the spatial rich models (SRM) used (SRM.npy).
Running the Model

To run this model, follow these steps: 
1.	Load all the eight files (SpatialStegoDectect.ipynb, X_training.npy, y_training.npy, X_testing.npy, y_testing.npy, X_validating.npy, y_validating.npy, and SRM.npy) in a same folder.
2.	Create the necessary folder for the trained models and other necessary folders as given in “FuzConvTrue.ipynb” ("/SpatialStegoDectect/trained_models/", and “/SpatialStegoDectect/Results/").
3.	To run the CNN, you need to have Python with MatplotLib, TensorFlow, Numpy, and OpenCV2 libraries and execute SpatialStegoDectect.ipynb

Thank you for considering our code and giving suggestions on how to make it more user-friendly!
