# Haze Airlight Color Estimation

Matlab code implementing the Airlight color estimation part of the [ICCP 2016 paper "Blind Dehazing Using Internal Patch Recurrence"](http://www.wisdom.weizmann.ac.il/~vision/BlindDehazing/blindDehazing_ICCP2016.pdf).
Tested on Matlab 2016a working on Windows.

If you find our work useful in your research or publication, please cite it:

```
@inproceedings{bahat2016blind,
  title={Blind dehazing using internal patch recurrence},
  author={Bahat, Yuval and Irani, Michal},
  booktitle={Computational Photography (ICCP), 2016 IEEE International Conference on},
  pages={1--9},
  year={2016},
  organization={IEEE}
}
```
----------

## Functions:
### 1. AirlightUsingPatchRecurrence(orgImage):
The core function estimating the haze airlight color in the input image. The input image should be in the form of 3 (RGB) channels (double, values in the range of [0,1]).
### 2. testAirlightEstimation:
A script for testing the airlight estimation. Estimated airlight is presented alongside the corresponding image. This script also displays a ground truth airlight color, either saved in the GTairlights folder or manually extracted by the user.
### 3. estimate4allImages:
Serially computing the airlight color for all images in the 'images' sub-folder. Comparing the estimated airlights with the pre-determined GT airlights for images that contain a region of sky, and presenting the mean distance in RGB space.


## Further comments:
The 'images' sub-folder contain example images, and the 'GTairlights' sub-folder contains ground truth airlight colors for images containing regions of sky.
This code uses the nearest neighbors search engine by Shai Bagon, named ANN. For convinience, it appears as-is in the ann_wrapper sub-folder.

The code is provided as-is for academic use only and without any guarantees. Please contact the author to report any bugs. Written by [Yuval Bahat](http://www.wisdom.weizmann.ac.il/~ybahat/).
