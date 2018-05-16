# mediatrix
Mediatrix Filamentation is an iterative algorithm which automaticaly extract features in thin and elongated objects, particulary suited for curved shapes. It can estimate the length, width and curvature center.
This repo contains the code used in the following papers:

https://arxiv.org/abs/1212.1799

https://arxiv.org/abs/1607.04644

For an example of usage see mediatrix_example.ipynb. The dependencies are:
   * numpy
   * aplpy
   * scipy
   * matplotlib
   * astropy
   * pythn-opencv (cv2, if you want support for png, tiff or any other non fits image format)

This is a standalone version that uses functions derived from bit and sltools libs developed by Clécio De Bom (cdebom) , Carlos H. Brandt (chbrandt) and Martin Makler (martinmakler). The authors would like to thank Cristina Furlanetto.

