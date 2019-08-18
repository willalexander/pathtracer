# pathtracer
Simple unidirectional path tracer written in C++.

![](/images/ddagonChrome.bmp){width:50%} ![](/images/dradgonGlass.bmp)

<img src="/images/dragonChrome.bmp"><img src="/images/dragonGlass.bmp">

### Supports:
* Area lights
* Sphere lights
* OBJ polygon geometry
* Lambert diffuse surfaces
* Specular Surfaces
* Homogeneous participating media
* Homogoenous subsurface scattering

### Dependencies:
* Boost
* TinyOBJ
* Embree

### Usage:
```
myRenderer -i scene.waff --psamples 64 --hsamples 1 --lsamples 1 -bmp image.bmp
```

Parameter | Meaning
----------|----------
i | The scene to render (a text file in .waff format. See exmaples in the repo)
psamples | Number of samples per pixel
hsamples | Number of hemisphere samples for indirect diffuse reflection
lsamples | Number of samples for each light source for direct lighting computation
bmp | Name of BMP image file for output
