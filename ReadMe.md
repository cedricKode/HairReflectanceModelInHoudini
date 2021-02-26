# A User Friendly Hair Reflectance Model in Houdini

High-quality realistic rendering of hair and fur is a long-standing goal in computer graphics. It includes a number of challenges, such as computation of accurate light scattering for individual fibres and handling the complex geometry of hair volume. Physically based approach in rendering hair is aims to achieve plausible result; however, its drawback is the lack of artist friendly controls. 

In this paper we present an implementation of a user-friendly hair reflectance model in Houdini based on proposals in (Chiang et al., 2016) and (Yan et al., 2017), who constructed physically based hair and fur shading models. The results of these approaches are illustrated in numerous renderings. All images were rendered using the physically plausible pipeline in Mantra.

## Requirements:

- SideFX Houdini (https://www.sidefx.com)
- Visual Studio Code (https://code.visualstudio.com)
or Sublime Text 3 (https://www.sublimetext.com/3)

## Folder Structure

- *cyHairConvert*. This folder includes QtCreator project with code for *cyHairConvert* command line tool. The compiled tool itself also may be found in this folder. The *models* folder include the initial binary hair models. 

- *cyHairRead*. This folder includes a scene with example of loading cyHair converted files via *cyHairLoad* asset. The internal Python code, which is used in *cyHairLoad* asset, is also presented as a separate *hairRead.py* file.

- *DKozlova_Master_Project_Thesis*. This is a pdf version of final thesis paper.

- *configs*. This is the main folder, which includes fur and hair shader nodes, and code for them. The code is stored in *vex/CVex* subfolder. *.vfl* files could be opened with any text editor, for example Visual Studio Code or Sublime Text.

- *index*. This is an *html* code for NCCA website.

- *ReadMe*. A file with explanation of project's structure.

## References

Chiang, M. J.-Y., Bitterli, B., Tappan, C., and Burley, B. 2016. A practical and controllable hair and fur model for production path tracing. In Proceedings of the 37th Annual Conference of the European Association for Computer Graphics, EG ’16, pp. 275–283, Goslar Germany, Germany. Eurographics Association.

Yan, L.-Q., Jensen, H. W., and Ramamoorthi, R. 2017. An efficient and practical near and far field fur reflectance model. ACM Trans. Graph. 36:67:1–67:13.
