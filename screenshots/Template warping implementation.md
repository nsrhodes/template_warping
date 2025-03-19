---
up:
  - "[[Projects]]"
prev:
  - "[[Template warping paper]]"
---
[[Linux processing PC]]


# Instructions
Tested only on Ubuntu 22.03

### 1. Install fsl
https://fsl.fmrib.ox.ac.uk/fsl/docs/#/
### 2. Install Matlab
Requires vision toolbox to write point clouds
tested using Matlab 2024b

Tested using fieldtrip-20250114 https://www.fieldtriptoolbox.org/download/
### 3. Install Meshlab
Tested using MeshLab 2023.12
 
Download AppImage from https://www.meshlab.net/#download

Open a terminal and make it executable:

```shell
chmod a+x path/to/MeshLab2023.12-linux.AppImage
```

and run it:
```shell
cd path/to
./MeshLab2023.12-linux.AppImage
```

If you get an error regarding missing "FUSE" requirement, install libfuse2 using
```shell
sudo apt install libfuse2
```

## User guide

1. Amend ##script## to include path to your fieldtrip installation and the directory containing the template MRI and the 3D mesh of the subject's head.
2. Amend any options (whether to pad the template MRI or crop the neck line, we recommend the defaults).
3. Run ##script##
4. When prompted, select the template MRI (here ./example/Adult_template.nii.gz)
	A head-segmentation and generation of the scalp mesh will follow. It is likely that you will have to specify the coordinate system during this step. When prompted in the command window, enter the correct axis directions. You will see a figure showing MRI slices .
	**specify_coordinate_sys.png**
	And the following prompt:
	```octave
	The coordinate system is not specified.
	Do you want to change the anatomical labels for the axes [Y, n]? 
	```
5. In this example we enter 'Y', 'r', 'a', 's', 'n', for the following prompts as it is a RAS coordinate system.
6. Wait for the scalp extraction to complete
7. You will be prompted with instructions on how to perform a rough alignment in Meshlab
8. Open a new terminal (ctrl+alt+T) and start Meshlab (see above)
9. Load the generated template MRI scalp
   (saved as <path_to_template_mri_file>\_template_points.stl)
   and the head mesh (here sub-101_head.ply)
	If the head mesh is in mm, the MRI mesh will be invisible
	**figure**
10. If necessary, rescale the **head mesh** to metres by right-clicking and selecting "Matrix: Set from translation/rotation/scale"
11. Apply a scale factor of 0.001 in x, y, and z. Select freeze matrix before applying the transform!
12. Press ctrl-H to reset the view
13. Crop the head mesh using the selector tool
14. Delete the selected points
15. Roughly align the meshes:
	1. Select the MRI scalp mesh
	2. Open the alignment tool (yellow A symbol)
	3. With the MRI scalp mesh selected, click Glue Here
	4. Now select the head mesh and click Point Based Gluing
	5. Select matching points on both meshes using double clicks and press OK
16. The meshes should be roughly aligned like so:
17. Freeze the head-mesh matrix.
18. Export the cropped and aligned head mesh using _File>Export Mesh As..._
    Here, we name the output sub-101_head_cropped.ply
19. Return to MATLAB and click OK
20. When prompted, select the cropped mesh.
21. Reload the MRI, if using the default, select the one with "padded" suffix in the filename
22. If you are prompted with
	```octave
	Completed FLIRT without errors
	```
	You have successfully obtained a warped MRI! (saved as  \<project dir>\/\<head mesh name>\_template_file_crg.nii.gz)
Check whether the results make sense by loading \*\_einscan_image.nii and ensuring it lines up with the generated MRI