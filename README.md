
# Instructions for the template warping procedure

> Tested only on Ubuntu 22.03

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
	![Fieldtrip MRI slice plot](https://github.com/nsrhodes/template_warping/blob/main/screenshots/set_coord_sys.png)
	And the following prompt:
	```octave
	The coordinate system is not specified.
	Do you want to change the anatomical labels for the axes [Y, n]? 
	```
5. In this example we enter 'Y', 'r', 'a', 's', 'n', for the following prompts as it is a RAS coordinate system.
6. Wait for the scalp extraction to complete
7. You will be prompted with instructions on how to perform a rough alignment in Meshlab
	![Meshlab instruction prompt](https://github.com/nsrhodes/template_warping/blob/main/screenshots/meshlab_prompts.png)
    
8. Open a new terminal (ctrl+alt+T) and start Meshlab (see above)
9. Load the generated template MRI scalp
   (saved as <path_to_template_mri_file>\_template_points.stl)
   and the head mesh (here sub-101_head.ply)
   ![Meshlab, Load MRI mesh and head scan](https://github.com/nsrhodes/template_warping/blob/main/screenshots/meshlab_select.png)
	If the head mesh is in mm, the MRI mesh will be invisible
	![Meshlab interface - only head mesh is visible](https://github.com/nsrhodes/template_warping/blob/main/screenshots/only_one_visible.png)
10. If necessary, rescale the **head mesh** to metre-units by right-clicking and selecting "Matrix: Set from translation/rotation/scale"
![Meshlab, Apply transform ](https://github.com/nsrhodes/template_warping/blob/main/screenshots/scale_mesh.png)
11. Apply a scale factor of 0.001 in x, y, and z. Select freeze matrix before applying the transform!
![Meshlab, Set scale to metres](https://github.com/nsrhodes/template_warping/blob/main/screenshots/scale_mesh2.png)
12. Press ctrl-H to reset the view.
![Meshlab, Both meshes appear and are correctly scaled](https://github.com/nsrhodes/template_warping/blob/main/screenshots/post_scale.png)
13. If the MRI-derived mesh appears dark as above, invert the faces orientation:
![Meshlab, Flip mesh orientations](https://github.com/nsrhodes/template_warping/blob/main/screenshots/invert_faces.png)
![Meshlab, Mesh face orientations are correct](https://github.com/nsrhodes/template_warping/blob/main/screenshots/post_flip.png)
14. Crop the head mesh using the selector tool (blue circle).
![Meshlab, Crop mesh](https://github.com/nsrhodes/template_warping/blob/main/screenshots/crop.png)
15. Delete the selected points (red circle above).
16. Roughly align the meshes:
	1. Select the MRI scalp mesh.
	2. Open the alignment tool (yellow A symbol).
    ![Open the alignment tool](https://github.com/nsrhodes/template_warping/blob/main/screenshots/align.png)
	3. With the MRI scalp mesh selected, click Glue Here.
    ![Glue here button](https://github.com/nsrhodes/template_warping/blob/main/screenshots/glue_here.png)
	4. Now select the head mesh and click Point Based Gluing.
    ![Points based gluing button](https://github.com/nsrhodes/template_warping/blob/main/screenshots/points_based_gluing.png)
	5. Select matching points on both meshes using double clicks and press OK.
    ![Point selection in alignment tool](https://github.com/nsrhodes/template_warping/blob/main/screenshots/points_selection.png)
17. The meshes should be roughly aligned like so:
    ![Rougly aligned meshes](https://github.com/nsrhodes/template_warping/blob/main/screenshots/nearly_aligned.png)
18. Freeze the head-mesh matrix.
    ![Freeze transform button](https://github.com/nsrhodes/template_warping/blob/main/screenshots/freeze_matrix.png)
19. Export the cropped and aligned head mesh using _File>Export Mesh As..._
    Here, we name the output sub-101_head_cropped.ply
20. Return to MATLAB and click OK
21. When prompted, select the cropped mesh.
22. Reload the MRI, if using the default, select the one with "padded" suffix in the filename
23. You will se a cropping tool, allowing you to remove the neck of the template MRI mesh to ensure that the aligned meshes are of a similar shape. Adjust the slope and intercept using the sliders and press confirm when you are done. You will see an outline of the template brain, the shape of the template head (turqoise, and the outline of the filled in head mesh (yellow). Adjust the red line until it covers similar features on the MRI shape to the neck line of the yellow head shape.)
![MRI cropping tool](https://github.com/nsrhodes/template_warping/blob/main/screenshots/crop_mri.png)
24. If you are prompted with
	```octave
	Completed FLIRT without errors
	```
	You have successfully obtained a warped MRI! (saved as  \<project dir>\/\<head mesh name>\_template_file_crg.nii.gz)
Check whether the results make sense by loading \*\_einscan_image.nii and ensuring it lines up with the generated MRI
![Checking the outputs in FSLEyes](https://github.com/nsrhodes/template_warping/blob/main/screenshots/result.png)