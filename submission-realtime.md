## Project 5: Realtime

The project handout can be found [here](https://cs1230.graphics/projects/realtime/1).

### Output Comparison

Run the program, open the specified `.json` file and follow the instructions to set the parameters.

If the output you are generating is an image, save the image with the specified file name using the "Save Image" button in the UI. Be sure to follow the instructions in the left column to set the file name. Once you save the images, they will appear in the table below.

If the output you are generating is a video, reference the [Screen Recording Guide](https://cs1230.graphics/docs/screen-recording) to capture a video of the specified interaction. Once you've recorded everything, navigate to this template file in Github, click edit, and either select or drag-and-drop each of your videos into the correct location. This will upload them to GitHub (but not include them in the repo) and automatically embed them into this Markdown file by providing a link. Make sure to double-check that they all show up properly in the preview.

We're **not** looking for your video to **exactly** match the expected output (i.e. draw the exact same shape). Just make sure to follow the instructions and verify that the recording fully demonstrates the functionality of that part of your project.

> [!NOTE]
> There may be issues with the image saving functionality, and the images may not save the exact same way as the expected outputs. This is acceptable, as long as your output demonstrates correct functionality.

|                                           File/Method To Produce Output                                            |                                                     Expected Output                                                     |                                                                     Your Output                                                                     |
| :----------------------------------------------------------------------------------------------------------------: | :---------------------------------------------------------------------------------------------------------------------: | :-------------------------------------------------------------------------------------------------------------------------------------------------: |
|                Input: `unit_cone.json`<br/>Output: `unit_cone.png`<br/>Parameters: (5, 5, 0.1, 100)                |      ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/unit_cone.png)      |            ![Place unit_cone.png in student_outputs/realtime/required folder](student_outputs/realtime/required/unit_cone.png)            |
|            Input: `unit_cone_cap.json`<br/>Output: `unit_cone_cap.png`<br/>Parameters: (5, 5, 0.1, 100)            |    ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/unit_cone_cap.png)    |        ![Place unit_cone_cap.png in student_outputs/realtime/required folder](student_outputs/realtime/required/unit_cone_cap.png)        |
|               Input: `unit_cube.json`<br/>Output: `unit_cube.png`<br/> Parameters: (5, 5, 0.1, 100)                |      ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/unit_cube.png)      |            ![Place unit_cube.png in student_outputs/realtime/required folder](student_outputs/realtime/required/unit_cube.png)            |
|            Input: `unit_cylinder.json`<br/>Output: `unit_cylinder.png`<br/>Parameters: (5, 5, 0.1, 100)            |    ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/unit_cylinder.png)    |        ![Place unit_cylinder.png in student_outputs/realtime/required folder](student_outputs/realtime/required/unit_cylinder.png)        |
|              Input: `unit_sphere.json`<br/>Output: `unit_sphere.png`<br/>Parameters: (5, 5, 0.1, 100)              |     ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/unit_sphere.png)     |          ![Place unit_sphere.png in student_outputs/realtime/required folder](student_outputs/realtime/required/unit_sphere.png)          |
|          Input: `unit_cone.json`<br/>Output: `unit_cone_min.png`<br/>Parameters: (**1**, **3**, 0.1, 100)          |      ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/unit_cone_min.png)      |        ![Place unit_cone_min.png in student_outputs/realtime/required folder](student_outputs/realtime/required/unit_cone_min.png)        |
|      Input: `unit_cone_cap.json`<br/>Output: `unit_cone_cap_min.png`<br/>Parameters: (**1**, **3**, 0.1, 100)      |    ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/unit_cone_cap_min.png)    |    ![Place unit_cone_cap_min.png in student_outputs/realtime/required folder](student_outputs/realtime/required/unit_cone_cap_min.png)    |
|          Input: `unit_cube.json`<br/>Output: `unit_cube_min.png`<br/>Parameters: (**1**, **1**, 0.1, 100)          |      ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/unit_cube_min.png)      |        ![Place unit_cube_min.png in student_outputs/realtime/required folder](student_outputs/realtime/required/unit_cube_min.png)        |
|      Input: `unit_cylinder.json`<br/>Output: `unit_cylinder_min.png`<br/>Parameters: (**1**, **3**, 0.1, 100)      |    ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/unit_cylinder_min.png)    |    ![Place unit_cylinder_min.png in student_outputs/realtime/required folder](student_outputs/realtime/required/unit_cylinder_min.png)    |
|        Input: `unit_sphere.json`<br/>Output: `unit_sphere_min.png`<br/>Parameters: (**2**, **3**, 0.1, 100)        |     ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/unit_sphere_min.png)     |      ![Place unit_sphere_min.png in student_outputs/realtime/required folder](student_outputs/realtime/required/unit_sphere_min.png)      |
|           Input: `parse_matrix.json`<br/>Output: `parse_matrix.png`<br/>Parameters: (**3**, 5, 0.1, 100)           |    ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/parse_matrix.png)     |         ![Place parse_matrix.png in student_outputs/realtime/required folder](student_outputs/realtime/required/parse_matrix.png)         |
|            Input: `ambient_total.json`<br/>Output: `ambient_total.png`<br/>Parameters: (5, 5, 0.1, 100)            |    ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/ambient_total.png)    |        ![Place ambient_total.png in student_outputs/realtime/required folder](student_outputs/realtime/required/ambient_total.png)        |
|            Input: `diffuse_total.json`<br/>Output: `diffuse_total.png`<br/>Parameters: (5, 5, 0.1, 100)            |    ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/diffuse_total.png)    |        ![Place diffuse_total.png in student_outputs/realtime/required folder](student_outputs/realtime/required/diffuse_total.png)        |
|           Input: `specular_total.json`<br/>Output: `specular_total.png`<br/>Parameters: (5, 5, 0.1, 100)           |   ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/specular_total.png)    |       ![Place specular_total.png in student_outputs/realtime/required folder](student_outputs/realtime/required/specular_total.png)       |
|              Input: `phong_total.json`<br/>Output: `phong_total.png`<br/>Parameters: (5, 5, 0.1, 100)              |     ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/phong_total.png)     |          ![Place phong_total.png in student_outputs/realtime/required folder](student_outputs/realtime/required/phong_total.png)          |
|      Input: `directional_light_1.json`<br/>Output: `directional_light_1.png`<br/>Parameters: (5, 5, 0.1, 100)      | ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/directional_light_1.png) |  ![Place directional_light_1.png in student_outputs/realtime/required folder](student_outputs/realtime/required/directional_light_1.png)  |
| Input: `directional_light_2.json`<br/>Output: `directional_light_2.png`<br/>Parameters: (**10**, **10**, 0.1, 100) | ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/directional_light_2.png) |  ![Place directional_light_2.png in student_outputs/realtime/required folder](student_outputs/realtime/required/directional_light_2.png)  |
|      Input: `phong_total.json`<br/>Output: `phong_total_near_far.png`<br/>Parameters: (5, 5, **9.5**, **12**)      | ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/phong_total_near_far.png) | ![Place phong_total_near_far.png in student_outputs/realtime/required folder](student_outputs/realtime/required/phong_total_near_far.png) |
|      Input: `directional_light_1.json`<br/>Output: `directional_light_1_near_far.png`<br/>Parameters: (**25**, **25**, **8**, **10**)      | ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/directional_light_1_near_far.png) | ![Place directional_light_1_near_far.png in student_outputs/realtime/required folder](student_outputs/realtime/required/directional_light_1_near_far.png) |
| Input: `point_light_1.json`<br/>Output: `point_light_1.png`<br/>Parameters: (5, 5, 0.1, 100) | ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/point_light_1.png) | ![Place point_light_1.png in student_outputs/realtime/required folder](student_outputs/realtime/required/point_light_1.png) |
| Input: `point_light_2.json`<br/>Output: `point_light_2.png`<br/>Parameters: (5, 5, 0.1, 100) | ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/point_light_2.png) | ![Place point_light_2.png in student_outputs/realtime/required folder](student_outputs/realtime/required/point_light_2.png) |
| Input: `spot_light_1.json`<br/>Output: `spot_light_1.png`<br/> Parameters: (5, 5, 0.1, 100)  | ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/spot_light_1.png)  |  ![Place spot_light_1.png in student_outputs/realtime/required folder](student_outputs/realtime/required/spot_light_1.png)  |
|  Input: `spot_light_2.json`<br/>Output: `spot_light_2.png`<br/>Parameters: (5, 5, 0.1, 100)  | ![](https://raw.githubusercontent.com/BrownCSCI1230/scenefiles/main/realtime/required_outputs/spot_light_2.png)  |  ![Place spot_light_2.png in student_outputs/realtime/required folder](student_outputs/realtime/required/spot_light_2.png)  |

#### Camera Translation

_Instructions:_ Load `chess.json`. For about 1 second each in this order, press:

- <kbd>W</kbd>, <kbd>A</kbd>, <kbd>S</kbd>, <kbd>D</kbd> to move in each direction by itself
- <kbd>W+A</kbd> to move diagonally forward and to the left
- <kbd>S+D</kbd> to move diagonally backward and to the right
- <kbd>Space</kbd> to move up
- <kbd>Cmd/Ctrl</kbd> to move down

##### Expected Output

https://github.com/BrownCSCI1230/projects_realtime_template/assets/45575415/710ff8b4-6db4-445b-811d-f6c838741e67

##### Your Output

<!---
Paste your output on top of this comment!
-->

https://github.com/user-attachments/assets/eaeaae48-3835-400c-8c25-612ded8657b9



#### Camera Rotation

_Instructions:_ Load `chess.json`. Take a look around!

##### Expected Output

https://github.com/BrownCSCI1230/projects_realtime_template/assets/45575415/a14f4d32-88ee-4f5f-9843-74dd5c89b9dd

##### Your Output

<!---
Paste your output on top of this comment!
-->

https://github.com/user-attachments/assets/bb8bcc2c-e7e0-4c76-a510-68e970387ad3

### Design Choices
In the code, the realtime file is managing what's going on in the OpenGL pipeline. I have scene parsing in the utils folder to read the geometry and shader settings. I also have separate files for the camera.

Separately, I have a shapes folder for all the different tesselations of the primitives, including the cylinder. I then have a ShapeMeshes files to port them over for the real time rendering in this project.

For shaders, I have the default.vert and .frag to implement the shading methods such as the Phong method.

In this way, the camera, the mesh generation, and the OpenGL codes are separate. This is helpful for if we want to add new shapes quickly.
### Collaboration/References

### Known Bugs
(I suspect that) due to the usage of implicit normals where applicable, such as on the curved surfaces of the sphere, cone, and cylinder, when the discretization resolution is very low, I note that for the cone, the (specular) highlight on a single triangular face has a 'turn.' Using the normal of the triangular faces themselves seem to remove the highlight entirely.

### Extra Credit
I chose the adaptive level of detail and completed both features (3 and 7 points respectively).

For the first part (adaptive based on number of objects), I use the total number of objects in the scene to make a global complexity factor. When this is small, the result is closer to whatever resolution we set. If this number gets bigger, the resolution decrease toward a minimum threshold. This applies to all objects uniformly, so the effect is much stronger when there are lots of objects in the scene.

For the second part (adaptive by distance), I compute the distance from the model space origin of each object (after it has been transformed and placed in the scene) to the position of the camera eye. I then map this to a factor from 1.0 (for being at the near plane) to a smaller number (at the far plane)

At the tesselation time, I multiple the two parameters by the factors derived from both parts. Of course, if only one part is enabled, it will not have an effect. I connected this to the first and second extra credit buttons so we can see what happens if either or both is on/off. For the test case, I chose a relatively crowded recursive spheres example because here, all the primitives are the same but spread out over a long distance, so we can more easily see the changes in resolution. I took the screenshots holding all the parameters the same, only changing the switches for these two extra credit features. When part 1 is on, we can see that the resolution of everything is very low. When part 2 is on, we can see that the resolution of spheres further away is smaller than what is closer to us. The caveat is that if the location of the far plane is much further than the distance spannd by the object, the differences in resolution can be less pronounced. In practice, if the far plane is not a user-defined parameter, we try have it based on objects that are furthest away, so that we do not have to see further than we have to. When both part 1 and part 2 are on, the resolution becomes very low, as we are multiplying their factors together -- the geometries now tend toward their lowest possible resolution as allowed by my implementation (i.e. lowest before we get degenerate shape with no volume)

The below images shows only the first part, only the second part, and then both parts enabled on the same mesh
<img width="1024" height="768" alt="recursive_sphere_4-ec1a" src="https://github.com/user-attachments/assets/58d0b3c5-6d34-415b-bfd4-ade0e6e4725b" />
<img width="1024" height="768" alt="recursive_sphere_4-ec1b" src="https://github.com/user-attachments/assets/1e8c5ff6-a90e-4fdc-ba3d-b47199502a87" />
<img width="1024" height="768" alt="recursive_sphere_4-ec1ab" src="https://github.com/user-attachments/assets/02762d71-d566-4f5c-8ac0-eaf853cae3ed" />



