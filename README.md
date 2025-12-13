# Final Project

Xiaoyi Liu (Jason), Veronika Grytsai

## Jason edit (selected feature: Offline Bezier Surface Tesselation)
for the 40n points of new features, I chose "Offline Bezier Surface Tesselation" which is four stars. Compared to the previous commit version (from project 6), I added a new class called BezierRing in which I create and tesselate bezier surface. 

Initially, I wanted to just replace the cylinders I currently have for each individual edge in the graph (that is the tree we grow using the Runions et al algorithm). I added perturbations to the conttrol points for the bezier surface so that each 'tube' becomes wobbly and (the hope is) more organic looking, which is helped by the fact that the bezier surface is smooth and we can just turn up the resolution.

In practice, because the tree have so many small segments, the visual effect of the above approach is not very obvious.

Therefore, I decided to take maximal degree-2 paths in the tree, so that we could group together short segments that belong to one 'section' of the tree. In addition to correctly retrieving these, the additional technical challenge is to ensure that the frames at select points over the resulting polyline paths are computed correctly, so that when we make a bezier surface 'tube' over these paths, we don't get incorrect cases where there is bad twisting in the control point topology. For frames at intermediate nodes, I take an average of the two edges that goes in and out of it, which helped maintain the thickness of the tube around many sharp corners. Then, I tesselate the bezier surface these control points create, so that we can get a much smoother branch (compared to the gaps between connecting cylinders in the previous approach). I think another advantage of this is that now we have fewer individual mesh elements overall, making it easier to do texturing and potentially more computationally efficient depending on the settings used. 

additional features:
- the auxin ('sources' in the code) is now in a sphere so the tree can look more realistic (growing straight up a bit before branching out as you said)
- there's now a floor
- since we're mainly going with the cylinders, I used an additional algorithm to make it so that the radius shrinks the 'newer' a branch is.

### Additional visualization to the video
replacing cylinders with bezier tubes with perturbed control points: <img width="1832" height="1372" alt="bezier-tubes-instead-of-cylinders" src="https://github.com/user-attachments/assets/f055bc35-8753-4ca3-a0d5-a7b3a7dc5317" />

using tubes spanning maximal degree-2 paths over the tree's underlying graph (lower meshing resolution): <img width="1832" height="1372" alt="bezier-tubes-over-multiple-edges-low-res" src="https://github.com/user-attachments/assets/abf7136c-015a-4ace-a369-4b506d5e9436" />

using tubes spanning maximal degree-2 paths over the tree's underlying graph (higher meshing resolution): <img width="1832" height="1372" alt="bezier-tubes-over-multiple-edges-high-res" src="https://github.com/user-attachments/assets/44208e39-ace6-4cd5-b13f-e207dede750e" />

the tree if we put auxins in a cube region some distance above the initial tree nodes: <img width="1832" height="1372" alt="auxin-in-cube-above-the-starter-nodes" src="https://github.com/user-attachments/assets/af425a29-e31d-47c3-9e43-e07bbbb6a158" />

## Veronika edit (selected feature: Bump Mapping)
I used the open source bump map of an oak tree to give our tree a nice texture. The Cylinder class was changed to calculate and store UVs, tangent, and bitangent, that are being later sent to shaders files. Shaders use bump maps and calculated tangent/bitangent to calculate new normal for each pixel using height map derivative method. This new normal is later calculated during Phong lighting model implementation.

The bump mapping appeared to give a really discreet change, that it was barely visible. I had to up the scale to make it more distinct, however from far away cylinders that create the tree are quite small, so it's hard to notice the texture.

additional features:
- I added a visualization option, which starts a timer and grows a tree step by step, so that is can be nicely visualized for the video
- played with background (gradient on a video), and a snowfall feature to add aesthetical quality to the video

