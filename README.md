# Final Project

Xiaoyi Liu (Jason), Veronika Grytsai

Jason edit: for the 40n points of new features, I chose "Offline Bezier Surface Tesselation" which is four stars. Compared to the previous commit version (from project 6), I added a new class called BezierRing in which I create and tesselate bezier surface. 

Initially, I wanted to just replace the cylinders I currently have for each individual edge in the graph (that is the tree we grow using the Runions et al algorithm). I added perturbations to the conttrol points for the bezier surface so that each 'tube' becomes wobbly and (the hope is) more organic looking, which is helped by the fact that the bezier surface is smooth and we can just turn up the resolution.

In practice, because the tree have so many small segments, the visual effect of the above approach is not very obvious.

Therefore, I decided to take maximal degree-2 paths in the tree, so that we could group together short segments that belong to one 'section' of the tree. In addition to correctly retrieving these, the additional technical challenge is to ensure that the frames at select points over the resulting polyline paths are computed correctly, so that when we make a bezier surface 'tube' over these paths, we don't get incorrect cases where there is bad twisting in the control point topology. For frames at intermediate nodes, I take an average of the two edges that goes in and out of it, which helped maintain the thickness of the tube around many sharp corners. Then, I tesselate the bezier surface these control points create, so that we can get a much smoother branch (compared to the gaps between connecting cylinders in the previous approach). I think another advantage of this is that now we have fewer individual mesh elements overall, making it easier to do texturing and potentially more computationally efficient depending on the settings used. 
