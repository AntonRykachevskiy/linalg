# linalg
## NLA 2016 project

Recovering Low Ranked Data


Team - West Coast Low Ranked
Anton Rykachevsky
Ilya Gukov
Anna Voevodskaya
Polina Sannikova


Background
    As a team, our attention was picked by research and applications related to uncovering low-rank structures in data (for example, TILT). A lot of data around us actually has low rank structure - if one looks around, many things are covered in repetitive patterns. Static video feeds could be considered a low rank structure, if we consider time as a dimension. We hope to use tools of numerical linear algebra and optimization to recover these structures.


Problem Formulation
    Currently our interest lies in two applications:


3D Reconstruction of objects by recovering low-rank textures on 2D data (photos)
Recovering an image of a feature (e.g. a tourist attraction) without extra “noise” if form of passer-bies, transport etc. (Separating static background and moving foreground from the series of images, stacked into one matrix)


In any case, the problem could be boiled down to basically representing our initial matrix A as a sum of low-ranked matrix and sparse (but not obligatory small in the sense of element size) error matrix. 


Possible formulation looks like this


minimize rank(L) + cE
subject to A = L + E


Data
In case of 3D reconstruction, we could google several images of a building, or just take some photos ourselves. All we need are objects of simple form, with low ranked structure on them.
In case of image recovery, we might take some footage ourselves as well.


Related Work
http://perception.csl.illinois.edu/matrix-rank/Files/iccv11_3d_workshop.pdf
http://perception.csl.illinois.edu/matrix-rank/Files/TILT-ACCV10.pdf
https://www.robots.ox.ac.uk/~vgg/rg/papers/iccv11_curved.pdf


Scope
The end result is hopefully the program that does the task we give it to :) 


Evaluation
During our implementation, we will compare it to other known results - when they succeed, when they fail. The measure of success is obvious in both cases: does it recover what it is supposed to recover, or not. In our work we will definitely encounter cases when our method does not work, but we do not consider existence of such cases a failure.
