# sparse-synthextures
A method for multispectral texture synthesis using sparse image modeling. Proof of concept.

Input texture | Synthesized texture
:-------------------------------------------------------:|:---------------------------------------------------------:
![Tea leaves texture sample](/out/tealeaves_original.jpg)  |  ![Reproduced tea leaves texture](/out/tealeaves.jpg)
![Cherries texture sample](/out/cherries_original.jpg)  |  ![Reproduced cherries texture](/out/cherries.jpg)
![Desert sand texture sample](/out/desert_original.jpg)  |  ![Reproduced desert sand texture](/out/desert.jpg)
![Input fingerprint texture sample](/out/fingerprint_original.jpg)  |  ![Reproduced fingerprint texture](/out/fingerprint.jpg)

## Problem formulation
Here, similarly to inpainting applications, texture images are represented by overlapping image patches, where each patch is a sparse linear additive model. 
In addition to small-scale (single patch-level) traits of the texture that sparse modeling is able to capture by itself, we attempt to bring out large-scale texture traits in the synthesized image, ones which express how the patches connect to each other. More precisely, we propose to minimize the difference between pixel values of the overlapping patches with the aim to synthesize an entirely new, visually convincing texture with minimally prononunced blocky artifacts.

We formalize the described requirements as the following convex optimization problem:
min \sum_{i,j} D[ **S1_{ij}****T****v_{i}**, **S2_{ij}****T****v_{j}**] ], where Euclidean distance is represented by D[.], _i_, _j_ are indeces of the patches composing the texture in a predefined order, **T** is dictionary with elements organised columnwise and **v_{i}** vector of sparse activation coefficients. For a pair of patches {_i_, _j_}, the overlapping region is described by two indicator matrices matrices **S1_{ij}** **S2_{ij}**. Also, note that sparsity penalization is used to promote sparsity of **v_{i}**.

The objective function is minimized by multiplicative updates under which nonnegativity of the activation coefficients is preserved.

If you have used this code in a publication, please mention this repository.

## License

Published under GPL-3.0 License.