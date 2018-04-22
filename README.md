# sparse-synthextures
A technique for multispectral texture synthesis when sparse overcomplete dictionary of the texture is known. Proof of concept.

<center>
Input texture | Synthesized texture
:-------------------------------------------------------:|:---------------------------------------------------------:
![Tea leaves texture sample](/out/tealeaves_original.jpg)  |  ![Reproduced tea leaves texture](/out/tealeaves.jpg)
![Cherries texture sample](/out/cherries_original.jpg)  |  ![Reproduced cherries texture](/out/cherries.jpg)
![Desert sand texture sample](/out/desert_original.jpg)  |  ![Reproduced desert sand texture](/out/desert.jpg)
![Input fingerprint texture sample](/out/fingerprint_original.jpg)  |  ![Reproduced fingerprint texture](/out/fingerprint.jpg)
</center>

Image texture is assumed to be composed of overlapping image patches where each of the patches is modeled as a linear superposition of a small number of underlying dictionary elements, but with an additional requirement - at the overlaps, we aim to minimize the difference between pixel values of the overlapping patches. Hope is that, next to small-scale traits of the texture captured by the sparse dictionary (used for e.g. image inpainting), this constraint will bring out large-scale texture traits which are needed to synthesize an entirely new, visually convincing texture with minimally prononunced blocky artifacts.

This novel requirement is formalized as minimization of the following objective function,
\sum_{i,j} D[ **S1_{ij}****T****v_{i}**, **S2_{ij}****T****v_{j}**] ], where Euclidean distance is represented by __D__, _i_, _j_ are indeces of the patches composing the texture in a predefined order, **T** is dictionary with elements organised columnwise and **v_{i}** vector of sparse activation coefficients. For a pair of patches {_i_, _j_}, the overlapping region is described by two indicator matrices matrices **S1_{ij}** **S2_{ij}**. Also, note that sparsity penalization is used to promote sparsity of **v_{i}**.

The objective function is minimized by multiplicative updates under which nonnegativity of the activation coefficients is preserved.

If you have used this code in a publication, please mention this repository.

## License

Published under GPL-3.0 License.
