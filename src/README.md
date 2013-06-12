# Hermitian Modular Forms for Fields of Low Discriminant

This is the source code directory which was developed along with the Diplom thesis. The initial code up to the hand over of the Diplom thesis was written by Albert Zeyer in 2013.

Both the code and the text can be found online at <https://github.com/albertz/diplom-thesis-math>.

This code implements the algorithm and all neccessary functions to calculate the vectorspace of Fourier expansions of Hermitian modular forms. All the theory is described in the text.

The main function to calculate the space is `herm_modform_space()` in the file `algo.py`.

At the time of writing, we need Sage 5.9 (5.8 might work) and these Sage patches:
- http://trac.sagemath.org/sage_trac/ticket/14240
- http://trac.sagemath.org/sage_trac/ticket/14497

Download the patch files, go to `$SAGEDIR/devel/sage`, do:

    hg import <patch_i>
    sage -b

In this directory, call `./compile.sh` and then `sage`. For example, you can do:

    import algo
    algo.herm_modform_space(D=-3, HermWeight=6, B_cF=7)

