"""
Functions for reduction of fourier indice of Hermitian modular forms.

AUTHOR :
    -- Dominic Gehre, Martin Raum (2009 - 08 - 20) Initial version
"""

# https://github.com/martinra/psage/blob/paper_computing_jacobi_forms/psage/modform/hermitianmodularforms/hermitianmodularformd2_fourierexpansion_cython.pyx

#===============================================================================
# 
# Copyright (C) 2009 Dominic Gehre, Martin Raum
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================


## to represent the semi-integral matrix [[a, b/2], [\bar b/2, c]] we choose the
## integer basis 1/\sqrt D, (1 + \sqrt D) / 2
class hermitian_form_with_character_evaluation :
    a = None
    b1 = None
    b2 = None
    c = None
    transposition_character = None
    determinant_character = None
    nu_character = None
    
#####################################################################
#####################################################################
#####################################################################

def reduce_GL(s, D) :
    """
    Reduce the positive semi-definite hermitian quatratic form `\lbrackt a,b,c \rbrackt`
    with respect to `\GL(2, \mathfrak{o}_{\Q(\sqrt{D})})` and the transposition.
    Return the evaluation of the characters `det` and `\nu` (see Dern where
    he calls this `\nu_\mathcal{p}`) as well as the number mod 2 of applied
    transpositions.
    
    INPUT:
        `s` -- A tuple `(a, b1, b2, c)`; We set `b = b1 / \sqrt{D} + b2 (1 + \sqrt{D})/2`.
        `D` -- An negative integer; The fundamental discriminant of the underlying
               imaginary quadratic number field.  
    
    OUTPUT:
        A pair `(s, (\mathrm{trans}, \mathrm{det}, \nu))` of a quadratic form `s`,
        which is represented by a four-tuple, and character evaluations, with
        `\nu` and `\mathrm{trans}` either `1` or `-1` and `\mathrm{det}` between
        `0` and `|\mathfrak{o}_{\Q(\sqrt{D})}^\times|`.
        
    TEST::
        sage: from hermitianmodularforms.hermitianmodularformd2_fourierexpansion_cython import reduce_GL
        sage: reduce_GL((2,1,1,2), -3)
        ((2, 1, 1, 2), (1, 0, 1))
        sage: reduce_GL((1,0,0,-1), -3)
        Traceback (most recent call last):
        ...
        NotImplementedError: only implemented for non-positive discriminants: (1, 0, 0, -1)
        sage: reduce_GL((-1,0,0,1), -3)
        Traceback (most recent call last):
        ...
        NotImplementedError: only implemented for non-positive discriminants: (-1, 0, 0, 1)
        sage: reduce_GL((-1,0,0,-1), -3)
        Traceback (most recent call last):
        ...
        NotImplementedError: only implemented for non-positive discriminants: (-1, 0, 0, -1)
        sage: reduce_GL((1,6,0,1), -3)
        Traceback (most recent call last):
        ...
        NotImplementedError: only implemented for non-positive discriminants: (1, 6, 0, 1)
        sage: reduce_GL((1,50,0,1), -7)
        Traceback (most recent call last):
        ...
        NotImplementedError: only implemented for non-positive discriminants: (1, 50, 0, 1)
    """
    (a, b1, b2, c) = s

    # the discriminant will be -D det(M) 
    # FIXME: the discriminant can become too big
    if a < 0 or c < 0 or - D * a * c - b1**2 - D * b1 * b2 - (D**2 - D) // 4 * b2**2 < 0 :
        raise NotImplementedError, "only implemented for non-positive discriminants: " + repr((a, b1, b2, c))

    res = _reduce_GL(a, b1, b2, c, D)

    return ( (res.a, res.b1, res.b2, res.c),
             (res.transposition_character, res.determinant_character, res.nu_character) )

def _reduce_GL(a, b1, b2, c, D) :
    """
    Reduce the positive semi-definite hermitian quatratic form
    `\lbrackt a,b,c \rbrackt` with respect to `\mathrm{GL}_2(\mathcal{o}_D)`.
    
    A reduced form `\lbrackt a,b,c \rbrackt` satisfies `a \le c`,
    `|\Im(b)| / \sqrt{D} \le a / 4`, `|\Re(b)| \le a / 2`. It satisfies `b2 \ge 0`. 
    If `D = -3` it also satisfies `-D b2 \ge 2 b1 \ge 2 b2`.
    If `D = -4` it also satisfies `2 b2 \ge b1 \ge b2`.
    In all other cases it also satisfies `2 b1 + D b2 \le 0`.
        
    TEST::
        sage: from hermitianmodularforms.hermitianmodularformd2_fourierexpansion_cython import reduce_GL
        sage: reduce_GL((1,0,0,2), -3) # indirect doctest
        ((1, 0, 0, 2), (1, 0, 1))
        sage: reduce_GL((2,0,0,1), -3) # indirect doctest
        ((1, 0, 0, 2), (1, 3, 1))
        sage: reduce_GL((2,1,2,2), -3) # indirect doctest
        ((1, 1, 1, 2), (-1, 1, 1))
        sage: reduce_GL((2,3,1,2), -3) # indirect doctest
        ((2, 3, 2, 2), (1, 1, 1))
        sage: reduce_GL((1,1,0,1), -3) # indirect doctest
        ((1, 1, 1, 1), (1, 2, 1))
        sage: reduce_GL((10,7, 8, 9), -3) # indirect doctest
        ((9, 10, 9, 10), (1, 0, 1))
        sage: reduce_GL((50, 76, 30, 19), -3) # indirect doctest
        ((19, 14, 11, 23), (1, 1, 1))
        sage: reduce_GL((1,0,0,2), -4) ## indirect doctest
        ((1, 0, 0, 2), (1, 0, 1))
        sage: reduce_GL((1,2,0,2), -4) ## indirect doctest
        ((1, 0, 0, 1), (1, 3, 1))
        sage: reduce_GL((1,3,1,1), -4) ## indirect doctest
        ((1, 1, 1, 1), (-1, 0, 1))
        sage: reduce_GL((17,0,5,10), -4) ## indirect doctest
        ((10, 15, 10, 17), (1, 3, 1))
        sage: reduce_GL((117,43,27,10), -4) ## indirect doctest
        ((10, 11, 9, 99), (-1, 1, 1))
        sage: reduce_GL((1,0,0,2), -7) ## indirect doctest
        ((1, 0, 0, 2), (1, 0, 1))
        sage: reduce_GL((1,2,0,2), -7) ## indirect doctest
        ((1, 2, 1, 2), (-1, 1, 1))
        sage: reduce_GL((17,0,5,10), -7) ## indirect doctest
        ((10, 0, 5, 17), (-1, 1, 1))
        sage: reduce_GL((117,43,27,10), -7) ## indirect doctest
        ((10, -6, 3, 65), (1, 0, 1))
    """
    res = hermitian_form_with_character_evaluation()
    twoa, q, r = [None]*3
    tmp = None
    ## for a describtion of det, trans and nu see the output documentation
    ## of reduce_GL
    det = 0
    trans = 1
    nu = 1

    
    ## K = QQ[\sqrt -1]
    if D == -4 :
        ## We want b to be of from the positive real axis by at most pi/4 in
        ## positive direction 
    
        ## if b is not in the right half plane
        if b2 < 0 :
            ## apply [[-1, 0], [0, 1]]
            b1 = -b1
            b2 = -b2
            
            det = det + 2
    
        ## if b is in the upper quarter of the right half plane
        if b1 < b2 :
            # apply [[i, 0], [0, 1]]
            tmp = -2 * b1 + 5 * b2
            b2 = -b1 + 2 * b2
            b1 = tmp
            
            det = det + 1
            
            ## apply transposition
            b1 = -b1 + 4 * b2
            
            trans = -trans

        ## if b is in the lower quarter of the right half plane
        elif b1 > 3 * b2 :
            ## apply [[-i, 0], [0, 1]]
            tmp = 2 * b1 - 5 * b2
            b2 = b1 - 2 * b2
            b1 = tmp
            
            det = det + 3
        
        ## b is in the second quarter counted from below of the right half plane
        elif b1 > 2 * b2 :
            ## apply transposition
            b1 = -b1 + 4 * b2
            
            trans = -trans

    ## K = QQ[\sqrt -3] 
    elif D == -3 :
        if b2 < 0 : 
            ## apply [[-1, 0], [0, 1]]
            b1 = -b1
            b2 = -b2
            
            det = det + 3

        ## TODO: Use one transformation for each sixth of the left half plane        
        ## check whether b is between 0 and pi/6 off from the real line in mathematical positive direction.
        ## b is in the upper quadrant
        if b1 < 0 :
            ## apply [[\bar eps, 0], [0, 1]]
            tmp = -b1 + 3*b2
            b2 = -b1 + 2*b2
            b1 = tmp
        
            det = det + 5

        ## b is in the last third of the lower quadrant 
        elif b1 > 3 * b2 :
            ## apply [[eps**2, 0], [0, 1]]
            tmp = -3 * b2 + b1
            b2 = b1 - 2 * b2
            b1 = tmp

            det = det + 2
        
        ## b is in the lower quadrant 
        elif 2 * b1 > 3 * b2 :
            ## apply [[eps, 0], [0, 1]]
            tmp = 2*b1 - 3*b2
            b2 = b1 - b2
            b1 = tmp
            
            det = det + 1
        
        ## now b is between 0 and pi/3 off from the real line in mathematical positive direction.
        ## check whether b is more than pi/6 off from the real in mathematical positive direction.
        if b1 < b2 :
            ## apply the transpose and
            ## apply [[eps, 0], [0, 1]]
            tmp = 3 * b2 - 2 * b1
            b2 = 2 * b2 - b1
            b1 = tmp

            trans = -trans
            det = det + 1
    #! elif D == -3
    else :
        # if b is not in the right half plane
        if b2 < 0 :
            ## apply [[-1, 0], [0, 1]]
            b1 = -b1
            b2 = -b2
            
            det = det + 1
        
        # if b is not in the upper half plane
        if 2 * b1 + D * b2 > 0 :
            b1 = - b1 - D * b2
            
            trans = -trans
            
    #! else D == -4

    ## TODO: Why is the _sig_on commented?
    #_sig_on

    ## while not GL-reduced
    while not (a<=c and abs(b2) <= a and abs(-4*b1 - 2 * D * b2) <= -D*a ) :
        ## apply symmetric step (reflection) if necessary
        if b1 < 0:
            ## apply [[-1, 0], [0, 1]]
            b1 = -b1
            b2 = -b2
            if D == -4 :
                det = det + 2
            elif D == -3 :
                det = det + 3
            else :
                det = det + 1
            
        ## abs(Im b / sqrt D) <= a/4 <=> -4*D abs(Im b) <= -D * a
        if abs(-4*b1 - 2 * D * b2) > -D*a :

            q = (-2*b1 - D * b2) / (-D * a)
            r = (-2*b1 - D * b2) % (-D * a)
            
            if r > (-D * a) / 2 :
                #r -= -D * a
                q += 1
            
            ## apply [[1, -q*(D + \sqrt D)/2], [0, 1]] 
            c = c + q**2 * a * (D**2 - D) / 4 + b1 * q
            b1 = b1 + q * a * (D-1) * D / 2
            b2 = b2 - q * a * D
    
        ## abs(Re b) <= a/2
        if abs (b2) > a :
            q = b2 / (2 * a)
            r = b2 % (2 * a)

            if r > a :
                r = r - 2 * a
                q = q + 1

            # apply [[1, -q], [0, 1]]
            c = c - b2 * q + a * q**2
            b1 = b1 + D * q * a 
            b2 = r # = b2 - 2 * q * a

        if a > c:
            ## apply [[0, 1], [1, 0]]
            tmp = c
            c = a
            a = tmp
            b1 = -b1 - D * b2
            if D == -4 :
                det = det + 2
            elif D == -3 :
                det = det + 3
            else :
                det = det + 1

        ## K = QQ[\sqrt -1]
        if D == -4 :
            ## We want b to be of from the positive real axis by at most pi/4 in
            ## positive direction 
        
            ## if b is not in the right half plane
            if b2 < 0 : 
                ## apply [[-1, 0], [0, 1]]
                b1 = -b1
                b2 = -b2
                
                det = det + 2
        
            ## if b is in the upper quarter of the right half plane
            if b1 < b2 :
                ## apply [[i, 0], [0, 1]]
                tmp = -2 * b1 + 5 * b2
                b2 = -b1 + 2 * b2
                b1 = tmp
                
                det = det + 1
                
                ## apply transposition
                b1 = -b1 + 4 * b2
                
                trans = -trans
                
            ## if b is in the lower quarter of the right half plane
            elif b1 > 3 * b2 :
                ## apply [[-i, 0], [0, 1]]
                tmp = 2 * b1 - 5 * b2
                b2 = b1 - 2 * b2
                b1 = tmp
                
                det = det + 3
            
            ## b is in the second quarter counted from below of the right half plane
            elif b1 > 2 * b2 :
                ## apply transposition
                b1 = -b1 + 4 * b2
                
                trans = -trans
            
        ## K = QQ[\sqrt -3] 
        elif D == -3 :
            if b2 < 0 : 
                ## apply [[-1, 0], [0, 1]]
                b1 = -b1
                b2 = -b2
                
                det = det + 3
    
            ## TODO: Use one transformation for each sixth of the left half plane        
            ## check whether b is between 0 and pi/6 off from the real line in mathematical positive direction.
            ## b is in the upper quadrant
            if b1 < 0 :
                ## apply [[\bar eps, 0], [0, 1]]
                tmp = -b1 + 3*b2
                b2 = -b1 + 2*b2
                b1 = tmp
                
                det = det + 5
                
            ## b is in the last third of the lower quadrant 
            elif b1 > 3 * b2 :
                ## apply [[eps**2, 0], [0, 1]]
                tmp = -3 * b2 + b1
                b2 = b1 - 2 * b2
                b1 = tmp
                
                det = det + 2
                
            ## b is in the lower quadrant 
            elif 2 * b1 > 3 * b2 :
                ## apply [[eps, 0], [0, 1]]
                tmp = 2*b1 - 3*b2
                b2 = b1 - b2
                b1 = tmp
                
                det = det + 1
            
            ## now b is between 0 and pi/3 off from the real line in mathematical positive direction.
            ## check whether b is more than pi/6 off from the real in mathematical positive direction.
            if b1 < b2 :
                ## apply the transpose and
                ## apply [[eps, 0], [0, 1]]
                tmp = 3 * b2 - 2 * b1
                b2 = 2 * b2 - b1
                b1 = tmp
                
                trans = -trans
                det = det + 1
        #! elif D == -3
        else :
            # if b is not in the right half plane
            if b2 < 0 :
                ## apply [[-1, 0], [0, 1]]
                b1 = -b1
                b2 = -b2
                
                det = det + 1
            
            # if b is not in the upper half plane
            if 2 * b1 + D * b2 > 0 :
                b1 = - b1 - D * b2
                
                trans = -trans
        #! else D == -4             
    #! while
        
    #_sig_off
    
    res.a = a
    res.b1 = b1
    res.b2 = b2
    res.c = c
    res.transposition_character = trans
    if D == -4 :
        res.determinant_character = det % 4
    elif D == -3 :
        res.determinant_character = det % 6
    else :
        res.determinant_character = det % 2
    res.nu_character = nu
            
    return res
