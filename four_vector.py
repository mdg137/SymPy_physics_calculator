import sympy as sp
from typing import Union
class Four_reference_frame:
    lorentz_metric=sp.Matrix([[1,0,0,0], 
                              [0,-1,0,0],
                              [0,0,-1,0],
                              [0,0,0,-1]]) 
    rot_matrices={
        'x':lambda a:sp.Matrix([[1,0,0,0],
                                [0,1,0,0],
                                [0,0,sp.cos(a),-sp.sin(a)],
                                [0,0,sp.sin(a),sp.cos(a)]]),
        'y':lambda a:sp.Matrix([[1,0,0,0],
                                [0,sp.cos(a),0,sp.sin(a)],
                                [0,0,1,0],
                                [0,-sp.sin(a),0,sp.cos(a)]]), 
        'z':lambda a:sp.Matrix([[1,0,0,0],
                                [0,sp.cos(a),-sp.sin(a),0],
                                [0,sp.sin(a),sp.cos(a),0],
                                [0,0,0,1]])
    }
    lor_matrices={ 
        'x':lambda b:sp.Matrix([[sp.sqrt(1/(1-b**2)),-b*sp.sqrt(1/(1-b**2)),0,0],
                                [-b*sp.sqrt(1/(1-b**2)),sp.sqrt(1/(1-b**2)),0,0],
                                [0,0,1,0],
                                [0,0,0,1]]),
        'y':lambda b:sp.Matrix([[sp.sqrt(1/(1-b**2)),0,-b*sp.sqrt(1/(1-b**2)),0],
                                [0,1,0,0],
                                [-b*sp.sqrt(1/(1-b**2)),0,sp.sqrt(1/(1-b**2)),0],
                                [0,0,0,1]]),
        'z':lambda b:sp.Matrix([[sp.sqrt(1/(1-b**2)),0,0,-b*sp.sqrt(1/(1-b**2))],
                                [0,1,0,0],
                                [0,0,1,0],
                                [-b*sp.sqrt(1/(1-b**2)),0,0,sp.sqrt(1/(1-b**2))]])
    }
    def __init__(self,l=sp.Matrix([1,0,0,0]),i=sp.Matrix([0,1,0,0]),j=sp.Matrix([0,0,1,0]),k=sp.Matrix([0,0,0,1])):
        """
        Initialize a new reference frame. Calling with no arguments initializes a reference frame with Standard Basis.

        Physically, this class represents inertial frames of reference. I am hopeful that non-inertial frames are possible,
        but right now that is beyond me. This code however, is intended to compute the symbolic mathematics involved in complex
        physics problems. The eventual conclusion will be the ability to solve Maxwell's equations using general boundary conditions
        and the presence of materials in spacetime.

        Right now, this class, together with class fourvector, function as a way to create vector spaces in Minkowski Spacetime,
        that is, 3+1D spacetime with no curvature, where the theory of Special Relativity resides. Because spacetime in 
        General Relativity has local lorentz invariance, this means that Special Relativity is relevant on small distance scales.
        Hence most of modern experimental and theoretical physics finds some basis in it. These classes allow you to create fourvectors 
        and look at their components in different reference frames. Examples shown below are currently time dilation, length contraction
        and relativistic velocity addition. 
        General functionality:

        R=Four_reference_frame() (calls __init__, assigns new reference frame to R.)
        print(R.i.components) (prints sp.Matrix([0,1,0,0]))

        print(R.transformation_matrix) (prints matrix which is used to transform to the standard basis. In this case it is the identity matrix)

        theta=sp.Symbol('theta',real=True)
        phi=sp.Symbol('phi',real=True)
        b=sp.Symbol('b',real=True)

        L=R.general_boost(theta,phi,b) (creates a new reference frame, boosted by b=v/c in a direction defined by declination angle (theta) and azimuthal angle (phi))
        print(L.transformation_matrix) (prints a lorentz transformation in a general direction defined by theta and phi. This matrix transforms vectors in R into vectors in L)

        v=fourvector([components...],R) (creates fourvector object associated with R reference frame)
        q=v.express(L) (creates a fourvector q which is the expression of v in reference frame L (in other words, vector v as observed by someone in reference frame L))
        """
        self.i=fourvector(i,self)
        self.j=fourvector(j,self)
        self.k=fourvector(k,self)
        self.l=fourvector(l,self)
        #transformation sp.Matrix that transforms components of vectors to reference frame from the previous frame
        self.transformation_matrix=l.row_join(i).row_join(j).row_join(k)        
    def general_boost(self,theta: Union[sp.Symbol,int,float], phi: Union[sp.Symbol,int,float], beta: Union[sp.Symbol,int,float] )->"Four_reference_frame":
        """A general Lorentz Boost works by rotating to orient the z-axis in the current frame with the orientation vector, then by boosting in the new z-direction with boost velocity beta=v/c,
        Then un-rotating the coordinate axes again. In Matrices: L_gen=R.T*L_z*R transforms the components of vectors in the current frame. The basis vectors transform with the inverse:
        L_gen**-1=R*L_z**-1*R.T. The function returns a reference frame with transformed basis vectors.

        Q=Four_reference_frame() (calls __init__, assigns new reference frame to R.)
        S=Q.general_boost(theta,phi,b) (theta,phi can be numbers or SymPy symbols, same for boost magnitude b=v/c, but note that b<1 always)
        """
        rotation_matrix=Four_reference_frame.rot_matrices['y'](-theta)*Four_reference_frame.rot_matrices['z'](-phi) 
        lorentz_matrix=Four_reference_frame.lor_matrices['z'](beta)
        transformation_matrix=sp.trigsimp((rotation_matrix.T*lorentz_matrix*rotation_matrix).T) # R*L*R.T
        return Four_reference_frame(transformation_matrix*self.l.components,transformation_matrix*self.i.components,transformation_matrix*self.j.components,transformation_matrix*self.k.components)
        #returns a reference frame with the transformed components as basis vectors
    def get_transformation(self,other: "Four_reference_frame")->sp.Matrix:
        """gets transformation between two reference frames (from self, to other), composed of 
           1) inverted self transformation into the standard basis and 
           2)transformation to other from standard basis

            This function is used to express vectors defined in one reference frame in another reference frame,
            inside fourvector.express(ref_frame) function
        """
        return sp.simplify(other.transformation_matrix*(self.transformation_matrix**-1))

class fourvector:
    def __init__(self,data: Union[list,sp.Matrix] ,ref_frame: Four_reference_frame)->None:
        """
        Fourvector class.

        A fourvector is a SymPy matrix containing its components, and a reference frame to which it 
        is associated. Fourvectors can be initialized using any existing reference frame in code.

        The fourvector.express function creates a new fourvector whose components are the components 
        of the current fourvector transformed into another reference frame. 

        Example:

        R=Four_reference_frame()
        L=R.general_boost(sp.Symbol('theta'),sp.Symbol('phi'),sp.Symbol('b'))
        v=fourvector([components],R)
        q=v.express(L)

        In this example, q and v are mathematically the same vector expressed in different
        reference frames. However, q and v are NOT the same object. In this code, q and v 
        are two DIFFERENT fourvector objects, which hold different components. Thus:

        q==v is True
        not(q is v) is also True

        The reason the first statement is still true is because equality between fourvectors of
        different reference frames is defined by transforming the second argument into the first's
        reference frame. So, regardless of whether q==v is true, the following is always true:

        (q==v)==(q==v.express(L))

        Addition, subtraction, and scalar multiplication all create new fourvector objects. 
        If addition or subtraction are done with two fourvectors in different reference frames, the convention
        for choosing a reference frame is to take the first argument's reference frame. So if we 
        want to compute q+v, the new fourvector object would be in q's reference frame. Similarly,
        v+q would be a fourvector in v's reference frame.

        When adding or subtracting fourvectors from different reference frames, the second argument is expressed
        in the first's reference frame. So, using the example above

        v+q is the same as v+q.express(R)
        q+v is the same as q+v.express(L)

        """
        self.components=sp.Matrix(data)
        self.rframe=ref_frame
    def __add__(self,other: 'fourvector')->'fourvector': #override addition, self-explanatory...
        if(self.rframe==other.rframe):
            return fourvector(self.components+other.components,self.rframe)
        else:
            return fourvector(self.components+other.express(self.rframe).components,self.rframe)            
    def __sub__(self,other: 'fourvector')->'fourvector':#...and subtraction. 
        if(self.rframe==other.rframe):
            return fourvector(self.components-other.components,self.rframe)
        else:
            return fourvector(self.components-other.express(self.rframe).components,self.rframe)            
    def __mul__(self,scalar: Union[int,float])->'fourvector': #scalar multiplication
        return fourvector(scalar*self.components,self.rframe)   
    def __rmul__(self,scalar: Union[int,float])->'fourvector': #scalar multiplication
        return fourvector(scalar*self.components,self.rframe)   
    def __eq__(self,other: 'fourvector')->'fourvector': #equality
        if(self.rframe==other.rframe):
            return self.components==other.components
        else:
            return self.components==other.express(self.rframe).components
    def __neg__(self)->'fourvector':
        return fourvector(-self.components,self.rframe)
    def dot(self,other: 'fourvector')->sp.Matrix:
        """
        Returns lorentz dot product of two fourvectors: 
        v.components.T*Four_reference_frame.lorentz_metric*q.components (a sp.Matrix object. To access the scalar directly, use v.dot(q)[0])

        if q and v are two fourvectors, then you would write:
            v.dot(q), or
            q.dot(v)
        
        The dot product of a vector with itself is preserved a lorentz transformation. In other words:
        Let R be the reference frame with the standard basis
        Let L=R.general_boost(theta,phi,b), I.E., a general lorentz boost with b=v/c
        Let v be a vector defined in R
        q=v.express(L), q is the same vector expressed in reference frame L's coodrinates
        q.dot(q)==v.dot(v) is always true

        dot products between fourvectors in different reference frames is not supported yet, as it 
        adds in complexity to SymPy's symbolic manipulation.
        """ 
        return sp.simplify(self.components.T*Four_reference_frame.lorentz_metric*other.components)
    def express(self,other_ref_frame: Four_reference_frame)->'fourvector':
        """Return a new fourvector that expresses self in another reference frame. 
        In other words, transform a vector as measured in one frame, to as measured in another.
        
        To express a fourvector in a different basis, the transformation is as follows:

        Let R be the reference frame with the standard basis and V and W be reference frames 
        that are lorentz transformations of R. These transformations could be rotations or boosts,
        and it is acceptable that either V or W could be identical to R. 

        In order to transform a vector expressed in space V to one expressed in space W, we transform using

        (WV**-1)v=w

        In other words, we transform back to the standard basis (V**-1) and then to frame W using W's transformation
        matrix.
           
        """
        transformation=self.rframe.get_transformation(other_ref_frame)
        out_vector_components=sp.simplify(transformation*self.components)
        return fourvector(out_vector_components,other_ref_frame)
    def get_four_velocity(self)->(sp.Matrix,sp.Matrix):
        """
        Rudimentary four velocity function. Not for use yet. 
        """
        improper_velocity=self.components.diff(t)
        boost=sp.sqrt((improper_velocity[1:4, 0].T*improper_velocity[1:4, 0])[0])
        gamma=sp.sqrt(1/(1-boost**2))
        proper_velocity=gamma*improper_velocity
        return proper_velocity,improper_velocity

        
