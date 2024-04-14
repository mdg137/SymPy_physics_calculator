from sympy import *
class fourvector():
    def __init__(self,data,ref_frame):
        if isinstance(data,Matrix):
            self.components=data
        else:
            self.components=Matrix(data)
        self.rframe=ref_frame
    def __add__(self,other): #override addition, self-explanatory...
        return fourvector(self.components+other.components,self.rframe)
    def __sub__(self,other):#...and subtraction. This is mostly an aesthetic choice though.
        return fourvector(self.components-other.components,self.rframe)
    def __mul__(self,scalar): #scalar multiplication
        return fourvector(scalar*self.components,self.rframe)   
    def __rmul__(self,scalar): #scalar multiplication
        return fourvector(scalar*self.components,self.rframe)   
    def __eq__(self,other): #equality
        return self.components==other.components
    def __neg__(self):
        return fourvector(-self.components,self.rframe)
    def dot(self,other): #dot product, employing the Lorentz Metric
        return simplify(self.components.T*Four_reference_frame.lorentz_metric*other.components)
    def express(self,other_ref_frame):
        #return a new fourvector that expresses self in another reference frame. In other words, transform a vector as measured in one frame, to as measured in another
        transformation=self.rframe.get_transformation(other_ref_frame)
        out_vector_components=simplify(transformation*self.components)
        return fourvector(out_vector_components,other_ref_frame)

class Four_reference_frame:
    lorentz_metric=Matrix([[1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,-1]]) 
    rot_matrices={
        'x':lambda a:Matrix([[1,0,0,0],[0,1,0,0],[0,0,cos(a),-sin(a)],[0,0,sin(a),cos(a)]]),
        'y':lambda a:Matrix([[1,0,0,0],[0,cos(a),0,sin(a)],[0,0,1,0],[0,-sin(a),0,cos(a)]]), 
        'z':lambda a:Matrix([[1,0,0,0],[0,cos(a),-sin(a),0],[0,sin(a),cos(a),0],[0,0,0,1]])
    }
    lor_matrices={
        'x':lambda b:Matrix([[sqrt(1/(1-b**2)),-b*sqrt(1/(1-b**2)),0,0],[-b*sqrt(1/(1-b**2)),sqrt(1/(1-b**2)),0,0],[0,0,1,0],[0,0,0,1]]),
        'y':lambda b:Matrix([[sqrt(1/(1-b**2)),0,-b*sqrt(1/(1-b**2)),0],[0,1,0,0],[-b*sqrt(1/(1-b**2)),0,sqrt(1/(1-b**2)),0],[0,0,0,1]]),
        'z':lambda b:Matrix([[sqrt(1/(1-b**2)),0,0,-b*sqrt(1/(1-b**2))],[0,1,0,0],[0,0,1,0],[-b*sqrt(1/(1-b**2)),0,0,sqrt(1/(1-b**2))]])
    }
    def __init__(self,l=Matrix([1,0,0,0]),i=Matrix([0,1,0,0]),j=Matrix([0,0,1,0]),k=Matrix([0,0,0,1])):
        self.i=fourvector(i,self)
        self.j=fourvector(j,self)
        self.k=fourvector(k,self)
        self.l=fourvector(l,self)
        #transformation matrix that transforms components of vectors to reference frame from the previous frame
        self.transformation_matrix=l.row_join(i).row_join(j).row_join(k)
        
    def general_boost(self,theta,phi,beta):
        #A general Lorentz Boost works by rotating to orient the z-axis in the current frame with the orientation vector, then by boosting in the new z-direction with boost velocity beta=v/c,
        #Then un-rotating the coordinate axes again. In Matrices: L_gen=R.T*L_z*R transforms the components of vectors in the current frame. The basis vectors transform with the inverse:
        # L_gen**-1=R*L_z**-1*R.T. The function returns a reference frame with transformed basis vectors.

        rotation_matrix=Four_reference_frame.rot_matrices['y'](-theta)*Four_reference_frame.rot_matrices['z'](-phi) 
        lorentz_matrix=Four_reference_frame.lor_matrices['z'](beta)
        transformation_matrix=trigsimp((rotation_matrix.T*lorentz_matrix*rotation_matrix).T) # R*L*R.T
        return Four_reference_frame(transformation_matrix*self.l.components,transformation_matrix*self.i.components,transformation_matrix*self.j.components,transformation_matrix*self.k.components)
        #returns a reference frame with the transformed components as basis vectors
    def get_transformation(self,other):
        #gets transformation between two reference frames, composed of 1) inverted self transformation into the standard basis and 2)transformation to other from standard basis
        return simplify(other.transformation_matrix*(self.transformation_matrix**-1))


R=Four_reference_frame()
