# SymPy Physics Calculator

The Sympy Physics Calculator is an ongoing project which, on completion, will solve graduate-level physics problems. Currently, it is two classes: A class for inertial reference frames, and a class for vectors with respect to a given inertial reference frame.

## Four_reference_frame

A Four_reference_frame object represents an inertial frame of reference in Special Relativity. An inertial frame of reference is one in which a physical object with no forces acting upon it moves with either constant or zero velocity. The reference frame can be likened to an all-powerful observer who can measure simultaneously all points in spacetime to determine where things are and how fast they move. Different observers who move at different (but constant) velocities relative to one another, are also inertial reference frames. A reference frame which undergoes acceleration is non-inertial.

In Special Relativity, a property that all inertial reference frames have is that the speed of light is the same. In other words, if I shine a flashlight at the wall and I measure the speed of light from the flashlight, I will find that number to be 3 x 10^8 m/s. If my friend travels alongside the beam of light at half the speed of light, he will measure the same light beam to be travelling at 3 x 10^8 m/s. In both of our frames of reference, we measure the same speed of light because our frames of reference are inertial.

### Basic Functionality

To initialize a new reference frame with the standard basis, call the constructor with no arguments:
```
R=Four_reference_frame()
```
The standard basis here refers to the set of basis vectors used to express vectors in this frame. The standard basis is `l=[1,0,0,0], i=[0,1,0,0], j=[0,0,1,0], k=[0,0,0,1]`. (In our convention, l,i,j,k are lists for the components of basis vectors in the t,x,y,z directions, respectivley.)Any vector expressed in this reference frame is a linear combination of these basis vectors. 

To create another reference frame, we can do so by either manually inputting the basis vectors or by transforming from this reference frame. The former is not recommended, as there is not yet checking for it the set of basis vectors is a valid set (i.e., is it a linearlly independent spanning set of the space). The latter is supported by the code. 

For example, if we wanted to create a reference frame travelling at velocity v, with `b=v/c` we would do:
```
L=R.general_boost(0,0,b)
```
This creates a new reference frame that is boosted in the z direction using the general_boost method. The first two arguments are the polar and azimuthal angles that describe the direction of the boost. Both being zero defines a direction pointing along the z-direction.

If we print the transformation matrix of L, it looks like this:
```
print(L.transformation_matrix)

Matrix([[sqrt(1/(1 - b**2)), 0, 0, -b*sqrt(1/(1 - b**2))],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [-b*sqrt(1/(1 - b**2)), 0, 0, sqrt(1/(1 - b**2))]])
```
This is a Lorentz boost in the z-direction, as expected. This transformation matrix represents the matrix that transforms vector components from the standard basis frame (R) to the new frame (L). 
If we create another frame Q, boosted from L by another velocity in another direction, Q.transformation_matrix will be the matrix composition of the transformations from R to L and from L to Q. 


## fourvector objects

Each reference frame in this framework has a transformation matrix formed out of it's basis vectors. The columns of the transformation matrix are the basis of the reference frame, expressed in terms of the standard basis. You can access these basis vectors directly using `R.l, R.i, R.j, R.k, etc.` 
Each basis vector in a reference frame is a fourvector object. To access their components you can say `R.i.components`, which will return a SymPy Matrix of shape (4,1) filled with the vector's components. 

To initialize a fourvector, you need a reference frame and a set of components to initialize it with: 
```
v=fourvector([components...],R)
w=fourvector([other components...],L)
```
Components can be either numeric types or SymPy symbols or expressions. They can be passed to the constructor either as a plain list or as a SymPy Matrix. The second argument is the reference frame to which this vector is associated. This can be accessed directly using `v.rframe`

Vectors can be added, subtracted, and multiplied by scalars. All of these return new fourvector objects with the appropriate components and reference frame. Furthermore, equality between vectors is defined by equality between the components of the vectors in the same reference frame.

Vectors can be expressed in different reference frames in the following way:
```
q=v.express(L)
q.rframe is L (returns True)
q==v (also True)
```
The express method takes a reference frame as an argument and returns a new fourvector object. This new object is a fourvector associated with the reference frame argument, whose components are the transformed version of v's components. In other words, they are v's components as expressed in reference frame L. 

Equality between reference frames is defined by the equality between the first vector's components, and the second vector's components transformed into the first's reference frame. This is because vectors expressed in different reference frames represent the same mathematical object. In other words, if I measure the position of a coconut to be at spacetime position `[t,x,y,z]` and my friend who is travelling at 0.5c in the z direction measures the position of the same coconut, the components of the two measurements will be different, but they represent the position of the same thing. 
In other words, the following two expressions are always equivalent in this framework:
```
q==v <==> q==v.express(q.rframe)
```

Addition and subtraction between different frames works in the same way. The second vector is transformed into the first's reference frame, the addition is performed, and the new fourvector is returned in the first argument's reference frame.

```
q+v==q+v.express(L) (always true)
x=q+v
x.rframe is q.rframe (always true)
```

The Lorentz dot product is computed using the dot method: `v.dot(v)` Dot products between vectors of different reference frames are not yet supported, due to the process adding complications to SymPy's symbolic manipulation.

Dot products are invariant under Lorentz transformations and so `v.dot(v)==q.dot(q)` because q is a Lorentz transformed version of v. 

## Future Developments

This is an ongoing project which I will continue to improve and refine over time. I plan to tackle scalar fields first, then vector fields. This will pave the way for being able to work with electromagnetic fields. Far in the future, I will write in code for covariant and contravariant indices, and relativistic tensors, which will pave the way for General Relativity. Stay tuned!


