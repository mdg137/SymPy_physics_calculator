import pytest
from four_vector import *
@pytest.fixture
def ref_frame():
    return Four_reference_frame()
@pytest.fixture
def vector(ref_frame):
    x=sp.Symbol('x',real=True)
    y=sp.Symbol('y',real=True)
    z=sp.Symbol('z',real=True)
    t=sp.Symbol('t',real=True)
    return fourvector([t,x,y,z],ref_frame)

#@pytest.mark.skip(reason="I know it passes and takes too long for now")
def test_general_boost_theta(ref_frame):
    R=ref_frame
    b=sp.Symbol('b',real=True)
    assert R.general_boost(0,0,b).transformation_matrix==Four_reference_frame.lor_matrices['z'](b)
    assert R.general_boost(sp.pi/2,0,b).transformation_matrix==Four_reference_frame.lor_matrices['x'](b)
    assert R.general_boost(sp.pi,0,b).transformation_matrix==Four_reference_frame.lor_matrices['z'](-b)
    assert R.general_boost(3*sp.pi/2,0,b).transformation_matrix==Four_reference_frame.lor_matrices['x'](-b)
    assert R.general_boost(2*sp.pi,0,b).transformation_matrix==Four_reference_frame.lor_matrices['z'](b)
    assert R.general_boost(-sp.pi/2,0,b).transformation_matrix==Four_reference_frame.lor_matrices['x'](-b)
    assert R.general_boost(-sp.pi,0,b).transformation_matrix==Four_reference_frame.lor_matrices['z'](-b)
    assert R.general_boost(-3*sp.pi/2,0,b).transformation_matrix==Four_reference_frame.lor_matrices['x'](b)
    assert R.general_boost(-2*sp.pi,0,b).transformation_matrix==Four_reference_frame.lor_matrices['z'](b)
#@pytest.mark.skip(reason="I know it passes and takes too long for now")
def test_general_boost_phi(ref_frame):
    R=ref_frame
    b=sp.Symbol('b',real=True)
    phi=sp.Symbol('phi',real=True)
    assert R.general_boost(0,phi,b).transformation_matrix==Four_reference_frame.lor_matrices['z'](b)
    assert R.general_boost(sp.pi,phi,b).transformation_matrix==Four_reference_frame.lor_matrices['z'](-b)
    assert R.general_boost(-sp.pi,phi,b).transformation_matrix==Four_reference_frame.lor_matrices['z'](-b)
    assert R.general_boost(2*sp.pi,phi,b).transformation_matrix==Four_reference_frame.lor_matrices['z'](b)    
    assert R.general_boost(-2*sp.pi,phi,b).transformation_matrix==Four_reference_frame.lor_matrices['z'](b)  
#@pytest.mark.skip(reason="I know it passes and takes too long for now")
def test_both(ref_frame):
    R=ref_frame
    b=sp.Symbol('b',real=True)
    theta=sp.Symbol('theta',real=True)
    phi=sp.Symbol('phi',real=True)
    assert R.general_boost(theta,phi,b).transformation_matrix==R.general_boost(theta+2*sp.pi,phi,b).transformation_matrix
    assert R.general_boost(theta,phi,b).transformation_matrix==R.general_boost(theta-2*sp.pi,phi,b).transformation_matrix
    assert R.general_boost(theta,phi,b).transformation_matrix==R.general_boost(theta,phi+2*sp.pi,b).transformation_matrix
    assert R.general_boost(theta,phi,b).transformation_matrix==R.general_boost(theta,phi-2*sp.pi,b).transformation_matrix
    assert R.general_boost(theta,phi,-b).transformation_matrix==R.general_boost(theta+sp.pi,phi,b).transformation_matrix
    assert R.general_boost(theta,phi,-b).transformation_matrix==R.general_boost(theta-sp.pi,phi,b).transformation_matrix
    assert R.general_boost(theta,phi,b).transformation_matrix==R.general_boost(theta,phi,b).transformation_matrix.T

def test_four_vector_transformations_theta(ref_frame,vector):
    R=ref_frame
    b=sp.Symbol('b',real=True)
    b1=sp.Symbol('b1',real=True)
    theta=sp.Symbol('theta',real=True)
    v=vector
    assert type(v)==fourvector
    assert v.rframe is R 
    L=R.general_boost(0,0,b)
    q=v.express(L)
    assert q.rframe is L 
    assert q.components==sp.simplify(L.transformation_matrix*v.components)
    W=L.general_boost(theta,0,b1)
    p=v.express(W)
    assert p.rframe is W
    print(sp.factor(sp.trigsimp(p.components)))
    print(sp.trigsimp(W.transformation_matrix*v.components))
    assert p.components==sp.simplify(W.transformation_matrix*v.components)    
    q2=q.express(W)
    assert q2.rframe is W
    assert q2.components==sp.simplify(sp.simplify(W.transformation_matrix*(L.transformation_matrix**-1))*q.components)
    assert v.dot(v)==q.dot(q)
    assert p.dot(p)==v.dot(v)
    assert q.dot(q)==p.dot(p)
    assert q2.dot(q2)==q.dot(q)

def test_four_vector_transformations_phi(ref_frame,vector):
    R=ref_frame
    v=vector
    b=sp.Symbol('b',real=True)
    b1=sp.Symbol('b1',real=True)
    phi=sp.Symbol('phi',real=True)
    assert type(v)==fourvector
    assert v.rframe is R 
    L=R.general_boost(0,0,b)
    q=v.express(L)
    assert q.rframe is L 
    assert q.components==sp.simplify(L.transformation_matrix*v.components)
    W=L.general_boost(sp.pi/2,phi,b1)
    p=v.express(W)
    assert p.rframe is W
    assert p.components==sp.simplify(W.transformation_matrix*v.components)    
    q2=q.express(W)
    assert q2.rframe is W
    assert q2.components==sp.simplify(sp.simplify(W.transformation_matrix*(L.transformation_matrix**-1))*q.components)
    assert v.dot(v)==q.dot(q)
    assert p.dot(p)==v.dot(v)
    assert q.dot(q)==p.dot(p)
    assert q2.dot(q2)==q.dot(q)

def test_four_vector_operations(ref_frame,vector):
    R=ref_frame
    v=vector
    assert v+v==v*2
    assert v-v==v*0
    q=vector
    assert v==q
    c=sp.Symbol('c',real=True)
    assert c*v==v*c
    q=fourvector([0,1,2,3],R)
    assert v+q==q+v
    assert v-q==-(q-v)

    



