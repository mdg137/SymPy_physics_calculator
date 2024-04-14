import pytest
from four_vector import *
R=Four_reference_frame()
x=Symbol('x',real=True)
y=Symbol('y',real=True)
z=Symbol('z',real=True)
t=Symbol('t',real=True)
b=Symbol('b',real=True)
b1=Symbol('b1',real=True)
b2=Symbol('b2',real=True)
v=fourvector([t,x,y,z],R)
#@pytest.mark.skip(reason="I know it passes and takes too long for now")
def test_general_boost_theta():
    print(R.general_boost(0,0,b).transformation_matrix)
    print(Four_reference_frame.lor_matrices['z'](b))
    assert R.general_boost(0,0,b).transformation_matrix==Four_reference_frame.lor_matrices['z'](b)
    assert R.general_boost(pi/2,0,b).transformation_matrix==Four_reference_frame.lor_matrices['x'](b)
    assert R.general_boost(pi,0,b).transformation_matrix==Four_reference_frame.lor_matrices['z'](-b)
    assert R.general_boost(3*pi/2,0,b).transformation_matrix==Four_reference_frame.lor_matrices['x'](-b)
    assert R.general_boost(2*pi,0,b).transformation_matrix==Four_reference_frame.lor_matrices['z'](b)
    assert R.general_boost(-pi/2,0,b).transformation_matrix==Four_reference_frame.lor_matrices['x'](-b)
    assert R.general_boost(-pi,0,b).transformation_matrix==Four_reference_frame.lor_matrices['z'](-b)
    assert R.general_boost(-3*pi/2,0,b).transformation_matrix==Four_reference_frame.lor_matrices['x'](b)
    assert R.general_boost(-2*pi,0,b).transformation_matrix==Four_reference_frame.lor_matrices['z'](b)
#@pytest.mark.skip(reason="I know it passes and takes too long for now")
def test_general_boost_phi():
    phi=Symbol('phi',real=True)
    assert R.general_boost(0,phi,b).transformation_matrix==Four_reference_frame.lor_matrices['z'](b)
    assert R.general_boost(pi,phi,b).transformation_matrix==Four_reference_frame.lor_matrices['z'](-b)
    assert R.general_boost(-pi,phi,b).transformation_matrix==Four_reference_frame.lor_matrices['z'](-b)
    assert R.general_boost(2*pi,phi,b).transformation_matrix==Four_reference_frame.lor_matrices['z'](b)    
    assert R.general_boost(-2*pi,phi,b).transformation_matrix==Four_reference_frame.lor_matrices['z'](b)  
    
#@pytest.mark.skip(reason="I know it passes and takes too long for now")
def test_both():
    phi=Symbol('phi',real=True)
    theta=Symbol('theta',real=True)
    assert R.general_boost(theta,phi,b).transformation_matrix==R.general_boost(theta+2*pi,phi,b).transformation_matrix
    assert R.general_boost(theta,phi,b).transformation_matrix==R.general_boost(theta-2*pi,phi,b).transformation_matrix
    assert R.general_boost(theta,phi,b).transformation_matrix==R.general_boost(theta,phi+2*pi,b).transformation_matrix
    assert R.general_boost(theta,phi,b).transformation_matrix==R.general_boost(theta,phi-2*pi,b).transformation_matrix
    assert R.general_boost(theta,phi,-b).transformation_matrix==R.general_boost(theta+pi,phi,b).transformation_matrix
    assert R.general_boost(theta,phi,-b).transformation_matrix==R.general_boost(theta-pi,phi,b).transformation_matrix
    assert R.general_boost(theta,phi,b).transformation_matrix==R.general_boost(theta,phi,b).transformation_matrix.T


def test_four_vector_transformations_theta():
    theta=Symbol('theta',real=True)
    assert type(v)==fourvector
    assert v.rframe is R 
    L=R.general_boost(0,0,b)
    q=v.express(L)
    assert q.rframe is L 
    assert q.components==simplify(L.transformation_matrix*v.components)
    W=L.general_boost(theta,0,b1)
    p=v.express(W)
    assert p.rframe is W
    print(factor(trigsimp(p.components)))
    print(trigsimp(W.transformation_matrix*v.components))
    assert p.components==simplify(W.transformation_matrix*v.components)    
    q2=q.express(W)
    assert q2.rframe is W
    assert q2.components==simplify(simplify(W.transformation_matrix*(L.transformation_matrix**-1))*q.components)
    assert v.dot(v)==q.dot(q)
    assert p.dot(p)==v.dot(v)
    assert q.dot(q)==p.dot(p)
    assert q2.dot(q2)==q.dot(q)
def test_four_vector_transformations_phi():
    phi=Symbol('phi',real=True)
    assert type(v)==fourvector
    assert v.rframe is R 
    L=R.general_boost(0,0,b)
    q=v.express(L)
    assert q.rframe is L 
    assert q.components==simplify(L.transformation_matrix*v.components)
    W=L.general_boost(pi/2,phi,b1)
    p=v.express(W)
    assert p.rframe is W
    print(factor(trigsimp(p.components)))
    print(trigsimp(W.transformation_matrix*v.components))
    assert p.components==simplify(W.transformation_matrix*v.components)    
    q2=q.express(W)
    assert q2.rframe is W
    assert q2.components==simplify(simplify(W.transformation_matrix*(L.transformation_matrix**-1))*q.components)
    assert v.dot(v)==q.dot(q)
    assert p.dot(p)==v.dot(v)
    assert q.dot(q)==p.dot(p)
    assert q2.dot(q2)==q.dot(q)

def test_four_vector_operations():
    assert v+v==v*2
    assert v-v==v*0
    q=fourvector([t,x,y,z],R)
    assert v==q
    assert not(v is q)
    c=Symbol('c',real=True)
    assert c*v==v*c
    q=fourvector([0,1,2,3],R)
    assert v+q==q+v
    assert v-q==-(q-v)



