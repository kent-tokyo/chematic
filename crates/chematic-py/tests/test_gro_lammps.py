from pathlib import Path
import math
import pytest
import chematic

FIXTURE=Path(__file__).resolve().parents[3]/'validation/openbabel_issues/gro-lammps/orthogonal.gro'

def row(id=0, velocity=False):
    value=f'{7:5d}{"ION":<5}{"NA":>5}{id:5d}'+''.join(f'{v:8.3f}' for v in [-.25,.5,1.75])
    if velocity:value+=''.join(f'{v:8.4f}' for v in [-.125,.25,.5])
    return value

def frame(box='2 3 4', velocity=False):
    return 'original geometry\n2\n'+row(0,velocity)+'\n'+row(0,velocity)+'\n'+box+'\n'

def convert(text,units='real',types=(1,2),masses=(22.99,35.45)):
    return chematic.gro_to_lammps_atomic(text,types,masses,units=units)

@pytest.mark.parametrize('units,factor',[('real',.01),('metal',10.)])
def test_positions_box_velocities_and_sidecar(units,factor):
    text=frame(velocity=True)
    result=convert(text,units)
    data=chematic.parse_lammps_data(result['data_text'],'atomic')
    assert data['box']=={'lo':[0.,0.,0.],'hi':[20.,30.,40.],'tilt':None}
    assert [[a[k] for k in ['x','y','z']] for a in data['atoms']]==[[-2.5,5.,17.5]]*2
    assert [data['velocities'][0][k] for k in ['vx','vy','vz']]==[-.125*factor,.25*factor,.5*factor]
    assert [a['id'] for a in data['atoms']]==[1,2]
    assert [a['atom_type'] for a in data['atoms']]==[1,2]
    assert [m['mass'] for m in data['masses']]==[22.99,35.45]
    assert data['bonds']==[]
    assert result['source_text']==text
    assert result['units']==units and result['atom_style']=='atomic'
    assert result['length_scale']==10. and result['velocity_scale']==factor
    assert [a['atom_number'] for a in result['atom_identity']]==[0,0]
    assert [a['lammps_id'] for a in result['atom_identity']]==[1,2]
    assert result['atom_identity'][0]['atom_name']=='NA'
    assert result['atom_identity'][0]['residue_name']=='ION'
    assert result['atom_identity'][0]['residue_number']==7

def test_original_probe_cell_is_not_recomputed_from_atom_extent():
    text=FIXTURE.read_text(encoding='utf-8')
    result=convert(text,types=[1]*3,masses=[12.011])
    data=chematic.parse_lammps_data(result['data_text'],'atomic')
    assert data['box']['hi']==[20.,30.,40.]
    assert [a['x'] for a in data['atoms']]==[1.,4.,7.]

def test_triclinic_tilt_and_no_coordinate_wrapping():
    result=convert(frame('2 3 4 0 0 0.5 0 -0.25 0.75'))
    data=chematic.parse_lammps_data(result['data_text'],'atomic')
    assert data['box']['hi']==[20.,30.,40.]
    assert data['box']['tilt']==[5.,-2.5,7.5]
    assert data['atoms'][0]['x']==-2.5
    assert data['velocities']==[]
    assert all(a['molecule_id'] is None and a['charge'] is None and a['image'] is None for a in data['atoms'])

@pytest.mark.parametrize('box',['0 3 4','2 -3 4','2 3 0','2 3 4 0.1 0 0 0 0 0','2 3 4 0 0.1 0 0 0 0','2 3 4 0 0 0 0.1 0 0','1e308 3 4','nan 3 4'])
def test_invalid_or_unsupported_cells(box):
    with pytest.raises(ValueError):convert(frame(box))

@pytest.mark.parametrize('types',[[],[1],[0,1],[-1,1],[3,1],[2**63-1,1]])
def test_invalid_type_assignments(types):
    with pytest.raises(ValueError):convert(frame(),types=types)

@pytest.mark.parametrize('masses',[[],[0.],[-1.],[math.nan],[math.inf],[12.]*4097])
def test_invalid_mass_definitions(masses):
    with pytest.raises(ValueError):convert(frame(),masses=masses)

@pytest.mark.parametrize('units',['','lj','si','REAL'])
def test_unknown_units_refused(units):
    with pytest.raises(ValueError,match='units'):convert(frame(),units)

def test_units_required_and_wrong_types_refused():
    with pytest.raises(TypeError):chematic.gro_to_lammps_atomic(frame(),[1,2],[12.,16.])
    with pytest.raises(TypeError):chematic.gro_to_lammps_atomic(frame(),[1,2],[12.,16.],units=None)

def test_mixed_velocities_refused():
    text='mixed\n2\n'+row(velocity=True)+'\n'+row()+'\n2 3 4\n'
    with pytest.raises(ValueError,match='velocities'):convert(text)

def test_exact_crlf_source_and_multiframe_refusal():
    text=frame().replace('\n','\r\n')
    assert convert(text)['source_text']==text
    with pytest.raises(ValueError):convert(frame()+frame())

def test_input_size_limit_and_empty_frame():
    with pytest.raises(ValueError):convert('x'*((16<<20)+1),types=[])
    result=convert('empty\n0\n2 3 4\n',types=[],masses=[12.])
    assert result['atom_identity']==[]
