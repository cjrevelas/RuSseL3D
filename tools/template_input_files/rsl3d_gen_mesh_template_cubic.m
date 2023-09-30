function model
% Liveling connection and model initialization
import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.component.create('comp1', false);
model.component('comp1').geom.create('geom1', 3);
model.component('comp1').defineLocalCoord(false);

% Set geometry and mesh generation parameters
boxLx = BOXLX;
boxLy = BOXLY;
boxLz = BOXLZ;

model.param.set('max_size_tri', '1.8');
model.param.set('min_size_tri', '1.4');
model.param.set('max_size_tet_in', '1.8');
model.param.set('min_size_tet_in', '1.4');
model.param.set('max_size_tet_out', '20');

% Build geometry
model.component('comp1').view('view1').set('transparency', true);

model.component('comp1').geom('geom1').repairTolType('relative');
model.component('comp1').geom('geom1').create('blk1', 'Block');
model.component('comp1').geom('geom1').feature('blk1').set('base', 'center');
model.component('comp1').geom('geom1').feature('blk1').set('size', [boxLx boxLy boxLz]);
%----------------------------INSERT NANOPARTICLES HERE-----------------------------------------------%
%----------------------------------------------------------------------------------------------------%
%----------------------------INSERT GRAFTING POINTS HERE---------------------------------------------%
%----------------------------------------------------------------------------------------------------%
model.component('comp1').geom('geom1').create('dif1', 'Difference');
model.component('comp1').geom('geom1').feature('dif1').selection('input').set({'blk1' 'sph2'});
model.component('comp1').geom('geom1').feature('dif1').selection('input2').set({'sph1'});
model.component('comp1').geom('geom1').feature('fin').set('repairtoltype', 'relative');
model.component('comp1').geom('geom1').create('ige1', 'IgnoreEdges');
model.component('comp1').geom('geom1').feature('ige1').selection('input').set('fin(1)', [9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32]);
model.component('comp1').geom('geom1').create('mcd1', 'MeshControlDomains');
model.component('comp1').geom('geom1').feature('mcd1').selection('input').set('ige1(1)', [2]);
model.component('comp1').geom('geom1').run;

% Set physics parameters and boundary conditions
model.component('comp1').physics.create('c', 'CoefficientFormPDE', 'geom1');
model.component('comp1').physics('c').prop('ShapeProperty').set('order', 1);
model.component('comp1').physics('c').feature('cfeq1').set('f', 0);
model.component('comp1').physics('c').feature('init1').set('u', 1);
model.component('comp1').physics('c').create('dir1', 'DirichletBoundary', 2);
model.component('comp1').physics('c').feature('dir1').selection.set([6]);
model.component('comp1').physics('c').create('pc1', 'PeriodicCondition', 2);
model.component('comp1').physics('c').feature('pc1').selection.set([3 4]);
model.component('comp1').physics('c').create('pc2', 'PeriodicCondition', 2);
model.component('comp1').physics('c').feature('pc2').selection.set([1 7]);
model.component('comp1').physics('c').create('pc3', 'PeriodicCondition', 2);
model.component('comp1').physics('c').feature('pc3').selection.set([2 5]);

% Create mesh and set the parameters of its size
model.component('comp1').mesh.create('mesh1');
model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri2', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri3', 'FreeTri');
model.component('comp1').mesh('mesh1').create('cpf1', 'CopyFace');
model.component('comp1').mesh('mesh1').create('cpf2', 'CopyFace');
model.component('comp1').mesh('mesh1').create('cpf3', 'CopyFace');
model.component('comp1').mesh('mesh1').create('ftet1', 'FreeTet');
model.component('comp1').mesh('mesh1').create('ftet2', 'FreeTet');
model.component('comp1').mesh('mesh1').feature('ftri1').selection.set([1]);
model.component('comp1').mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').selection.geom('geom1');
model.component('comp1').mesh('mesh1').feature('ftri2').selection.set([2]);
model.component('comp1').mesh('mesh1').feature('ftri2').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri3').selection.set([3]);
model.component('comp1').mesh('mesh1').feature('ftri3').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftet1').selection.geom('geom1', 3);
model.component('comp1').mesh('mesh1').feature('ftet1').selection.set([2]);
model.component('comp1').mesh('mesh1').feature('ftet1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftet2').selection.geom('geom1', 3);
model.component('comp1').mesh('mesh1').feature('ftet2').selection.set([1]);
model.component('comp1').mesh('mesh1').feature('ftet2').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').selection.set([2]);

model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hmax', 'max_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hmin', 'min_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hcurve', 0.3);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hcurveactive', true);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hnarrow', 0.85);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hnarrowactive', true);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hgrad', 1.35);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hgradactive', true);
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hmax', 'max_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hmin', 'min_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hcurve', 0.3);
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hcurveactive', true);
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hnarrow', 0.85);
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hnarrowactive', true);
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hgrad', 1.35);
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hgradactive', true);
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('hmax', 'max_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('hmin', 'min_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('hcurve', 0.3);
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('hcurveactive', true);
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('hnarrow', 0.85);
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('hnarrowactive', true);
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('hgrad', 1.35);
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('hgradactive', true);
model.component('comp1').mesh('mesh1').feature('cpf1').selection('source').set([1]);
model.component('comp1').mesh('mesh1').feature('cpf1').selection('destination').set([7]);
model.component('comp1').mesh('mesh1').feature('cpf2').selection('source').set([2]);
model.component('comp1').mesh('mesh1').feature('cpf2').selection('destination').set([5]);
model.component('comp1').mesh('mesh1').feature('cpf3').selection('source').set([3]);
model.component('comp1').mesh('mesh1').feature('cpf3').selection('destination').set([4]);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hmax', 'max_size_tet_in');
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hmin', 'min_size_tet_in');
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hcurve', 0.3);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hcurveactive', true);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hnarrow', 0.85);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hnarrowactive', true);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hgrad', 1.35);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hgradactive', true);
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hmax', 'max_size_tet_out');
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hmin', 'max_size_tet_in');
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hcurve', 0.3);
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hcurveactive', true);
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hnarrow', 0.85);
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hnarrowactive', true);
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hgrad', 1.35);
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hgradactive', true);
model.component('comp1').mesh('mesh1').run;

% Solve the model
model.study.create('std1');
model.study('std1').create('time', 'Transient');
model.study('std1').feature('time').set('tlist', 'range(0,0.2,1)');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature.remove('fcDef');
model.sol('sol1').feature('v1').set('clist', {'range(0,0.2,1)' '0.001[s]'});
model.sol('sol1').feature('t1').set('tlist', 'range(0,0.2,1)');
model.sol('sol1').runAll;

% Export mesh in ascii format
model.result.export.create('mesh1', 'Mesh');
model.result.export('mesh1').set('type3D', 'mphascii');
model.result.export('mesh1').set('filename', 'in.mesh.mphtxt');
model.result.export('mesh1').run;

% Save model in comsol mph format
DIR = './';
FF  = 'model';

filename = strcat(FF,'.mph');
file     = strcat(DIR,filename);
model.save(file);
