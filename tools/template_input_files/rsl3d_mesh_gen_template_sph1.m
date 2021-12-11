function rsl3d_mesh_gen
%livelink and model initialization
import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.component.create('comp1', false);
model.component('comp1').geom.create('geom1', 3);
model.component('comp1').defineLocalCoord(false);
%set geometry and mesh generation parameters
boxLx      = BOXLX;
boxLy      = BOXLY;
boxLz      = BOXLZ;
Rnp_eff    = RNP_EFF;
x_center_sph1 = X_CENTER_SPH1;
y_center_sph1 = Y_CENTER_SPH1;
z_center_sph1 = Z_CENTER_SPH1;
dRnp       = 5.0;
dh_min_in  = 0.8;  %1.0;
dh_max_in  = 1.0;  %1.2; %benchmark defined value was equal to 2.4
dh_min_out = dh_max_in;
dh_max_out = 20.0;
%build geometry
model.component('comp1').geom('geom1').repairTolType('relative');
model.component('comp1').geom('geom1').create('blk1', 'Block');
model.component('comp1').geom('geom1').feature('blk1').set('base', 'center');
model.component('comp1').geom('geom1').feature('blk1').set('size', [boxLx boxLy boxLz]);
model.component('comp1').geom('geom1').create('sph1', 'Sphere');
model.component('comp1').geom('geom1').feature('sph1').set('pos', [x_center_sph1 y_center_sph1 z_center_sph1]);
model.component('comp1').geom('geom1').feature('sph1').set('r', Rnp_eff);
model.component('comp1').geom('geom1').create('sph2', 'Sphere');
model.component('comp1').geom('geom1').feature('sph2').set('pos', [x_center_sph1 y_center_sph1 z_center_sph1]);
model.component('comp1').geom('geom1').feature('sph2').set('r', Rnp_eff+dRnp);
model.component('comp1').geom('geom1').create('dif2', 'Difference');
model.component('comp1').geom('geom1').feature('dif2').selection('input').set({'blk1' 'sph2'});
model.component('comp1').geom('geom1').feature('dif2').selection('input2').set({'sph1'});
%----------------------------INSERT GRAFTING POINTS HERE---------------------------------------------%
%----------------------------------------------------------------------------------------------------%
model.component('comp1').geom('geom1').feature('fin').set('repairtoltype', 'relative');
model.component('comp1').geom('geom1').create('ige1', 'IgnoreEdges');
model.component('comp1').geom('geom1').feature('ige1').selection('input').set('fin(1)', [9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32]);
model.component('comp1').geom('geom1').create('mcd1', 'MeshControlDomains');
model.component('comp1').geom('geom1').feature('mcd1').selection('input').set('ige1(1)', [2]);
model.component('comp1').geom('geom1').run;

model.component('comp1').selection.create('c_dst_pc1', 'Explicit');
model.component('comp1').selection('c_dst_pc1').geom('geom1', 2);
model.component('comp1').selection.create('c_dst_pc2', 'Explicit');
model.component('comp1').selection('c_dst_pc2').geom('geom1', 2);
%set boundary conditions
model.component('comp1').physics.create('c', 'CoefficientFormPDE', 'geom1');
model.component('comp1').physics('c').create('dir1', 'DirichletBoundary', 2);
model.component('comp1').physics('c').feature('dir1').selection.set([6]);
model.component('comp1').physics('c').create('pc1', 'PeriodicCondition', 2);
model.component('comp1').physics('c').feature('pc1').selection.set([1 7]);
model.component('comp1').physics('c').create('pc2', 'PeriodicCondition', 2);
model.component('comp1').physics('c').feature('pc2').selection.set([2 5]);
model.component('comp1').physics('c').create('pc3', 'PeriodicCondition', 2);
model.component('comp1').physics('c').feature('pc3').selection.set([3 4]);
%set physics and mesh size parameters
model.component('comp1').mesh.create('mesh1');
model.component('comp1').mesh('mesh1').create('ftet2', 'FreeTet');
model.component('comp1').mesh('mesh1').create('ftet1', 'FreeTet');
model.component('comp1').mesh('mesh1').feature('ftet2').selection.geom('geom1', 3);
model.component('comp1').mesh('mesh1').feature('ftet2').selection.set([2]);
model.component('comp1').mesh('mesh1').feature('ftet2').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftet1').selection.geom('geom1', 3);
model.component('comp1').mesh('mesh1').feature('ftet1').selection.set([1]);

model.component('comp1').view('view1').set('transparency', true);

model.component('comp1').physics('c').prop('ShapeProperty').set('order', 1);
model.component('comp1').physics('c').prop('Units').set('CustomSourceTermUnit', 1);
model.component('comp1').physics('c').feature('cfeq1').set('c', {'17.98^2' '0' '0' '0' '17.98^2' '0' '0' '0' '17.98^2'});
model.component('comp1').physics('c').feature('cfeq1').set('f', 0);
model.component('comp1').physics('c').feature('init1').set('u', 1);

model.component('comp1').mesh('mesh1').feature('size').set('hauto', 6);
model.component('comp1').mesh('mesh1').feature('size').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size').set('hmax', dh_max_out);
model.component('comp1').mesh('mesh1').feature('size').set('hmin', dh_min_out);
model.component('comp1').mesh('mesh1').feature('size').set('hcurve', 0.3);
model.component('comp1').mesh('mesh1').feature('size').set('hnarrow', 0.85);
model.component('comp1').mesh('mesh1').feature('size').set('hgrad', 1.35);
model.component('comp1').mesh('mesh1').feature('ftet2').set('method', 'dellegacy52');
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hauto', 3);
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hmax', dh_max_in);
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hmin', dh_min_in);
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftet1').set('method', 'dellegacy52');
model.component('comp1').mesh('mesh1').run;
%solve the PDE
model.study.create('std1');
model.study('std1').create('time', 'Transient');
model.study('std1').feature('time').set('activate', {'c' 'on'});

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature.remove('fcDef');

model.result.create('pg1', 'PlotGroup3D');
model.result('pg1').create('slc1', 'Slice');
model.result.export.create('mesh1', 'Mesh');

model.study('std1').feature('time').set('tlist', 'range(0,1,10)');
model.study('std1').feature('time').set('discretization', {'c' 'physics'});

model.sol('sol1').attach('std1');
model.sol('sol1').label('Solver 1');
model.sol('sol1').feature('v1').set('clist', {'range(0,1,10)' '0.01[s]'});
model.sol('sol1').feature('t1').set('tlist', 'range(0,1,10)');
model.sol('sol1').feature('t1').feature('dDef').set('mumpsrreorder', false);
model.sol('sol1').feature('t1').feature('dDef').set('reusereorder', false);
model.sol('sol1').runAll;
%export the generated mesh and the comsol model file for inspection
model.result('pg1').feature('slc1').set('resolution', 'normal');
model.result.export('mesh1').set('type3D', 'mphascii');
model.result.export('mesh1').set('filename', 'mesh.in.mphtxt');
model.result.export('mesh1').run;

DIR = './';
FF  = 'tmp_rsl3d_mesh_gen';

filename = strcat(FF,'.mph');
file     = strcat(DIR,filename);
model.save(file);
