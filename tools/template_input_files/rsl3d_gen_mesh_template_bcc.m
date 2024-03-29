function model
% Livelink connection and model initialization
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
model.component('comp1').geom('geom1').feature('dif1').selection('input').set({'blk1' 'sph18'});
model.component('comp1').geom('geom1').feature('dif1').selection('input2').set({'sph1' 'sph17' 'sph2' 'sph3' 'sph4' 'sph5' 'sph6' 'sph7' 'sph8'});
model.component('comp1').geom('geom1').create('co1', 'Compose');
model.component('comp1').geom('geom1').feature('co1').set('formula', 'dif1+dif1*(sph9+sph10+sph11+sph12+sph13+sph14+sph15+sph16)');

model.component('comp1').geom('geom1').feature('fin').set('repairtoltype', 'relative');
model.component('comp1').geom('geom1').create('ige1', 'IgnoreEdges');
model.component('comp1').geom('geom1').feature('ige1').selection('input').set('fin(1)', [45 46 47 48 56 57 58 61 62 63 64 67]);
model.component('comp1').geom('geom1').create('mcd1', 'MeshControlDomains');
model.component('comp1').geom('geom1').feature('mcd1').selection('input').set('ige1(1)', [1 3 4 5 6 7 8 9 10]);
model.component('comp1').geom('geom1').run;

% Set physics parameters and boundary conditions
model.component('comp1').physics.create('c', 'CoefficientFormPDE', 'geom1');
model.component('comp1').physics('c').prop('ShapeProperty').set('order', 1);
model.component('comp1').physics('c').feature('cfeq1').set('f', 0);
model.component('comp1').physics('c').feature('init1').set('u', 1);
model.component('comp1').physics('c').create('dir1', 'DirichletBoundary', 2);
model.component('comp1').physics('c').feature('dir1').selection.set([3 4 7 8 10 11 12 13 14]);
model.component('comp1').physics('c').create('pc1', 'PeriodicCondition', 2);
model.component('comp1').physics('c').feature('pc1').selection.set([2 9]);
model.component('comp1').physics('c').create('pc2', 'PeriodicCondition', 2);
model.component('comp1').physics('c').feature('pc2').selection.set([1 15]);
model.component('comp1').physics('c').create('pc3', 'PeriodicCondition', 2);
model.component('comp1').physics('c').feature('pc3').selection.set([5 6]);

% Create mesh and set the parameters of its size
model.component('comp1').mesh.create('mesh1');
model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri2', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri3', 'FreeTri');
model.component('comp1').mesh('mesh1').create('cpf1', 'CopyFace');
model.component('comp1').mesh('mesh1').create('cpf2', 'CopyFace');
model.component('comp1').mesh('mesh1').create('cpf3', 'CopyFace');
model.component('comp1').mesh('mesh1').create('ftri4', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri5', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri6', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri7', 'FreeTri');
model.component('comp1').mesh('mesh1').create('cpf4', 'CopyFace');
model.component('comp1').mesh('mesh1').create('cpf5', 'CopyFace');
model.component('comp1').mesh('mesh1').create('cpf6', 'CopyFace');
model.component('comp1').mesh('mesh1').create('cpf7', 'CopyFace');
model.component('comp1').mesh('mesh1').create('ftri8', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri9', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri10', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri11', 'FreeTri');
model.component('comp1').mesh('mesh1').create('cpf8', 'CopyFace');
model.component('comp1').mesh('mesh1').create('cpf9', 'CopyFace');
model.component('comp1').mesh('mesh1').create('cpf10', 'CopyFace');
model.component('comp1').mesh('mesh1').create('cpf11', 'CopyFace');
model.component('comp1').mesh('mesh1').create('ftri12', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri13', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri14', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri15', 'FreeTri');
model.component('comp1').mesh('mesh1').create('cpf12', 'CopyFace');
model.component('comp1').mesh('mesh1').create('cpf13', 'CopyFace');
model.component('comp1').mesh('mesh1').create('cpf14', 'CopyFace');
model.component('comp1').mesh('mesh1').create('cpf15', 'CopyFace');
model.component('comp1').mesh('mesh1').create('ftet1', 'FreeTet');
model.component('comp1').mesh('mesh1').create('ftet2', 'FreeTet');
model.component('comp1').mesh('mesh1').feature('ftri1').selection.set([16]);
model.component('comp1').mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri2').selection.set([32]);
model.component('comp1').mesh('mesh1').feature('ftri2').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri3').selection.set([24]);
model.component('comp1').mesh('mesh1').feature('ftri3').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri4').selection.set([19]);
model.component('comp1').mesh('mesh1').feature('ftri4').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri5').selection.set([17]);
model.component('comp1').mesh('mesh1').feature('ftri5').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri6').selection.set([1]);
model.component('comp1').mesh('mesh1').feature('ftri6').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri7').selection.set([18]);
model.component('comp1').mesh('mesh1').feature('ftri7').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri8').selection.set([21]);
model.component('comp1').mesh('mesh1').feature('ftri8').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri9').selection.set([23]);
model.component('comp1').mesh('mesh1').feature('ftri9').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri10').selection.set([22]);
model.component('comp1').mesh('mesh1').feature('ftri10').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri11').selection.set([2]);
model.component('comp1').mesh('mesh1').feature('ftri11').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri12').selection.set([30]);
model.component('comp1').mesh('mesh1').feature('ftri12').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri13').selection.set([29]);
model.component('comp1').mesh('mesh1').feature('ftri13').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri14').selection.set([6]);
model.component('comp1').mesh('mesh1').feature('ftri14').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri15').selection.set([31]);
model.component('comp1').mesh('mesh1').feature('ftri15').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftet1').selection.geom('geom1', 3);
model.component('comp1').mesh('mesh1').feature('ftet1').selection.set([1 3 4 5 6 7 8 9 10]);
model.component('comp1').mesh('mesh1').feature('ftet1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftet2').selection.geom('geom1', 3);
model.component('comp1').mesh('mesh1').feature('ftet2').selection.set([2]);
model.component('comp1').mesh('mesh1').feature('ftet2').create('size1', 'Size');

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
model.component('comp1').mesh('mesh1').feature('cpf2').selection('source').set([32]);
model.component('comp1').mesh('mesh1').feature('cpf2').selection('destination').set([20]);
model.component('comp1').mesh('mesh1').feature('cpf3').selection('source').set([24]);
model.component('comp1').mesh('mesh1').feature('cpf3').selection('destination').set([28]);
model.component('comp1').mesh('mesh1').feature('cpf1').selection('source').set([16]);
model.component('comp1').mesh('mesh1').feature('cpf1').selection('destination').set([36]);
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hmax', 'max_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hmin', 'min_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hcurve', 0.3);
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hcurveactive', true);
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hnarrow', 0.85);
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hnarrowactive', true);
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hgrad', 1.35);
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hgradactive', true);
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('hmax', 'max_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('hmin', 'min_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('hcurve', 0.3);
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('hcurveactive', true);
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('hnarrow', 0.85);
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('hnarrowactive', true);
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('hgrad', 1.35);
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('hgradactive', true);
model.component('comp1').mesh('mesh1').feature('ftri6').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri6').feature('size1').set('hmax', 'max_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri6').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri6').feature('size1').set('hmin', 'min_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri6').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftri6').feature('size1').set('hcurve', 0.3);
model.component('comp1').mesh('mesh1').feature('ftri6').feature('size1').set('hcurveactive', true);
model.component('comp1').mesh('mesh1').feature('ftri6').feature('size1').set('hnarrow', 0.85);
model.component('comp1').mesh('mesh1').feature('ftri6').feature('size1').set('hnarrowactive', true);
model.component('comp1').mesh('mesh1').feature('ftri6').feature('size1').set('hgrad', 1.35);
model.component('comp1').mesh('mesh1').feature('ftri6').feature('size1').set('hgradactive', true);
model.component('comp1').mesh('mesh1').feature('ftri7').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri7').feature('size1').set('hmax', 'max_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri7').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri7').feature('size1').set('hmin', 'min_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri7').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftri7').feature('size1').set('hcurve', 0.3);
model.component('comp1').mesh('mesh1').feature('ftri7').feature('size1').set('hcurveactive', true);
model.component('comp1').mesh('mesh1').feature('ftri7').feature('size1').set('hnarrow', 0.85);
model.component('comp1').mesh('mesh1').feature('ftri7').feature('size1').set('hnarrowactive', true);
model.component('comp1').mesh('mesh1').feature('ftri7').feature('size1').set('hgrad', 1.35);
model.component('comp1').mesh('mesh1').feature('ftri7').feature('size1').set('hgradactive', true);
model.component('comp1').mesh('mesh1').feature('cpf4').selection('source').set([19]);
model.component('comp1').mesh('mesh1').feature('cpf4').selection('destination').set([39]);
model.component('comp1').mesh('mesh1').feature('cpf5').selection('source').set([17]);
model.component('comp1').mesh('mesh1').feature('cpf5').selection('destination').set([37]);
model.component('comp1').mesh('mesh1').feature('cpf6').selection('source').set([1]);
model.component('comp1').mesh('mesh1').feature('cpf6').selection('destination').set([15]);
model.component('comp1').mesh('mesh1').feature('cpf7').selection('source').set([18]);
model.component('comp1').mesh('mesh1').feature('cpf7').selection('destination').set([38]);
model.component('comp1').mesh('mesh1').feature('ftri8').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri8').feature('size1').set('hmax', 'max_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri8').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri8').feature('size1').set('hmin', 'min_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri8').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftri8').feature('size1').set('hcurve', 0.3);
model.component('comp1').mesh('mesh1').feature('ftri8').feature('size1').set('hcurveactive', true);
model.component('comp1').mesh('mesh1').feature('ftri8').feature('size1').set('hnarrow', 0.85);
model.component('comp1').mesh('mesh1').feature('ftri8').feature('size1').set('hnarrowactive', true);
model.component('comp1').mesh('mesh1').feature('ftri8').feature('size1').set('hgrad', 1.35);
model.component('comp1').mesh('mesh1').feature('ftri8').feature('size1').set('hgradactive', true);
model.component('comp1').mesh('mesh1').feature('ftri9').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri9').feature('size1').set('hmax', 'max_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri9').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri9').feature('size1').set('hmin', 'min_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri9').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftri9').feature('size1').set('hcurve', 0.3);
model.component('comp1').mesh('mesh1').feature('ftri9').feature('size1').set('hcurveactive', true);
model.component('comp1').mesh('mesh1').feature('ftri9').feature('size1').set('hnarrow', 0.85);
model.component('comp1').mesh('mesh1').feature('ftri9').feature('size1').set('hnarrowactive', true);
model.component('comp1').mesh('mesh1').feature('ftri9').feature('size1').set('hgrad', 1.35);
model.component('comp1').mesh('mesh1').feature('ftri9').feature('size1').set('hgradactive', true);
model.component('comp1').mesh('mesh1').feature('ftri10').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri10').feature('size1').set('hmax', 'max_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri10').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri10').feature('size1').set('hmin', 'min_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri10').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftri10').feature('size1').set('hcurve', 0.3);
model.component('comp1').mesh('mesh1').feature('ftri10').feature('size1').set('hcurveactive', true);
model.component('comp1').mesh('mesh1').feature('ftri10').feature('size1').set('hnarrow', 0.85);
model.component('comp1').mesh('mesh1').feature('ftri10').feature('size1').set('hnarrowactive', true);
model.component('comp1').mesh('mesh1').feature('ftri10').feature('size1').set('hgrad', 1.35);
model.component('comp1').mesh('mesh1').feature('ftri10').feature('size1').set('hgradactive', true);
model.component('comp1').mesh('mesh1').feature('ftri11').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri11').feature('size1').set('hmax', 'max_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri11').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri11').feature('size1').set('hmin', 'min_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri11').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftri11').feature('size1').set('hcurve', 0.3);
model.component('comp1').mesh('mesh1').feature('ftri11').feature('size1').set('hcurveactive', true);
model.component('comp1').mesh('mesh1').feature('ftri11').feature('size1').set('hnarrow', 0.85);
model.component('comp1').mesh('mesh1').feature('ftri11').feature('size1').set('hnarrowactive', true);
model.component('comp1').mesh('mesh1').feature('ftri11').feature('size1').set('hgrad', 1.35);
model.component('comp1').mesh('mesh1').feature('ftri11').feature('size1').set('hgradactive', true);
model.component('comp1').mesh('mesh1').feature('cpf8').selection('source').set([21]);
model.component('comp1').mesh('mesh1').feature('cpf8').selection('destination').set([33]);
model.component('comp1').mesh('mesh1').feature('cpf9').selection('source').set([23]);
model.component('comp1').mesh('mesh1').feature('cpf9').selection('destination').set([35]);
model.component('comp1').mesh('mesh1').feature('cpf10').selection('source').set([22]);
model.component('comp1').mesh('mesh1').feature('cpf10').selection('destination').set([34]);
model.component('comp1').mesh('mesh1').feature('cpf11').selection('source').set([2]);
model.component('comp1').mesh('mesh1').feature('cpf11').selection('destination').set([9]);
model.component('comp1').mesh('mesh1').feature('ftri12').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri12').feature('size1').set('hmax', 'max_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri12').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri12').feature('size1').set('hmin', 'min_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri12').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftri12').feature('size1').set('hcurve', 0.3);
model.component('comp1').mesh('mesh1').feature('ftri12').feature('size1').set('hcurveactive', true);
model.component('comp1').mesh('mesh1').feature('ftri12').feature('size1').set('hnarrow', 0.85);
model.component('comp1').mesh('mesh1').feature('ftri12').feature('size1').set('hnarrowactive', true);
model.component('comp1').mesh('mesh1').feature('ftri12').feature('size1').set('hgrad', 1.35);
model.component('comp1').mesh('mesh1').feature('ftri12').feature('size1').set('hgradactive', true);
model.component('comp1').mesh('mesh1').feature('ftri13').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri13').feature('size1').set('hmax', 'max_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri13').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri13').feature('size1').set('hmin', 'min_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri13').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftri13').feature('size1').set('hcurve', 0.3);
model.component('comp1').mesh('mesh1').feature('ftri13').feature('size1').set('hcurveactive', true);
model.component('comp1').mesh('mesh1').feature('ftri13').feature('size1').set('hnarrow', 0.85);
model.component('comp1').mesh('mesh1').feature('ftri13').feature('size1').set('hnarrowactive', true);
model.component('comp1').mesh('mesh1').feature('ftri13').feature('size1').set('hgrad', 1.35);
model.component('comp1').mesh('mesh1').feature('ftri13').feature('size1').set('hgradactive', true);
model.component('comp1').mesh('mesh1').feature('ftri14').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri14').feature('size1').set('hmax', 'max_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri14').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri14').feature('size1').set('hmin', 1.4);
model.component('comp1').mesh('mesh1').feature('ftri14').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftri14').feature('size1').set('hcurve', 0.3);
model.component('comp1').mesh('mesh1').feature('ftri14').feature('size1').set('hcurveactive', true);
model.component('comp1').mesh('mesh1').feature('ftri14').feature('size1').set('hnarrow', 0.85);
model.component('comp1').mesh('mesh1').feature('ftri14').feature('size1').set('hnarrowactive', true);
model.component('comp1').mesh('mesh1').feature('ftri14').feature('size1').set('hgrad', 1.35);
model.component('comp1').mesh('mesh1').feature('ftri14').feature('size1').set('hgradactive', true);
model.component('comp1').mesh('mesh1').feature('ftri15').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri15').feature('size1').set('hmax', 'max_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri15').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri15').feature('size1').set('hmin', 'min_size_tri');
model.component('comp1').mesh('mesh1').feature('ftri15').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftri15').feature('size1').set('hcurve', 0.3);
model.component('comp1').mesh('mesh1').feature('ftri15').feature('size1').set('hcurveactive', true);
model.component('comp1').mesh('mesh1').feature('ftri15').feature('size1').set('hnarrow', 0.85);
model.component('comp1').mesh('mesh1').feature('ftri15').feature('size1').set('hnarrowactive', true);
model.component('comp1').mesh('mesh1').feature('ftri15').feature('size1').set('hgrad', 1.35);
model.component('comp1').mesh('mesh1').feature('ftri15').feature('size1').set('hgradactive', true);
model.component('comp1').mesh('mesh1').feature('cpf12').selection('source').set([30]);
model.component('comp1').mesh('mesh1').feature('cpf12').selection('destination').set([26]);
model.component('comp1').mesh('mesh1').feature('cpf13').selection('source').set([29]);
model.component('comp1').mesh('mesh1').feature('cpf13').selection('destination').set([25]);
model.component('comp1').mesh('mesh1').feature('cpf14').selection('source').set([6]);
model.component('comp1').mesh('mesh1').feature('cpf14').selection('destination').set([5]);
model.component('comp1').mesh('mesh1').feature('cpf15').selection('source').set([31]);
model.component('comp1').mesh('mesh1').feature('cpf15').selection('destination').set([27]);
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
