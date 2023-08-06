function model
% Livelink and model initialization
import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.component.create('comp1', false);
model.component('comp1').geom.create('geom1', 3);
model.component('comp1').defineLocalCoord(false);

% Set geometry and mesh generation parameters
boxLx      = BOXLX;
boxLy      = BOXLY;
boxLz      = BOXLZ;
dh_min_in  = 0.8;
dh_max_in  = 1.0;
dh_min_out = dh_max_in;
dh_max_out = 20.0;

% Build geometry
model.component('comp1').geom('geom1').repairTolType('relative');
model.component('comp1').geom('geom1').create('blk1', 'Block');
model.component('comp1').geom('geom1').feature('blk1').set('base', 'center');
model.component('comp1').geom('geom1').feature('blk1').set('size', [boxLx boxLy boxLz]);
%----------------------------INSERT NANOPARTICLES HERE-----------------------------------------------%
%----------------------------------------------------------------------------------------------------%
%----------------------------INSERT GRAFTING POINTS HERE---------------------------------------------%
%----------------------------------------------------------------------------------------------------%
model.component('comp1').geom('geom1').feature('fin').set('repairtoltype', 'relative');
model.component('comp1').geom('geom1').run;

% Create comsol mph model file
DIR = './';
FF  = 'model';

filename = strcat(FF,'.mph');
file     = strcat(DIR,filename);
model.save(file);
