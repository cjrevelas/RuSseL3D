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

model.param.set('max_size_tri', '10.0');
model.param.set('min_size_tri', '0.5');
model.param.set('max_size_tet_in', '2.4');
model.param.set('min_size_tet_in', '0.5');
model.param.set('max_size_tet_out', '20.0');
model.param.set('min_size_tet_out', '0.5');

% Build geometry
model.component('comp1').view('view1').set('transparency',true);

model.component('comp1').geom('geom1').repairTolType('relative');
model.component('comp1').geom('geom1').create('blk1', 'Block');
model.component('comp1').geom('geom1').feature('blk1').set('base', 'center');
model.component('comp1').geom('geom1').feature('blk1').set('size', [boxLx boxLy boxLz]);
%----------------------------INSERT NANOPARTICLES HERE-----------------------------------------------%
%----------------------------------------------------------------------------------------------------%
%----------------------------INSERT GRAFTING POINTS HERE---------------------------------------------%
%----------------------------------------------------------------------------------------------------%
model.component('comp1').geom('geom1').run;



% Save model in comsol mph format
DIR = './';
FF  = 'model';

filename = strcat(FF,'.mph');
file     = strcat(DIR,filename);
model.save(file);
