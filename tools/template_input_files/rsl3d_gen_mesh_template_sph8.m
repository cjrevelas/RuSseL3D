function rsl3d_mesh_gen
%liveling and model initialization
import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.component.create('comp1', false);
model.component('comp1').geom.create('geom1', 3);
model.component('comp1').defineLocalCoord(false);

%set geometry and mesh generation parameters
boxLx         = BOXLX;
boxLy         = BOXLY;
boxLz         = BOXLZ;
Rnp_eff       = RNP_EFF;
dRnp          = 5.0;
dh_min_in     = 0.8;
dh_max_in     = 1.0;
dh_min_out    = dh_max_in;
dh_max_out    = 20.0;
x_center_sph1 = X_CENTER_SPH1;
y_center_sph1 = Y_CENTER_SPH1;
z_center_sph1 = Z_CENTER_SPH1;

x_center_sph2 = X_CENTER_SPH2;
y_center_sph2 = Y_CENTER_SPH2;
z_center_sph2 = Z_CENTER_SPH2;

x_center_sph3 = X_CENTER_SPH3;
y_center_sph3 = Y_CENTER_SPH3;
z_center_sph3 = Z_CENTER_SPH3;

x_center_sph4 = X_CENTER_SPH4;
y_center_sph4 = Y_CENTER_SPH4;
z_center_sph4 = Z_CENTER_SPH4;

x_center_sph5 = X_CENTER_SPH5;
y_center_sph5 = Y_CENTER_SPH5;
z_center_sph5 = Z_CENTER_SPH5;

x_center_sph6 = X_CENTER_SPH6;
y_center_sph6 = Y_CENTER_SPH6;
z_center_sph6 = Z_CENTER_SPH6;

x_center_sph7 = X_CENTER_SPH7;
y_center_sph7 = Y_CENTER_SPH7;
z_center_sph7 = Z_CENTER_SPH7;

x_center_sph8 = X_CENTER_SPH8;
y_center_sph8 = Y_CENTER_SPH8;
z_center_sph8 = Z_CENTER_SPH8;
%build geometry
model.component('comp1').geom('geom1').repairTolType('relative');
model.component('comp1').geom('geom1').create('blk1', 'Block');
model.component('comp1').geom('geom1').feature('blk1').set('base', 'center');
model.component('comp1').geom('geom1').feature('blk1').set('size', [boxLx boxLy boxLz]);
%----------------------------INSERT NANOPARTICLES HERE-----------------------------------------------%
%----------------------------------------------------------------------------------------------------%
model.component('comp1').geom('geom1').create('sph1', 'Sphere');
model.component('comp1').geom('geom1').feature('sph1').set('pos', [x_center_sph1 y_center_sph1 z_center_sph1]);
model.component('comp1').geom('geom1').feature('sph1').set('r', Rnp_eff);
model.component('comp1').geom('geom1').create('sph2', 'Sphere');
model.component('comp1').geom('geom1').feature('sph2').set('pos', [x_center_sph2 y_center_sph2 z_center_sph2]);
model.component('comp1').geom('geom1').feature('sph2').set('r', Rnp_eff);
model.component('comp1').geom('geom1').create('sph3', 'Sphere');
model.component('comp1').geom('geom1').feature('sph3').set('pos', [x_center_sph3 y_center_sph3 z_center_sph3]);
model.component('comp1').geom('geom1').feature('sph3').set('r', Rnp_eff);
model.component('comp1').geom('geom1').create('sph4', 'Sphere');
model.component('comp1').geom('geom1').feature('sph4').set('pos', [x_center_sph4 y_center_sph4 z_center_sph4]);
model.component('comp1').geom('geom1').feature('sph4').set('r', Rnp_eff);
model.component('comp1').geom('geom1').create('sph5', 'Sphere');
model.component('comp1').geom('geom1').feature('sph5').set('pos', [x_center_sph5 y_center_sph5 z_center_sph5]);
model.component('comp1').geom('geom1').feature('sph5').set('r', Rnp_eff);
model.component('comp1').geom('geom1').create('sph6', 'Sphere');
model.component('comp1').geom('geom1').feature('sph6').set('pos', [x_center_sph6 y_center_sph6 z_center_sph6]);
model.component('comp1').geom('geom1').feature('sph6').set('r', Rnp_eff);
model.component('comp1').geom('geom1').create('sph7', 'Sphere');
model.component('comp1').geom('geom1').feature('sph7').set('pos', [x_center_sph7 y_center_sph7 z_center_sph7]);
model.component('comp1').geom('geom1').feature('sph7').set('r', Rnp_eff);
model.component('comp1').geom('geom1').create('sph8', 'Sphere');
model.component('comp1').geom('geom1').feature('sph8').set('pos', [x_center_sph8 y_center_sph8 z_center_sph8]);
model.component('comp1').geom('geom1').feature('sph8').set('r', Rnp_eff);
model.component('comp1').geom('geom1').create('sph9', 'Sphere');
model.component('comp1').geom('geom1').feature('sph9').set('pos', [x_center_sph1 y_center_sph1 z_center_sph1]);
model.component('comp1').geom('geom1').feature('sph9').set('r', Rnp_eff+dRnp);
model.component('comp1').geom('geom1').create('sph10', 'Sphere');
model.component('comp1').geom('geom1').feature('sph10').set('pos', [x_center_sph2 y_center_sph2 z_center_sph2]);
model.component('comp1').geom('geom1').feature('sph10').set('r', Rnp_eff+dRnp);
model.component('comp1').geom('geom1').create('sph11', 'Sphere');
model.component('comp1').geom('geom1').feature('sph11').set('pos', [x_center_sph3 y_center_sph3 z_center_sph3]);
model.component('comp1').geom('geom1').feature('sph11').set('r', Rnp_eff+dRnp);
model.component('comp1').geom('geom1').create('sph12', 'Sphere');
model.component('comp1').geom('geom1').feature('sph12').set('pos', [x_center_sph4 y_center_sph4 z_center_sph4]);
model.component('comp1').geom('geom1').feature('sph12').set('r', Rnp_eff+dRnp);
model.component('comp1').geom('geom1').create('sph13', 'Sphere');
model.component('comp1').geom('geom1').feature('sph13').set('pos', [x_center_sph5 y_center_sph5 z_center_sph5]);
model.component('comp1').geom('geom1').feature('sph13').set('r', Rnp_eff+dRnp);
model.component('comp1').geom('geom1').create('sph14', 'Sphere');
model.component('comp1').geom('geom1').feature('sph14').set('pos', [x_center_sph6 y_center_sph6 z_center_sph6]);
model.component('comp1').geom('geom1').feature('sph14').set('r', Rnp_eff+dRnp);
model.component('comp1').geom('geom1').create('sph15', 'Sphere');
model.component('comp1').geom('geom1').feature('sph15').set('pos', [x_center_sph7 y_center_sph7 z_center_sph7]);
model.component('comp1').geom('geom1').feature('sph15').set('r', Rnp_eff+dRnp);
model.component('comp1').geom('geom1').create('sph16', 'Sphere');
model.component('comp1').geom('geom1').feature('sph16').set('pos', [x_center_sph8 y_center_sph8 z_center_sph8]);
model.component('comp1').geom('geom1').feature('sph16').set('r', Rnp_eff+dRnp);
%----------------------------INSERT GRAFTING POINTS HERE---------------------------------------------%
%----------------------------------------------------------------------------------------------------%
model.component('comp1').geom('geom1').feature('fin').set('repairtoltype', 'relative');
model.component('comp1').geom('geom1').run;


DIR = './';
FF  = 'model';

filename = strcat(FF,'.mph');
file     = strcat(DIR,filename);
model.save(file);
