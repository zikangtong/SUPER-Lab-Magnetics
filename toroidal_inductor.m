%%  TOROIDAL INDUCTOR AUTOMATION SCRIPT IN COMSOL
%
%   Author: Zikang Tong
%   Affiliation: Stanford Universiy
%   Description: This script determines the inductance and ac loss of a
%   toroidal inductor in COMSOL. The input variables include the geometric
%   parameters such as radius and height as well as the number of turns for
%   each toroid. All units are in mm. 
%   
%   Variable Defintions: 
%   ri = inner radius
%   ro = outer radius
%   h = height
%   t = copper thickness
%   separation = gap separation between adjacent windings

mph_fileloc = 'G:\My Drive\Magnetics Project\COMSOL\';
mph_filename = 'delete.mph';

ri = 7.5; 
ro = 15;
h = 12.5;
t = 1.5;
separation = 1.5;
turns = 10;

I_prim = 1;
freq = 20e6;%units in Hz

%%%%%%%%%%%%%%%%%%%%%%%%%

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');
model.geom.create('geom1', 3);
model.geom('geom1').lengthUnit('mm');

%%%%%%%%%%%%%%%%%%%%%%%%%

model.geom('geom1').feature.create('wp1', 'WorkPlane');
model.geom('geom1').feature('wp1').set('planetype', 'quick');
model.geom('geom1').feature('wp1').set('quickplane', 'yz');

%   Create outer rectangle
model.geom('geom1').feature('wp1').geom.feature.create('outer_rec', 'Rectangle');
model.geom('geom1').feature('wp1').geom.feature('outer_rec').set('size', [ro-ri+t, h+t]);
model.geom('geom1').feature('wp1').geom.feature('outer_rec').set('pos', [ri-t/2, -h/2-t/2]);

%   Create inner rectangle
model.geom('geom1').feature('wp1').geom.feature.create('inner_rec', 'Rectangle');
model.geom('geom1').feature('wp1').geom.feature('inner_rec').set('size', [ro-ri-t, h-t]);
model.geom('geom1').feature('wp1').geom.feature('inner_rec').set('pos', [ri+t/2, -h/2+t/2]);

%   Creates difference of 2 unions
model.geom('geom1').feature('wp1').geom.feature.create('dif1', 'Difference');
model.geom('geom1').feature('wp1').geom.feature('dif1').selection('input').set({'outer_rec'});
model.geom('geom1').feature('wp1').geom.feature('dif1').selection('input2').set({'inner_rec'});

model.geom('geom1').feature.create('rev1', 'Revolve');
model.geom('geom1').feature('rev1').set('workplane', 'wp1');
model.geom('geom1').feature('rev1').selection('input').set({'wp1'});

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Creates the angled slicers
model.geom('geom1').feature.create('wp2', 'WorkPlane');
model.geom('geom1').feature('wp2').set('planetype', 'quick');
model.geom('geom1').feature('wp2').set('quickplane', 'xy');
model.geom('geom1').feature('wp2').set('quickz', -h/2-t/2);

model.geom('geom1').feature('wp2').geom.feature.create('r1', 'Rectangle');
model.geom('geom1').feature('wp2').geom.feature('r1').set('size', [sqrt((ro*cos(pi/turns)-ri)^2 + (ro*sin(pi/turns))^2), separation/2]);
model.geom('geom1').feature('wp2').geom.feature('r1').set('pos', [ri, 0]);
model.geom('geom1').feature('wp2').geom.feature('r1').set('rot', atan(ro*sin(pi/turns)/(ro*cos(pi/turns)-ri))*180/pi);
model.geom('geom1').feature('wp2').geom.feature.create('r2', 'Rectangle');
model.geom('geom1').feature('wp2').geom.feature('r2').set('size', [sqrt((ro*cos(pi/turns)-ri)^2 + (ro*sin(pi/turns))^2), separation/2]);
model.geom('geom1').feature('wp2').geom.feature('r2').set('pos', [ri+separation/2*cos(pi/2 - atan(ro*sin(pi/turns)/(ro*cos(pi/turns)-ri))), -1*separation/2*sin(pi/2 - atan(ro*sin(pi/turns)/(ro*cos(pi/turns)-ri)))]);
model.geom('geom1').feature('wp2').geom.feature('r2').set('rot', atan(ro*sin(pi/turns)/(ro*cos(pi/turns)-ri))*180/pi);
model.geom('geom1').feature('wp2').geom.create('uni1', 'Union');
model.geom('geom1').feature('wp2').geom.feature('uni1').selection('input').set({'r1' 'r2'});
model.geom('geom1').feature('wp2').geom.feature('uni1').set('intbnd', 'off');

for i = 1:turns-1
    model.geom('geom1').feature('wp2').geom.create(strcat('copy',int2str(i)), 'Copy');
    model.geom('geom1').feature('wp2').geom.feature(strcat('copy',int2str(i))).selection('input').set({'uni1'});
    model.geom('geom1').feature('wp2').geom.create(strcat('rot',int2str(i)), 'Rotate');
    model.geom('geom1').feature('wp2').geom.feature(strcat('rot',int2str(i))).selection('input').set(strcat('copy',int2str(i)));
    model.geom('geom1').feature('wp2').geom.feature(strcat('rot',int2str(i))).set('rot', 360/turns*i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Extrude
model.geom('geom1').feature.create('angled_slicer', 'Extrude');
model.geom('geom1').feature('angled_slicer').selection('input').set({'wp2'});
model.geom('geom1').feature('angled_slicer').setIndex('distance',h, 0);

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Rotate the angled slicers
model.geom('geom1').create('angled_slicer_copy', 'Copy');
model.geom('geom1').feature('angled_slicer_copy').selection('input').set('angled_slicer');
model.geom('geom1').create('angled_slicer_rot', 'Rotate');
model.geom('geom1').feature('angled_slicer_rot').selection('input').set('angled_slicer_copy');
model.geom('geom1').feature('angled_slicer_rot').set('axistype', 'x');
model.geom('geom1').feature('angled_slicer_rot').set('rot', 180);

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Create inner posts
post_list = {};
model.geom('geom1').create('post1', 'Block');
model.geom('geom1').feature('post1').set('size', [t*1.2, separation, h+t]);
model.geom('geom1').feature('post1').set('base', 'center');
model.geom('geom1').feature('post1').set('pos', [ri, 0, 0]);
post_list{end+1} = 'post1';
for i = 1:turns-1
    model.geom('geom1').create(strcat('post',int2str(i+1)), 'Copy');
    model.geom('geom1').feature(strcat('post',int2str(i+1))).selection('input').set('post1');
    model.geom('geom1').create(strcat('post_rot',int2str(i+1)), 'Rotate');
    model.geom('geom1').feature(strcat('post_rot',int2str(i+1))).selection('input').set(strcat('post',int2str(i+1)));
    model.geom('geom1').feature(strcat('post_rot',int2str(i+1))).set('rot', 360/turns*i);
    post_list{end+1} = strcat('post_rot',int2str(i+1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Create outter posts
model.geom('geom1').create('outpost1', 'Block');
model.geom('geom1').feature('outpost1').set('size', [t*1.2, separation, h+t]);
model.geom('geom1').feature('outpost1').set('base', 'center');
model.geom('geom1').feature('outpost1').set('pos', [ro, 0, 0]);
model.geom('geom1').create('outpost1_rot', 'Rotate');
model.geom('geom1').feature('outpost1_rot').selection('input').set('outpost1');
model.geom('geom1').feature('outpost1_rot').set('rot', 360/turns/2);
post_list{end+1} = 'outpost1_rot';

for i = 1:turns-1
    model.geom('geom1').create(strcat('outpost',int2str(i+1)), 'Copy');
    model.geom('geom1').feature(strcat('outpost',int2str(i+1))).selection('input').set('outpost1_rot');
    model.geom('geom1').create(strcat('outpost_rot',int2str(i+1)), 'Rotate');
    model.geom('geom1').feature(strcat('outpost_rot',int2str(i+1))).selection('input').set(strcat('outpost',int2str(i+1)));
    model.geom('geom1').feature(strcat('outpost_rot',int2str(i+1))).set('rot', 360/turns*i);
    post_list{end+1} = strcat('outpost_rot',int2str(i+1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Create openning
model.geom('geom1').feature.create('wp3', 'WorkPlane');
model.geom('geom1').feature('wp3').set('planetype', 'quick');
model.geom('geom1').feature('wp3').set('quickplane', 'xy');
model.geom('geom1').feature('wp3').set('quickz', -h/2+t/2);
model.geom('geom1').feature('wp3').geom.feature.create('c1', 'Circle');
model.geom('geom1').feature('wp3').geom.feature('c1').set('r', ro+t/2);
model.geom('geom1').feature('wp3').geom.feature('c1').set('angle', 360/turns);
model.geom('geom1').feature('wp3').geom.feature('c1').set('rot', 180/turns);
model.geom('geom1').feature('wp3').geom.feature.create('c2', 'Circle');
model.geom('geom1').feature('wp3').geom.feature('c2').set('r', ro-t/2);
model.geom('geom1').feature('wp3').geom.feature('c2').set('angle', 360/turns);
model.geom('geom1').feature('wp3').geom.feature('c2').set('rot', 180/turns);
model.geom('geom1').feature('wp3').geom.feature.create('dif1', 'Difference');
model.geom('geom1').feature('wp3').geom.feature('dif1').selection('input').set({'c1'});
model.geom('geom1').feature('wp3').geom.feature('dif1').selection('input2').set({'c2'});
model.geom('geom1').feature.create('openning', 'Extrude');
model.geom('geom1').feature('openning').selection('input').set({'wp3'});
model.geom('geom1').feature('openning').setIndex('distance', h-t, 0);

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Unions the slicer
model.geom('geom1').create('prim_slicer', 'Union');
post_list{end+1} = 'angled_slicer';
post_list{end+1} = 'angled_slicer_rot';
model.geom('geom1').feature('prim_slicer').selection('input').set(post_list);

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Takes difference
model.geom('geom1').create('dif1', 'Difference');
model.geom('geom1').feature('dif1').selection('input').set('rev1');
model.geom('geom1').feature('dif1').selection('input2').set({'prim_slicer', 'openning'});
model.geom('geom1').create('Final_rot', 'Rotate');
model.geom('geom1').feature('Final_rot').selection('input').set('dif1');
model.geom('geom1').feature('Final_rot').set('rot', -360/turns);

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Creates the lumped ports and the boundary
model.geom('geom1').create('blk3', 'Block');
model.geom('geom1').feature('blk3').set('base', 'center');
model.geom('geom1').feature('blk3').set('size', [1,1,h-t]);
model.geom('geom1').feature('blk3').set('pos', [ro,0,0]);
model.geom('geom1').create('sph1', 'Sphere');
model.geom('geom1').feature('sph1').set('r', ro*2);
model.geom('geom1').run;

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Creates the selections
%   Creating geometric definitions
model.selection.create('Winding', 'Explicit');
model.selection('Winding').set([2,3]);
model.selection.create('Non_conducting', 'Explicit');
model.selection('Non_conducting').set([1,4]);
model.selection.create('Magnetic_Field', 'Explicit');
model.selection('Magnetic_Field').set([1,4]);
model.selection.create('Conducting_Boundaries', 'Explicit');
model.selection('Conducting_Boundaries').geom('geom1', 3, 2, {'exterior'});
model.selection('Conducting_Boundaries').set([2,3]);
model.selection.create('Port1', 'Explicit');
model.selection('Port1').geom('geom1', 3, 2, {'exterior'});
model.selection('Port1').set(4);

% %%%%%%%%%%%%%%%%%%%%%%%%%
% %   Creates the materials
model.material.create('mat1', 'Common', 'mod1');
model.material('mat1').label('Air');
model.material('mat1').set('family', 'air');
model.material('mat1').propertyGroup('def').set('relpermeability', '1');
model.material('mat1').propertyGroup('def').set('relpermittivity', '1');
model.material('mat1').propertyGroup('def').set('dynamicviscosity', 'eta(T[1/K])[Pa*s]');
model.material('mat1').propertyGroup('def').set('ratioofspecificheat', '1.4');
model.material('mat1').propertyGroup('def').set('electricconductivity', '0[S/m]');
model.material('mat1').propertyGroup('def').set('heatcapacity', 'Cp(T[1/K])[J/(kg*K)]');
model.material('mat1').propertyGroup('def').set('density', 'rho(pA[1/Pa],T[1/K])[kg/m^3]');
model.material('mat1').propertyGroup('def').set('thermalconductivity', 'k(T[1/K])[W/(m*K)]');
model.material('mat1').propertyGroup('def').set('soundspeed', 'cs(T[1/K])[m/s]');
model.material('mat1').propertyGroup('def').func.create('eta', 'Piecewise');
model.material('mat1').propertyGroup('def').func('eta').set('funcname', 'eta');
model.material('mat1').propertyGroup('def').func('eta').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('eta').set('extrap', 'constant');
model.material('mat1').propertyGroup('def').func('eta').set('pieces', {'200.0' '1600.0' '-8.38278E-7+8.35717342E-8*T^1-7.69429583E-11*T^2+4.6437266E-14*T^3-1.06585607E-17*T^4'});
model.material('mat1').propertyGroup('def').func.create('Cp', 'Piecewise');
model.material('mat1').propertyGroup('def').func('Cp').set('funcname', 'Cp');
model.material('mat1').propertyGroup('def').func('Cp').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('Cp').set('extrap', 'constant');
model.material('mat1').propertyGroup('def').func('Cp').set('pieces', {'200.0' '1600.0' '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
model.material('mat1').propertyGroup('def').func.create('rho', 'Analytic');
model.material('mat1').propertyGroup('def').func('rho').set('funcname', 'rho');
model.material('mat1').propertyGroup('def').func('rho').set('args', {'pA' 'T'});
model.material('mat1').propertyGroup('def').func('rho').set('expr', 'pA*0.02897/8.314/T');
model.material('mat1').propertyGroup('def').func('rho').set('dermethod', 'manual');
model.material('mat1').propertyGroup('def').func('rho').set('argders', {'pA' 'd(pA*0.02897/8.314/T,pA)'; 'T' 'd(pA*0.02897/8.314/T,T)'});
model.material('mat1').propertyGroup('def').func.create('k', 'Piecewise');
model.material('mat1').propertyGroup('def').func('k').set('funcname', 'k');
model.material('mat1').propertyGroup('def').func('k').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('k').set('extrap', 'constant');
model.material('mat1').propertyGroup('def').func('k').set('pieces', {'200.0' '1600.0' '-0.00227583562+1.15480022E-4*T^1-7.90252856E-8*T^2+4.11702505E-11*T^3-7.43864331E-15*T^4'});
model.material('mat1').propertyGroup('def').func.create('cs', 'Analytic');
model.material('mat1').propertyGroup('def').func('cs').set('funcname', 'cs');
model.material('mat1').propertyGroup('def').func('cs').set('args', {'T'});
model.material('mat1').propertyGroup('def').func('cs').set('expr', 'sqrt(1.4*287*T)');
model.material('mat1').propertyGroup('def').func('cs').set('dermethod', 'manual');
model.material('mat1').propertyGroup('def').func('cs').set('argders', {'T' 'd(sqrt(1.4*287*T),T)'});
model.material('mat1').propertyGroup('def').addInput('temperature');
model.material('mat1').propertyGroup('def').addInput('pressure');
model.material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.material('mat1').propertyGroup('RefractiveIndex').set('n', '1');
model.material('mat1').set('family', 'air');
model.material.create('mat2', 'Common', 'mod1');
model.material('mat2').label('Copper');
model.material('mat2').set('family', 'copper');
model.material('mat2').propertyGroup('def').set('relpermeability', '1');
model.material('mat2').propertyGroup('def').set('electricconductivity', '5.998e7[S/m]');
model.material('mat2').propertyGroup('def').set('heatcapacity', '385[J/(kg*K)]');
model.material('mat2').propertyGroup('def').set('relpermittivity', '1');
model.material('mat2').propertyGroup('def').set('emissivity', '0.5');
model.material('mat2').propertyGroup('def').set('density', '8700[kg/m^3]');
model.material('mat2').propertyGroup('def').set('thermalconductivity', '400[W/(m*K)]');
model.material('mat2').propertyGroup.create('linzRes', 'Linearized resistivity');
model.material('mat2').propertyGroup('linzRes').set('alpha', '3.9e-3[1/K]');
model.material('mat2').propertyGroup('linzRes').set('rho0', '1.72e-8[ohm*m]');
model.material('mat2').propertyGroup('linzRes').set('Tref', '273.15[K]');
model.material('mat2').set('family', 'copper');
model.material.create('mat3', 'Common', 'mod1');
model.material('mat3').label('Copper 1');
model.material('mat3').set('family', 'copper');
model.material('mat3').propertyGroup('def').set('relpermeability', '1');
model.material('mat3').propertyGroup('def').set('electricconductivity', '5.998e7[S/m]');
model.material('mat3').propertyGroup('def').set('heatcapacity', '385[J/(kg*K)]');
model.material('mat3').propertyGroup('def').set('relpermittivity', '1');
model.material('mat3').propertyGroup('def').set('emissivity', '0.5');
model.material('mat3').propertyGroup('def').set('density', '8700[kg/m^3]');
model.material('mat3').propertyGroup('def').set('thermalconductivity', '400[W/(m*K)]');
model.material('mat3').propertyGroup.create('linzRes', 'Linearized resistivity');
model.material('mat3').propertyGroup('linzRes').set('alpha', '3.9e-3[1/K]');
model.material('mat3').propertyGroup('linzRes').set('rho0', '1.72e-8[ohm*m]');
model.material('mat3').propertyGroup('linzRes').set('Tref', '273.15[K]');
model.material('mat3').set('family', 'copper');
model.material('mat1').selection.named('Non_conducting');
model.material('mat2').selection.named('Winding');
model.material('mat3').selection.named('Conducting_Boundaries');

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Creates the physics
model.physics.create('mf', 'InductionCurrents', 'geom1');
model.physics('mf').feature.create('imp1', 'Impedance', 2);
model.physics('mf').feature('imp1').selection.named('Conducting_Boundaries');
model.physics('mf').feature.remove('imp1');
model.physics('mf').selection.set([]);
model.physics('mf').selection.named('Magnetic_Field');
model.physics('mf').feature.create('imp1', 'Impedance', 2);
model.physics('mf').feature('imp1').selection.named('Conducting_Boundaries');
model.physics('mf').feature.create('lport1', 'LumpedPort', 2);
model.physics('mf').feature('lport1').set('PortType', 1, 'UserDefined');
model.physics('mf').feature('lport1').set('ahPort', {'0' '0' '1'});
model.physics('mf').feature('lport1').set('TerminalType', 1, 'Current');
model.physics('mf').feature('lport1').set('hPort', 1, (h-t)/1000);%units in m
model.physics('mf').feature('lport1').set('wPort', 1, 0.004);%units in m
model.physics('mf').prop('MeshControl').set('EnableMeshControl', true);
model.physics('mf').feature('lport1').selection.named('Port1');
model.physics('mf').feature('lport1').set('I0', I_prim);

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Create mesh
model.mesh.create('mesh1', 'geom1');
model.mesh('mesh1').automatic(false);
model.mesh('mesh1').feature('size').set('custom', 'on');
model.mesh('mesh1').feature('size').set('hmin', '0.01');
model.mesh('mesh1').run;

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Create study
model.study.create('std1');
model.study('std1').feature.create('freq', 'Frequency');
model.study('std1').feature('freq').activate('mf', true);
model.study('std1').feature('freq').set('plist', freq);

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Run solution
model.sol.create('sol1');
model.sol('sol1').study('std1');
model.study('std1').feature('freq').set('notlistsolnum', 1);
model.study('std1').feature('freq').set('notsolnum', '1');
model.study('std1').feature('freq').set('listsolnum', 1);
model.study('std1').feature('freq').set('solnum', '1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'freq');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'freq');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('p1', 'Parametric');
model.sol('sol1').feature('s1').feature.remove('pDef');
model.sol('sol1').feature('s1').feature('p1').set('pname', {'freq'});
model.sol('sol1').feature('s1').feature('p1').set('plistarr', {'20000000'});
model.sol('sol1').feature('s1').feature('p1').set('punit', {'Hz'});
model.sol('sol1').feature('s1').feature('p1').set('pcontinuationmode', 'no');
model.sol('sol1').feature('s1').feature('p1').set('preusesol', 'auto');
model.sol('sol1').feature('s1').feature('p1').set('plot', 'off');
model.sol('sol1').feature('s1').feature('p1').set('plotgroup', 'Default');
model.sol('sol1').feature('s1').feature('p1').set('probesel', 'all');
model.sol('sol1').feature('s1').feature('p1').set('probes', {});
model.sol('sol1').feature('s1').feature('p1').set('control', 'freq');
model.sol('sol1').feature('s1').set('control', 'freq');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').create('i1', 'Iterative');
model.sol('sol1').feature('s1').feature('i1').set('linsolver', 'bicgstab');
model.sol('sol1').feature('s1').feature('i1').set('prefuntype', 'right');
model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').create('sv1', 'SORVector');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sv1').set('prefun', 'sorvec');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sv1').set('sorvecdof', {'mod1_A'});
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').create('sv1', 'SORVector');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sv1').set('prefun', 'soruvec');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sv1').set('sorvecdof', {'mod1_A'});
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'i1');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').attach('std1');

model.result.create('pg1', 'PlotGroup3D');
model.result('pg1').label('Magnetic Flux Density Norm (mf)');
model.result('pg1').set('data', 'dset1');
model.result('pg1').feature.create('mslc1', 'Multislice');
model.result('pg1').feature('mslc1').set('data', 'parent');

model.sol('sol1').runAll;

model.result('pg1').run;

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Get impedances
model.result.numerical.create('gev1', 'EvalGlobal');
model.result.numerical('gev1').setIndex('expr', 'mf.Zport_1', 1);
model.result.table.create('tbl1', 'Table');
model.result.numerical('gev1').set('table', 'tbl1');
model.result.numerical('gev1').setResult;
table = model.result.table('tbl1').getReal;
Z_real = table(2);
table = model.result.table('tbl1').getImag;
Z_imag = table(2);

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Save file
mphsave(model, strcat(mph_fileloc,mph_filename) )