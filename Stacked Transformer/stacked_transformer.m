%%  STACKED TOROIDAL TRANSFORMER AUTOMATION SCRIPT IN COMSOL
%
%   Author: Zikang Tong
%   Affiliation: Stanford Universiy
%   Description: This script determines the inductance and ac loss of an
%   interleaved toroidal transformer in COMSOL. The input variables include
%   the geometric parameters such as radius and height as well as the 
%   number of turns for each toroid. All units are in mm. 
%   
%   Step 1: Start with defining the toroid parameters 
%   Variable Definitions: 
%   ri = inner radius
%   ro = outer radius
%   h = height
%   t = copper thickness
%   separation = gap separation between adjacent windings

ri = 8;
ro = 15;
h = 4;
t = 1.5;
separation = 1;

%   Step 2: provide a list of number of turns 
%   Primary:
prim_list = [4, 4, 4, 6, 6];
%   Secondary:
sec_list = [7, 9, 9, 9, 9];
toroid_spacing = 1;

I_prim=1;
I_sec=0;
freq = 20e6;

mph_fileloc = 'G:\My Drive\Magnetics Project\MATLAB\';
mph_filename = 'delete1.mph';

%%%%%%%%%%%%%%%%%%%%%%%%%

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');
model.geom.create('geom1', 3);
model.geom('geom1').lengthUnit('mm');

prim_union_list = {};
for i = 1:length(prim_list)
    toroid_generator(model, ri,ro,h,t,separation,prim_list(i),strcat('prim',int2str(i)), i);
    model.geom('geom1').create(strcat('prim_mov',int2str(i)), 'Move');
    prim_union_list{end+1} = strcat('prim_mov',int2str(i));
    model.geom('geom1').feature(strcat('prim_mov',int2str(i))).selection('input').set(strcat('prim',int2str(i)));
    model.geom('geom1').feature(strcat('prim_mov',int2str(i))).set('displz', ((h+t+toroid_spacing)*2*(i-1)));
    model.geom('geom1').create(strcat('a_tab','prim',int2str(i)), 'Block');
    prim_union_list{end+1} = strcat('a_tab','prim',int2str(i));
    model.geom('geom1').feature(strcat('a_tab','prim',int2str(i))).set('size', [separation*4, separation*4, t]);
    model.geom('geom1').feature(strcat('a_tab','prim',int2str(i))).set('base', 'center');
    model.geom('geom1').feature(strcat('a_tab','prim',int2str(i))).set('pos', [ro+t/2+separation/2, 0, (-h/2)+((h+t+toroid_spacing)*2*(i-1))]);
    model.geom('geom1').create(strcat('b_tab','prim',int2str(i)), 'Block');
    prim_union_list{end+1} = strcat('b_tab','prim',int2str(i));
    model.geom('geom1').feature(strcat('b_tab','prim',int2str(i))).set('size', [separation*4, separation*4, t]);
    model.geom('geom1').feature(strcat('b_tab','prim',int2str(i))).set('base', 'center');
    model.geom('geom1').feature(strcat('b_tab','prim',int2str(i))).set('pos', [ro+t/2+separation/2, 0, (h/2)+((h+t+toroid_spacing)*2*(i-1))]);
    
    if i~=1
        model.geom('geom1').create(strcat('short','prim',int2str(i)), 'Block');
        prim_union_list{end+1} = strcat('short','prim',int2str(i));
        model.geom('geom1').feature(strcat('short','prim',int2str(i))).set('size', [separation*2, separation*2, 2*toroid_spacing+h+t]);
        model.geom('geom1').feature(strcat('short','prim',int2str(i))).set('base', 'center');
        model.geom('geom1').feature(strcat('short','prim',int2str(i))).set('pos', [ro+t/2+1.5*separation, 0, (t+h+toroid_spacing)*(2*(i-1)-1)]); 
    end
end
%   extra tabs
model.geom('geom1').create(strcat('bot_tab','prim',int2str(i)), 'Block');
prim_union_list{end+1} = strcat('bot_tab','prim',int2str(i));
model.geom('geom1').feature(strcat('bot_tab','prim',int2str(i))).set('size', [separation*4, separation*4, t]);
model.geom('geom1').feature(strcat('bot_tab','prim',int2str(i))).set('base', 'center');
model.geom('geom1').feature(strcat('bot_tab','prim',int2str(i))).set('pos', [ro+t/2+4.5*separation, 0, (-h/2)]);
model.geom('geom1').create(strcat('top_tab','prim',int2str(i)), 'Block');
prim_union_list{end+1} = strcat('top_tab','prim',int2str(i));
model.geom('geom1').feature(strcat('top_tab','prim',int2str(i))).set('size', [separation*4, separation*4, t]);
model.geom('geom1').feature(strcat('top_tab','prim',int2str(i))).set('base', 'center');
model.geom('geom1').feature(strcat('top_tab','prim',int2str(i))).set('pos', [ro+t/2+4.5*separation, 0, (h/2)+((h+t+toroid_spacing)*2*(length(prim_list)-1))]);
%   union all pieces
model.geom('geom1').create('uni_prim', 'Union');
model.geom('geom1').feature('uni_prim').selection('input').set(prim_union_list);
model.geom('geom1').feature('uni_prim').set('intbnd', 'off');
model.geom('geom1').feature('uni_prim').set('selresult', 'off');
model.geom('geom1').feature('uni_prim').set('keep', 'off');
model.geom('geom1').feature('uni_prim').set('selresult', 'on');
model.geom('geom1').feature('uni_prim').set('selresultshow', 'off');

%   port creation
model.geom('geom1').create('port1', 'Block');
model.geom('geom1').feature('port1').set('size', [separation, separation, 2*(h+t+toroid_spacing)*4+h-t]);
model.geom('geom1').feature('port1').set('base', 'corner');
model.geom('geom1').feature('port1').set('pos', [ro+t/2+5.5*separation, -separation/2, -h/2+t/2]);

model.geom('geom1').create('rot_prim', 'Rotate');
model.geom('geom1').feature('rot_prim').selection('input').set({'uni_prim', 'port1'});
model.geom('geom1').feature('rot_prim').set('rot', 30);


sec_union_list = {};
for i = 1:length(sec_list)
    toroid_generator(model, ri,ro,h,t,separation,sec_list(i),strcat('sec',int2str(i)), i);
    model.geom('geom1').create(strcat('sec_mov',int2str(i)), 'Move');
    sec_union_list{end+1} = strcat('sec_mov',int2str(i));
    model.geom('geom1').feature(strcat('sec_mov',int2str(i))).selection('input').set(strcat('sec',int2str(i)));
    model.geom('geom1').feature(strcat('sec_mov',int2str(i))).set('displz', ((h+t+toroid_spacing)*(2*i-1)));
    model.geom('geom1').create(strcat('a_tab','sec',int2str(i)), 'Block');
    sec_union_list{end+1} = strcat('a_tab','sec',int2str(i));
    model.geom('geom1').feature(strcat('a_tab','sec',int2str(i))).set('size', [separation*4, separation*4, t]);
    model.geom('geom1').feature(strcat('a_tab','sec',int2str(i))).set('base', 'center');
    model.geom('geom1').feature(strcat('a_tab','sec',int2str(i))).set('pos', [ro+t/2+separation/2, 0, (-h/2)+((h+t+toroid_spacing)*(2*i-1))]);
    model.geom('geom1').create(strcat('b_tab','sec',int2str(i)), 'Block');
    sec_union_list{end+1} = strcat('b_tab','sec',int2str(i));
    model.geom('geom1').feature(strcat('b_tab','sec',int2str(i))).set('size', [separation*4, separation*4, t]);
    model.geom('geom1').feature(strcat('b_tab','sec',int2str(i))).set('base', 'center');
    model.geom('geom1').feature(strcat('b_tab','sec',int2str(i))).set('pos', [ro+t/2+separation/2, 0, (h/2)+((h+t+toroid_spacing)*(2*i-1))]);
    if i~=1
        model.geom('geom1').create(strcat('short','sec',int2str(i)), 'Block');
        sec_union_list{end+1} = strcat('short','sec',int2str(i));
        model.geom('geom1').feature(strcat('short','sec',int2str(i))).set('size', [separation*2, separation*2, 2*toroid_spacing+h+t]);
        model.geom('geom1').feature(strcat('short','sec',int2str(i))).set('base', 'center');
        model.geom('geom1').feature(strcat('short','sec',int2str(i))).set('pos', [ro+t/2+1.5*separation, 0, (t+h+toroid_spacing)*(2*(i-1)-1)+h+t+toroid_spacing]); 
    end
end
%   extra tabs
model.geom('geom1').create(strcat('bot_tab','sec',int2str(i)), 'Block');
sec_union_list{end+1} = strcat('bot_tab','sec',int2str(i));
model.geom('geom1').feature(strcat('bot_tab','sec',int2str(i))).set('size', [separation*4, separation*4, t]);
model.geom('geom1').feature(strcat('bot_tab','sec',int2str(i))).set('base', 'center');
model.geom('geom1').feature(strcat('bot_tab','sec',int2str(i))).set('pos', [ro+t/2+4.5*separation, 0, (-h/2)+h+t+toroid_spacing]);
model.geom('geom1').create(strcat('top_tab','sec',int2str(i)), 'Block');
sec_union_list{end+1} = strcat('top_tab','sec',int2str(i));
model.geom('geom1').feature(strcat('top_tab','sec',int2str(i))).set('size', [separation*4, separation*4, t]);
model.geom('geom1').feature(strcat('top_tab','sec',int2str(i))).set('base', 'center');
model.geom('geom1').feature(strcat('top_tab','sec',int2str(i))).set('pos', [ro+t/2+4.5*separation, 0, (h/2)+((h+t+toroid_spacing)*2*(length(prim_list)-1))+h+t+toroid_spacing]);
%   union all pieces
model.geom('geom1').create('uni_sec', 'Union');
model.geom('geom1').feature('uni_sec').selection('input').set(sec_union_list);
model.geom('geom1').feature('uni_sec').set('intbnd', 'off');
model.geom('geom1').feature('uni_sec').set('selresult', 'off');
model.geom('geom1').feature('uni_sec').set('keep', 'off');
model.geom('geom1').feature('uni_sec').set('selresult', 'on');
model.geom('geom1').feature('uni_sec').set('selresultshow', 'off');

%   port2 creation
model.geom('geom1').create('port2', 'Block');
model.geom('geom1').feature('port2').set('size', [separation, separation, 2*(h+t+toroid_spacing)*4+h-t]);
model.geom('geom1').feature('port2').set('base', 'corner');
model.geom('geom1').feature('port2').set('pos', [ro+t/2+5.5*separation, -separation/2, h/2+3/2*t+toroid_spacing]);

model.geom('geom1').create('rot_sec', 'Rotate');
model.geom('geom1').feature('rot_sec').selection('input').set({'uni_sec','port2'});
model.geom('geom1').feature('rot_sec').set('rot', -30);

% create boundary
model.geom('geom1').create('sph1', 'Sphere');
model.geom('geom1').feature('sph1').set('r', max(length(prim_list), length(sec_list))*h*6);
model.geom('geom1').run;


%%%%%%%%%%%%%%%%%%%%%%%%%
%   Creates the selections
%   Creating geometric definitions
model.selection.create('Winding', 'Explicit');
model.selection('Winding').set([2,3]);
model.selection.create('Non_conducting', 'Explicit');
model.selection('Non_conducting').set([1,4,5]);
model.selection.create('Magnetic_Field', 'Explicit');
model.selection('Magnetic_Field').set([1,4,5]);
model.selection.create('Conducting_Boundaries', 'Explicit');
model.selection('Conducting_Boundaries').geom('geom1', 3, 2, {'exterior'});
model.selection('Conducting_Boundaries').set([2,3]);
model.selection.create('Port1', 'Explicit');
model.selection('Port1').geom('geom1', 3, 2, {'exterior'});
model.selection('Port1').set(5);
model.selection.create('Port2', 'Explicit');
model.selection('Port2').geom('geom1', 3, 2, {'exterior'});
model.selection('Port2').set(4);


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
model.physics('mf').feature('lport1').set('wPort', 1, 0.001*t/2*4);%units in m
model.physics('mf').feature.create('lport2', 'LumpedPort', 2);
model.physics('mf').feature('lport2').set('PortType', 1, 'UserDefined');
model.physics('mf').feature('lport2').set('ahPort', {'0' '0' '1'});
model.physics('mf').feature('lport2').set('TerminalType', 1, 'Current');
model.physics('mf').feature('lport2').set('hPort', 1, (h-t)/1000);%units in m
model.physics('mf').feature('lport2').set('wPort', 1, 0.001*t/2*4);%units in m
model.physics('mf').prop('MeshControl').set('EnableMeshControl', true);
model.physics('mf').feature('lport1').selection.named('Port1');
model.physics('mf').feature('lport2').selection.named('Port2');
model.physics('mf').feature('lport1').set('I0', I_prim);
model.physics('mf').feature('lport2').set('I0', I_sec);

flag = 'physics and geometry finished'

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Create mesh
model.mesh.create('mesh1', 'geom1');
model.mesh('mesh1').automatic(false);
model.mesh('mesh1').feature('size').set('custom', 'on');
model.mesh('mesh1').feature('size').set('hmin', t/100);
model.mesh('mesh1').run;

flag = 'mesh finished'

% %%%%%%%%%%%%%%%%%%%%%%%%%
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
model.result.numerical('gev1').setIndex('expr', 'mf.Zport_2', 1);
model.result.table.create('tbl1', 'Table');
model.result.numerical('gev1').set('table', 'tbl1');
model.result.numerical('gev1').setResult;
table = model.result.table('tbl1').getReal;
Z_prim_real = table(2);
Z_sec_real = table(3);
table = model.result.table('tbl1').getImag;
Z_prim_imag = table(2);
Z_sec_imag = table(3);

flag = 'simulation finished'

mphsave(model, strcat(mph_fileloc,mph_filename));

flag = 'file saved'
