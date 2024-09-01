
lon1 = ncread('data path', 'lon');  % 251
lat1 = ncread('data path', 'lat');  % 501
lon = lon1(1:find(lon1==316));
lat = lat1(1:find(lat1==42.0));
d = ncread('data_H path', 'depth');  % 40
p = d(1:38);  % 0.0 2.0 4.0 6.0 8.0 10.0 12.0 15.0 20.0 25.0 30.0 35.0 40.0 45.0 50.0 60.0 70.0 80.0  
    % 90.0 100.0 125.0 150.0 200.0 250.0 300.0 350.0 400.0 500.0 600.0 700.0 800.0 900.0 1000.0
a = length(lon); 
b = length(lat); 
c = length(p); 
% Practical Salinity, SP; in-situ temperature, t; pressure, p
SP = ncread('data_H path', 'salinity');  % [251, 501, 40]
SP = SP(1:a, 1:b, 1:c);
SP = fillmissing(SP, 'linear');  % 选择的区域有陆地，不能进行填充
t = ncread('data_H path', 'water_temp');  % [251, 501, 40]
t = t(1:a, 1:b, 1:c);
t = fillmissing(t, 'linear');
    
SA = zeros(a, b, c);
CT = zeros(a, b, c);
for k = 1:c
    SA(:, :, k) = gsw_SA_from_SP(SP(:, :, k), p(k),lon,lat);
    CT(:, :, k) = gsw_CT_from_t(SA(:, :, k), t(:, :, k) ,p(k));
end
% compute potential density
pd = zeros(a, b, c);
for i = 1:c
    pd(:, :, i) = gsw_rho(SA(:, :, i),CT(:, :, i),p(i));
end
%% 获取密度异常数据+海水涩度数据
spiciness = spi(SA,CT);
% 读取 isQG 密度异常 和 HYCOM 密度异常
rho_i = ncread('data_i_time', 'rhot');
%%
t_i = zeros(a, b, c);
s_i = zeros(a, b, c);
SA_i = zeros(a, b, c);
CT_i = zeros(a, b, c);
for j = 1:c
    z = p(j);
    rho0 = nanmean(pd(:, :, j), 'all');
    rhoa_i = rho_i(:, :, j)+(rho0-1000)*ones(size(rho_i(:, :, j)));
    [SAi, CTi] = SA_CT(rhoa_i,spiciness(:, :, j), z);
    SA_i(:, :, j) = SAi;
    CT_i(:, :, j) = CTi;
    pr = z*ones(size(SAi));
    t_i(:, :, j) = gsw_t_from_CT(SAi,CTi,pr);
    s_i(:, :, j) = gsw_SP_from_SA(SAi,z,lon,lat);
end
t_i = fillmissing(t_i, 'linear');
s_i = fillmissing(s_i, 'linear');

%% caulate SA and CT
function [SA,CT] = SA_CT(sigma1,spiciness1,z)
%
% INPUT:
%  sigma1 =  density anomaly of a seawater sample (e.g 26 kg/m^3)
%            referenced to 500 dbar                             [ kg/m^3 ]
%   Note. This input has had 1000 kg/m^3 subtracted from it.  That is, 
%   it is 'density anomaly', not 'density'.
%
%  spiciness1 = spiciness appropriate to 500 dbar, as given by the paper
%               of McDougall and Krzysik (2015), and the GSW algorithm 
%               gsw_spiciness1(SA,CT)                            [ kg/m^3 ]
%  z = depth
%  sigma1 and spiciness1 must have the same dimensions.  
%
% OUTPUT:
%  SA  =  Absolute Salinity.                                       [ g/kg ]
%    Note, SA is expressed on the Reference-Composition 
%    Salinity Scale of Millero et al. (2008). 
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., and O.A. Krzysik, 2015: Spiciness. Journal of Marine 
%   Research, 73, 141-152.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling, 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org

if ~(nargin==3)
   error('gsw_SA_CT_from_sigma1_spiciness1:  Requires two inputs')
end 

[md,nd] = size(sigma1);
[msp,nsp] = size(spiciness1);

if (msp ~= md || nsp ~= nd)
    error('gsw_SA_CT_from_sigma1_spiciness1: sigma1 and spiciness1 must have same dimensions')
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------
rho = 1000.0 + sigma1; 

% initial estimate of SA from the polynomial of SA(sigma1,spiciness1)
SA = gsw_SA_poly_spiceness1(sigma1,spiciness1); 

[CT,~] = gsw_CT_from_rho(rho,SA,z); % The second solution is ignored.

%--------------------------------------------------------------------------
% Begin the modified Newton-Raphson iterative procedure of 
% McDougall and Wotherspoon (2014)
%--------------------------------------------------------------------------
for Number_of_iterations = 1:7
    delta_spiciness = spiciness1 - gsw_spiciness1(SA,CT);   
    derivative = gsw_deriv_SA_poly_spiciness1(sigma1,spiciness1);    
    SA_old = SA;
    SA = SA_old + delta_spiciness.*derivative;
    [CT, ~] = gsw_CT_from_rho(rho,SA,z); % The second solution is ignored. 
end

CT(SA < 0 | SA > 42) = NaN; 
SA(SA < 0 | SA > 42) = NaN;
CT(CT < -5 | CT > 40) = NaN; 
SA(CT < -5 |CT > 40) = NaN;

%--------------------------------------------------------------------------
% Note that this algorithm returns only one set of values of [SA,CT].
% At low salinities where the TMD is larger (warmer) than the freezing
% temperature there can be two solutions for the same input values 
% (sigma1,spiciness1).  
%--------------------------------------------------------------------------

end

function SA = gsw_SA_poly_spiceness1(sigma1,spiciness1)

 SApoly =[13.695625022104206
   0.627353321828843
   0.000074159905817
   0.711767497207624
  -0.000682109830188
   0.000726104526580
  -0.000091267411622
   0.000099817118989
   0.000022649012391
  -0.000039992686627];
 
SA = SApoly(1) + sigma1.*(SApoly(2) + SApoly(6).*spiciness1 + sigma1.*(SApoly(3) ...
    + SApoly(7).*spiciness1 + SApoly(9)*sigma1)) + spiciness1.*(SApoly(4) ...
    + spiciness1.*(SApoly(5) + SApoly(8)*sigma1 + SApoly(10)*spiciness1));

end

function SA_poly_spiciness1_derivative = gsw_deriv_SA_poly_spiciness1(sigma1,spiciness1)

 SApoly =[13.695625022104206
   0.627353321828843
   0.000074159905817
   0.711767497207624
  -0.000682109830188
   0.000726104526580
  -0.000091267411622
   0.000099817118989
   0.000022649012391
  -0.000039992686627];

 SA_poly_spiciness1_derivative = SApoly(4) + sigma1.*(SApoly(6) + SApoly(7)*sigma1) + ...
    spiciness1.*(2.*(SApoly(5) + SApoly(8)*sigma1) + 3.*SApoly(10)*spiciness1);

end

function spiciness = spi(SA,CT)
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
%  SA & CT need to have the same dimensions.
%
% OUTPUT:
%  spiciness1  =  spiciness referenced to a pressure of 1000 dbar 
%                                                                [ kg/m^3 ]
%
% AUTHOR: 
%  Oliver Krzysik and Trevor McDougall                 [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
%   Accurate and computationally efficient algorithms for potential 
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  McDougall, T.J., and O.A. Krzysik, 2015: Spiciness. Journal of Marine 
%   Research, 73, 141-152.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%   
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_spiciness1:  Requires two inputs')
end 

[ms,ns] = size(SA);
[mt,nt] = size(CT);

if (mt ~= ms | nt ~= ns)
    error('gsw_spiciness1: SA and CT must have same dimensions')
end

if ms == 1
    SA = SA';
    CT = CT';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

%deltaS = 24;
sfac = 0.0248826675584615;                   % sfac = 1/(40*(35.16504/35)).
offset = 5.971840214030754e-1;                      % offset = deltaS*sfac.

x2 = sfac.*SA;
xs = sqrt(x2 + offset);
ys = CT.*0.025;

s01 = -9.19874584868912e1;
s02 = -1.33517268529408e1;
s03 =  2.18352211648107e1;
s04 = -2.01491744114173e1;
s05 =  3.70004204355132e1;
s06 = -3.78831543226261e1;
s07 =  1.76337834294554e1;
s08 =  2.87838842773396e2;
s09 =  2.14531420554522e1;
s10 =  3.14679705198796e1;
s11 = -4.04398864750692e1;
s12 = -7.70796428950487e1;
s13 =  1.36783833820955e2;
s14 = -7.36834317044850e1;
s15 = -6.41753415180701e2;
s16 =  1.33701981685590;
s17 = -1.75289327948412e2;
s18 =  2.42666160657536e2;
s19 =  3.17062400799114e1;
s20 = -2.28131490440865e2;
s21 =  1.39564245068468e2;
s22 =  8.27747934506435e2;
s23 = -3.50901590694775e1;
s24 =  2.87473907262029e2;
s25 = -4.00227341144928e2;
s26 =  6.48307189919433e1;
s27 =  2.16433334701578e2;
s28 = -1.48273032774305e2;
s29 = -5.74545648799754e2;
s30 =  4.50446431127421e1;
s31 = -2.30714981343772e2;
s32 =  3.15958389253065e2;
s33 = -8.60635313930106e1;
s34 = -1.22978455069097e2;
s35 =  9.18287282626261e1;
s36 =  2.12120473062203e2;
s37 = -2.21528216973820e1;
s38 =  9.19013417923270e1;
s39 = -1.24400776026014e2;
s40 =  4.08512871163839e1;
s41 =  3.91127352213516e1;
s42 = -3.10508021853093e1;
s43 = -3.24790035899152e1;
s44 =  3.91029016556786;
s45 = -1.45362719385412e1;
s46 =  1.96136194246355e1;
s47 = -7.06035474689088;
s48 = -5.36884688614009;
s49 =  4.43247303092448;
 
spiciness = s01 + ys.*(s02 + ys.*(s03 + ys.*(s04 + ys.*(s05 + ys.*(s06 + s07*ys))))) ...
    + xs.*(s08 + ys.*(s09 + ys.*(s10 + ys.*(s11 + ys.*(s12 + ys.*(s13 + s14*ys)))))...
    + xs.*(s15 + ys.*(s16 + ys.*(s17 + ys.*(s18 + ys.*(s19 + ys.*(s20 + s21*ys))))) ...
    + xs.*(s22 + ys.*(s23 + ys.*(s24 + ys.*(s25 + ys.*(s26 + ys.*(s27 + s28*ys))))) ...
    + xs.*(s29 + ys.*(s30 + ys.*(s31 + ys.*(s32 + ys.*(s33 + ys.*(s34 + s35*ys))))) ...
    + xs.*(s36 + ys.*(s37 + ys.*(s38 + ys.*(s39 + ys.*(s40 + ys.*(s41 + s42*ys))))) ...
    + xs.*(s43 + ys.*(s44 + ys.*(s45 + ys.*(s46 + ys.*(s47 + ys.*(s48 + s49*ys)))))))))));

if transposed
    spiciness = spiciness.';
end
end