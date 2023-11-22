function []=pl3d(pl,r_pl,figure,scale)
%
% pl3d.m - function useful to plot in 3d the planets included in the
%          interplanetary mission
%
% INPUTS:
%  pl      [1]    ID values for Mercury, Earth, Saturn and the Sun  [-]
%  r_pl    [3x1]  Position vector of the planet to plot             [km]
%  figure  [1]    Standard input, set f=figure in the script        [-]
%  scale   [1]    Desired size of the planet to plot                [-]
%
% OUTPUTS:
%  []      [-]    NaN                                               [-]
%
% First Edition: 01/07/2022
%
% Authors: 
% Rocco Larocca, Matteo Mascelloni, Davide Lanza, Afaq Shakeel

% just our planets
switch pl
    case 1
        pic='Mercury.jpg';        
        R_pl=astroConstants(21);        
    case 3
        pic='earth.jpg';
        R_pl=astroConstants(23);
    case 6
        pic='Saturn.png';
        R_pl=astroConstants(26);
    case 0
        pic='Sun.jpg';
        R_pl=astroConstants(3);  
end

axis_handle=gca(figure);             
[x,y,z]=sphere(100);                     

X=scale*R_pl*x+r_pl(1);
Y=scale*R_pl*y+r_pl(2);
Z=scale*R_pl*z+r_pl(3);

planet=surf(axis_handle,X,Y,Z,'FaceColor','texturemap','EdgeColor','none');        
A=imread(pic);                      
planet.CData=flip(A);