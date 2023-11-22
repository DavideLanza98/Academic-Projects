function plot_Earth_2D()
% 
% plot_Earth_2D.m - plots the 2D representation of Earth's surface for
% Ground Track plotting and add settings
% 
% PROTOTYPE:
%  plot_Earth_2D()
%
% CONTRIBUTORS:
%  Lanza Davide
%  Larocca Rocco
%  Mascelloni Matteo 
%  Shakeel Afaq
% 
% VERSION:
%  07-01-2022

figure('Name', 'Ground Tracks');
hold on;
axis equal;
set(gca,'XTick',-180:30:180,'XTickMode','manual');
set(gca,'YTick',-90:30:90,'YTickMode','manual');
xlim([-180,180]); ylim([-90,90]);
      
image_file = 'earth.jpg';
cdata = flip(imread(image_file));
imagesc([-180,180],[-90, 90],cdata);

end

